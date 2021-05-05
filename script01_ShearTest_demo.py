from yade import pack,plot,export
import math

sp         = pack.SpherePack()
O.periodic = True

# dimensions of sample (fixed by particle size such as L/D~X - length unit is mm)
DIAMETER  = 1.e-1
RADIUS    = 0.5*DIAMETER
length    = 15*(DIAMETER)*2
height    = length 
width     = length
thickness = length/100.

# microproperties
DENS     = 2600
E        = 1e9
P        = 0.25
compFRIC = 1.
FRIC     = 30.
TENS     = 0.
COH      = 0.

# boundary conditions
PI = 1.e5
conf = 0.
if conf<PI:
	PC = PI
else:
	PC = conf
SN = 5.e6 # normal stress for shearing phase
RATE = 0.02 # shearing rate

# simulation control
DAMPSHEAR = 0.
ITER      = 1e6
VTK       = 10
OUT       = 'l20h20w20_e1e9p025cf1f30_pi1e5pc0sn5e6_v0.02'

#### create sample and loading boxes

O.cell.hSize = Matrix3(length,0,0,0,3*height,0,0,0,width)

O.materials.append(CohFrictMat(isCohesive=True,density=DENS,young=E,poisson=P,frictionAngle=radians(0.),normalCohesion=1e100,shearCohesion=1e100,label='boxMat'))
O.materials.append(CohFrictMat(isCohesive=True,density=DENS,young=E,poisson=P,frictionAngle=radians(compFRIC),normalCohesion=TENS,shearCohesion=COH,label='sphereMat'))

upBox  = utils.box(center=(length/2.0,2*height+thickness/2.0,width/2.0),orientation=Quaternion(1,0,0,0),extents=(length,thickness/2.,width),fixed=1,wire=False,color=(1,0,0),material='boxMat') 
lowBox = utils.box(center=(length/2.0,height-thickness/2.0,width/2.0),orientation=Quaternion(1,0,0,0),extents=(length,thickness/2.,width),fixed=1,wire=False,color=(1,0,0),material='boxMat')
O.bodies.append([upBox,lowBox])

sp.makeCloud((0,height+1.5*RADIUS,0),(length,2*height-1.5*RADIUS,width),rMean=RADIUS,rRelFuzz=0.2,periodic=True)
O.bodies.append([utils.sphere(s[0],s[1],color=(0,0,1),material='sphereMat') for s in sp])

effCellVol = (O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2]
volRatio   = (O.cell.hSize[0,0]*O.cell.hSize[1,1]*O.cell.hSize[2,2])/effCellVol

#### engines

flow = PeriodicFlowEngine(
		isActivated         = 0, 
		useSolver          = 3,
		defTolerance       = -1, 
		meshUpdateInterval = 1000,
		duplicateThreshold = 0.5,
		boundaryUseMaxMin  = [0,0,1,1,0,0],
		wallIds            = [-1,-1,1,0,-1,-1],
		wallThickness      = thickness,
		bndCondIsPressure  = [0,0,0,0,0,0],
		bndCondValue       = [0,0,0,0,0,0],
		permeabilityFactor = 1,
		viscosity          = 1,
		fluidBulkModulus   = 2.2e9,
		label              = 'flowEng'
	)

O.engines = [
	ForceResetter()
	,InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb()],verletDist=-0.1,allowBiggerThanPeriod=True)
	,InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D(),Ig2_Box_Sphere_ScGeom6D()],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(traceEnergy=True,label='contactLaw')]
	)
	,flow
	,PeriTriaxController(
		dynCell       = True,
		mass          = 10,
		maxUnbalanced = 1e-3,
		relStressTol  = 1e-4,
		stressMask    = 7,
		goal          = (-PI/volRatio,-PI/volRatio,-PI/volRatio),
		globUpdate    = 1,
		maxStrainRate = (1,1,1),
		doneHook      = 'triaxDone()',
		label         = 'triax'
		)
	,GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8,defaultDt=utils.PWaveTimeStep(),label='timeStepper')
	,NewtonIntegrator(damping=0.3,label='newton')
	,PyRunner(command='dataRecorder()',iterPeriod=10,label='recData',dead=True)
	,PyRunner(iterPeriod=1,command='saveFlowVTK()',label='saveFluid',dead=True)
	,VTKRecorder(fileName=OUT+'.',iterPeriod=1,skipNondynamic=1,recorders=['spheres','colors','velocity','bstresses'],label='saveSolid',dead=True)
]

def saveFlowVTK():
	flow.saveVtk(folder='vtkFiles_'+OUT)

def dataRecorder():
	h             = vol=vol_s=nb_s=0.
	h             = O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
	vol           = h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
	contactStress = getStress(vol)
	for o in O.bodies:
			if isinstance(o.shape,Sphere) and o.shape.color[0]!=1:
				nb_s += 1
				vol_s += 4.*pi/3.*(o.shape.radius)**3
	n           = 1-vol_s/vol
	nbFrictCont = 0.
	for i in O.interactions:
		if i.isReal and i.phys.cohesionBroken: 
			nbFrictCont+=1
	plot.addData(
		iter             = O.iter,
		stress_upWall0  = abs(O.forces.f(0)[0]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])),
		stress_upWall1  = abs(O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])),
		stress_upWall2  = abs(O.forces.f(0)[2]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])),
		contactStress00 = (contactStress[0,0]),
		contactStress01 = (contactStress[0,1]),
		contactStress02 = (contactStress[0,2]),
		contactStress10 = (contactStress[1,0]),
		contactStress11 = (contactStress[1,1]),
		contactStress12 = (contactStress[1,2]),
		contactStress20 = (contactStress[2,0]),
		contactStress21 = (contactStress[2,1]),
		contactStress22 = (contactStress[2,2]),
		xW              = O.bodies[0].state.pos[0],
		height          = h,
		volume          = vol,
		porosity        = n,
		k               = 2.0*nbFrictCont/nb_s,
		Ek              = kineticEnergy(),
		unbF            = unbalancedForce()
		)

phase=0
def triaxDone():
	global phase
	volRatio=(O.cell.hSize[0,0]*O.cell.hSize[1,1]*O.cell.hSize[2,2])/((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])
	if phase==0:
		h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
		vol=h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
		contactStress=getStress(vol)
		vol_s=Rmean=Rmax=nbSph=0
		Rmin=1e6
		for o in O.bodies:
			if isinstance(o.shape,Sphere):
				nbSph+=1
				Rmean+=o.shape.radius
				if o.shape.radius>Rmax: Rmax=o.shape.radius
				if o.shape.radius<Rmin: Rmin=o.shape.radius
				vol_s += 4.*pi/3.*(o.shape.radius)**3
		Rmean=Rmean/nbSph
		n = 1-vol_s/vol
		print('DONE! iter =',O.iter,'| sample generated: nb spheres =',nbSph,', Rmean =',Rmean,', Rratio =',Rmax/Rmin,', porosity =',n)
		print('Changing contact properties now')
		utils.setContactFriction(radians(FRIC))
		print('APPLYING CONFINING PRESSURE: sx,sy and sz will go to PC =',PC)
		triax.goal=(-PC/volRatio,-PC/volRatio,-PC/volRatio)
		phase+=1
	elif phase==1:
		print('DONE! iter =',O.iter,'| isotropic confinement done: stresses =',volRatio*triax.stress)
		triax.dead=True
		O.pause()

#### Initialization
print('SAMPLE PREPARATION!')

#recData.dead=False # uncomment if you want to record what is happening during preparation of sample (isotropic compaction) 
O.run(1000000,1)
saveSolid.dead=False
O.step()
saveSolid.dead=True
O.save(OUT+'_isoConfined_'+str(O.iter)+'.yade.gz')

print('Normal stress (platen) =',O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
print('Normal stress (contacts) =',getStress((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])[1,1])

#### Applying normal stress
print('NORMAL LOADING! iter =',O.iter)

stage = 0
stiff = fnPlaten=currentSN=0.
def servo():
	global stage,stiff,fnPlaten,currentSN
	if stage==0:
		currentSN=O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])
		unbF=unbalancedForce()
		boundaryVel=copysign(min(0.1,abs(0.5*(currentSN-SN))),currentSN-SN)
		O.bodies[0].state.vel[1]=boundaryVel
		if ( (abs(currentSN-SN)/SN)<0.001 and unbF<0.001 ):
			stage+=1
			fnPlaten=O.forces.f(0)[1]
			print('Normal stress =',currentSN,' | unbF =',unbF)
			for i in O.interactions.withBody(O.bodies[0].id):
				stiff+=i.phys.kn
			print('Normal stiffness =',stiff)
			print('DONE! iter =',O.iter)
			O.pause()
	if stage==1:
		fnDesired = SN*(O.cell.hSize[0,0]*O.cell.hSize[2,2])
		boundaryVel = copysign(min(0.1,abs(0.35*(O.forces.f(0)[1]-fnDesired)/stiff/O.dt)),O.forces.f(0)[1]-fnDesired)
		O.bodies[0].state.vel[1]=boundaryVel

O.engines = O.engines[:4]+[PyRunner(command='servo()',iterPeriod=1,label='servo')]+O.engines[4:]

O.trackEnergy = True

O.run(1000000,1)
print('STABILIZING! iter=',O.iter)
O.run(1000,1)
dxi  = 4*(2*RADIUS)
n    = int(length/dxi)
xmin = 1e6
for i in range(0,n):
	for o in O.bodies:
		if o.id>1:
			if o.state.pos[0]<xmin: xmin=o.state.pos[0]
			if (o.state.pos[0]>=i*dxi) and (o.state.pos[0]<((i+0.5)*dxi)):
				o.shape.color[1]=1
markerSpheres=[]
for o in O.bodies:
	if (o.state.pos[0]>(xmin+0.5*O.cell.hSize[0,0]-1*RADIUS)) and (o.state.pos[0]<(xmin+0.5*O.cell.hSize[0,0]+2*RADIUS)):
		o.shape.color[2]=0
		markerSpheres.append(o)

def markerRecorder():
	global intOnFracPlane
	inFile=open(OUT+'_markerPosAndVel_'+str(O.iter),'a')
	for s in markerSpheres:
		inFile.write( str(s.state.pos[0]) + '\t' + str(s.state.pos[1]) + '\t' + str(s.state.pos[2]) + '\t' + str(s.state.vel[0]) + '\t' + str(s.state.vel[1]) + '\t' + str(s.state.vel[2]) + '\n')
	inFile.close()

O.engines = O.engines+[PyRunner(command='markerRecorder()',iterPeriod=1,label='recMark',dead=True)]

print('Normal stress (platen) = ',O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
print('Normal stress (contacts) = ',getStress((O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1])*O.cell.hSize[0,0]*O.cell.hSize[2,2])[1,1])

#### preparing for shearing
Gl1_Sphere.stripes=1

## gluing particles in contact with the walls
for i in O.interactions:
		if i.isReal:
			if isinstance(O.bodies[i.id1].shape,Box):
				O.bodies[i.id2].shape.color[0]=1
			if isinstance(O.bodies[i.id2].shape,Box):
				O.bodies[i.id1].shape.color[0]=1

for i in O.interactions:
		if i.isReal and ( O.bodies[i.id1].shape.color[0]==1 and O.bodies[i.id2].shape.color[0]==1 ):
			O.bodies[i.id1].mat.normalCohesion=O.bodies[i.id1].mat.normalCohesion
			O.bodies[i.id2].mat.normalCohesion=O.bodies[i.id1].mat.normalCohesion
			O.bodies[i.id1].mat.shearCohesion=O.bodies[i.id1].mat.shearCohesion
			O.bodies[i.id2].mat.shearCohesion=O.bodies[i.id1].mat.shearCohesion
			i.phys.initCohesion=True

saveSolid.dead=False
recMark.dead=False
O.step()
saveSolid.dead=True
recMark.dead=True
O.save(OUT+'_normallyLoaded_'+str(O.iter)+'.yade.gz')
if recData.dead==False: plot.saveDataTxt(OUT) 

#### shearing

recData.dead=False
recData.iterPeriod=int(ITER/10000.)
recMark.dead=False
recMark.iterPeriod=int(100000)
newton.damping=DAMPSHEAR
saveSolid.dead=False
saveSolid.iterPeriod=int(ITER/VTK)
shearVel=0
iterShear=recIter=O.iter
while 1:
	if ( O.iter >= int(recIter+100000) ):
		recIter=O.iter
		plot.saveDataTxt(OUT)
		O.save(OUT+'_sheared_'+str(O.iter)+'.yade.gz')
	if ( O.iter >= int(iterShear+ITER) ):
		print('iter =',O.iter,' -> END!')
		plot.saveDataTxt(OUT)
		O.save(OUT+'_sheared_'+str(O.iter)+'.yade.gz')
		O.save('laststate.yade.gz')
		sys.exit(0)
	O.run(int(10),1)
	if shearVel<RATE:
		shearVel+=(RATE/100.)
	O.bodies[0].state.vel[0]=shearVel
