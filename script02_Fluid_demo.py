from yade import pack,plot,export

prevSim= 'laststate'
O.load(prevSim+'.yade.gz')

sigmaN = 5.e6
rTau   = 0.8
output = prevSim+'_output'

## time stepper (needed if time control required)
timer        = utils.typedEngine('GlobalStiffnessTimeStepper')
timer.active = False
timer.dead   = True
O.dt         = 0.5*utils.PWaveTimeStep()

## vtk recorder
saveSolid.dead           = True
saveSolid.iterPeriod     = int(1)
saveSolid.fileName       = output+'.'
saveSolid.skipNondynamic = 1
saveSolid.recorders      = ['spheres','colors','velocity','bstresses']

def saveFlowVTK():
 flowEng.saveVtk(folder=output+'_vtkFluid')
saveFluid.dead       = True
saveFluid.iterPeriod = int(1)

## data recorder
inputP=0.
def dataRecorder():
	global inputP
	h=vol=vol_s=nb_s=0.
	h=O.bodies[0].state.pos[1]-O.bodies[1].state.pos[1]
	vol=h*O.cell.hSize[0,0]*O.cell.hSize[2,2]
	contactStress=getStress(vol)
	for o in O.bodies:
			if isinstance(o.shape,Sphere) and o.shape.color[0]!=1:
				nb_s += 1
				vol_s += 4.*pi/3.*(o.shape.radius)**3
	n = 1-vol_s/vol
	nbFrictCont=0.
	for i in O.interactions:
		if i.isReal and i.phys.cohesionBroken: 
			nbFrictCont+=1
	plot.addData(
		iter             = O.iter
		,stress_upWall0  = abs(O.forces.f(0)[0]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
		,stress_upWall1  = abs(O.forces.f(0)[1]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
		,stress_upWall2  = abs(O.forces.f(0)[2]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
		,contactStress00 = (contactStress[0,0])
		,contactStress01 = (contactStress[0,1])
		,contactStress02 = (contactStress[0,2])
		,contactStress10 = (contactStress[1,0])
		,contactStress11 = (contactStress[1,1])
		,contactStress12 = (contactStress[1,2])
		,contactStress20 = (contactStress[2,0])
		,contactStress21 = (contactStress[2,1])
		,contactStress22 = (contactStress[2,2])
		,xW              = O.bodies[0].state.pos[0]
		,height          = h
		,volume          = vol
		,porosity        = n
		,k               = 2.0*nbFrictCont/nb_s
		,Ek              = kineticEnergy()
		,unbF            = unbalancedForce()
	)

## positions recorder
markerSpheres=[]
for o in O.bodies:
	if o.shape.color[2]==0:
		markerSpheres.append(o)

print('nb or markers=',len(markerSpheres))

def markerRecorder():
	global intOnFracPlane
	inFile=open(output+'_markerPosAndVel_'+str(O.iter),'a')
	for s in markerSpheres:
		inFile.write( str(s.state.pos[0]) + '\t' + str(s.state.pos[1]) + '\t' + str(s.state.pos[2]) + '\t' + str(s.state.vel[0]) + '\t' + str(s.state.vel[1]) + '\t' + str(s.state.vel[2]) + '\n')
	inFile.close()
	
recMark.dead=True
recMark.iterPeriod=1

#### servocontrol of top plate:

# we need one step to have a force on platen
for o in O.bodies:
	o.dynamic=False
saveSolid.dead=False
recMark.dead=False
O.run(1,1)
saveSolid.dead=True
recMark.dead=True
O.save(output+'_'+str(O.iter)+'.yade.gz')
for o in O.bodies:
	if isinstance(o.shape,Sphere): o.dynamic=True

nStiff=0
for i in O.interactions.withBody(O.bodies[0].id):
	nStiff+=i.phys.kn
print('normal stiffness=',nStiff)
initFnPlaten=O.forces.f(0)[1]
initSigmaN=initFnPlaten/(O.cell.hSize[0,0]*O.cell.hSize[2,2])
print('normal stress (platen) =',initSigmaN)

def servo():
	fnDesired=sigmaN*(O.cell.hSize[0,0]*O.cell.hSize[2,2])
	nBoundaryVel=copysign(min(0.1,abs(0.35*(O.forces.f(0)[1]-fnDesired)/nStiff/O.dt)),O.forces.f(0)[1]-fnDesired)
	#nBoundaryVel=copysign(abs(0.35*(O.forces.f(0)[1]-fnDesired)/nStiff/O.dt),O.forces.f(0)[1]-fnDesired)
	O.bodies[0].state.vel[1]=nBoundaryVel

#### Stabilization
print('stabilizing | iter=',O.iter)
newton.damping=0.3
O.bodies[0].state.vel[0]=0
recData.dead=False
recData.iterPeriod=1

O.run(int(1e3),1)
saveSolid.dead=False
recMark.dead=False
O.run(1,1)
saveSolid.dead=True
recMark.dead=True
plot.saveDataTxt(output)
O.save(output+'_'+str(O.iter)+'.yade.gz')

#### shear stress control

sStiff=0
for i in O.interactions.withBody(O.bodies[0].id):
	sStiff+=i.phys.ks
print('shear stiffness=',sStiff)
initFsPlaten=O.forces.f(0)[0]
initTau=initFsPlaten/(O.cell.hSize[0,0]*O.cell.hSize[2,2])
print('shear stress (platen) =',initTau)

maxSBVel=0.1
def servoShear():
	global maxSBVel
	fsDesired=rTau*initTau*(O.cell.hSize[0,0]*O.cell.hSize[2,2])
	sBoundaryVel=copysign(min(maxSBVel,abs(0.35*(O.forces.f(0)[0]-fsDesired)/sStiff/O.dt)),O.forces.f(0)[0]-fsDesired)
	O.bodies[0].state.vel[0]=sBoundaryVel

O.engines = O.engines[:5]+[PyRunner(command='servoShear()',iterPeriod=1,label='servoShear')]+O.engines[5:]

print('shear control now | iter=',O.iter)

while 1:
	O.run(100,1)
	if (unbalancedForce()<0.001) and ( ((abs(O.forces.f(0)[0]/(O.cell.hSize[0,0]*O.cell.hSize[2,2])) - abs(rTau*initTau))/abs(rTau*initTau))<0.001 ):
		print ('stress state reached| iter=',O.iter,' | shear stress (platen) =',O.forces.f(0)[0]/(O.cell.hSize[0,0]*O.cell.hSize[2,2]))
		print('stabilizing')
		O.run(int(1e3),1)
		plot.saveDataTxt(output)
		saveSolid.dead=False
		recMark.dead=False
		O.run(1,1)
		saveSolid.dead=True
		recMark.dead=True
		O.save(output+'_'+str(O.iter)+'.yade.gz')
		break

#### injection
print('FLUID NOW | iter=',O.iter)
maxSBVel                      = 1e3
flowEng.isActivated           = 1
flowEng.bndCondIsPressure     = [0,0,1,0,0,0]
deltaP                        = 1e5
inputP                        = deltaP
flowEng.bndCondValue          = [0.,0.,inputP,0.,0.,0.]
flowEng.fluidBulkModulus      = 2.2e9
flowEng.permeabilityFactor    = 1
flowEng.viscosity             = 1
newton.damping                = 0.

iterInit      = injectIter=recIter=O.iter
deltaIter     = 3000
iterMax       = 20*deltaIter
runningIter   = O.iter

while 1 :
	O.run(1,1)
	if ( O.iter >= int(injectIter+deltaIter) ):
		injectIter = O.iter
		inputP     +=deltaP
		flowEng.bndCondValue=[0.,0.,inputP,0.,0.,0.]
		flowEng.updateBCs()
		print('updateBCs! inputP=',inputP)
	if ( O.iter >= int(recIter+int(deltaIter/2.)) ):
		recIter        = O.iter
		saveSolid.dead = False
		saveFluid.dead = False
		recMark.dead   = False
		O.run(1,1)
		saveSolid.dead = True
		saveFluid.dead = True
		recMark.dead   = True
		plot.saveDataTxt(output)
		O.save(output+'_'+str(O.iter)+'.yade.gz')
		print('saving data!')
	if ( O.iter >= int(iterInit+iterMax) ):
		print('iter=',O.iter,' -> END!')
		plot.saveDataTxt(output)
		O.save(output+'_'+str(O.iter)+'.yade.gz')
		break