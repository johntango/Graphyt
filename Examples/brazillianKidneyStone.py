
# Import the necessary modules
import sys
import graphyt
import pyck
import VTKWriter
import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

# Other Info
# dt = 2.2e-7s
# Num Particles = 7944

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.0005
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [0.40,0.40,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
R = 0.027
R2 = R -  0.0015
R3 = R2 - 0.0015
R4 = R3 - 0.0015
R5 = R4 - 0.0015
R6 = R5 - 0.0015
R7 = R6 - 0.0015
R8 = R7 - 0.0015
R9 = R8 -  0.0015
R10 = R9 - 0.0015
R11 = R10 - 0.0015
R12 = R11 - 0.0015
R13 = R12 - 0.0015
R14 = R13 - 0.0015
R15 = R14 - 0.0015
solid1 = pyck.Sphere(3,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R)
# For stratified model
# solid2 = pyck.Sphere(4,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R2)
# solid3 = pyck.Sphere(5,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R3)
# solid4 = pyck.Sphere(6,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R4)
# solid5 = pyck.Sphere(7,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R5)
# solid6 = pyck.Sphere(8,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R6)
# solid7 = pyck.Sphere(9,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R7)
# solid8 = pyck.Sphere(10,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R8)
# solid9 = pyck.Sphere(11,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R9)
# solid10 = pyck.Sphere(12,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R10)
# solid11 = pyck.Sphere(13,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R11)
# solid12 = pyck.Sphere(14,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R12)
# solid13 = pyck.Sphere(15,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R13)
# solid14 = pyck.Sphere(16,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R14)
# solid15 = pyck.Sphere(17,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R15)
# For marbled model
R16 = R15 - 0.003
solid2 = pyck.Sphere(4,[40e-3+R/2.0 + 0.0075, 40e-3+R/2.0 + 0.005,0.0],R16)
solid3 = pyck.Sphere(5,[40e-3+R/2.0 + 0.0012, 40e-3+R/2.0 - 0.001,0.0],R16+0.001)
solid4 = pyck.Sphere(6,[40e-3+R/2.0 - 0.005, 40e-3+R/2.0 - 0.001,0.0],R16-0.001)
solid5 = pyck.Sphere(7,[40e-3+R/2.0, 40e-3+R/2.0 + 0.01,0.0],R16)
solid6 = pyck.Sphere(8,[40e-3+R/2.0 - 0.007, 40e-3+R/2.0 + 0.005,0.0],R16+0.0005)
solid7 = pyck.Sphere(9,[40e-3+R/2.0 + 0.013, 40e-3+R/2.0 - 0.003,0.0],R16-0.0005)
solid8 = pyck.Sphere(10,[40e-3+R/2.0 - 0.001, 40e-3+R/2.0 - 0.01,0.0],R16+0.0002)
solid9 = pyck.Sphere(11,[40e-3+R/2.0 - 0.013, 40e-3+R/2.0 - 0.003,0.0],R16-0.0005)

topPlaten_lowerLeft = [0.04+R/2.0-3*cellsize,0.04+R/2.0+R,0.0]
topPlaten_upperRight = [0.04+R/2.0+3*cellsize,0.04+R/2.0+R+0.0025,0.0] 
topPlaten = pyck.Cuboid(1,topPlaten_lowerLeft,topPlaten_upperRight)
bottomPlaten_lowerLeft = [0.04+R/2.0-3*cellsize,0.04+R/2.0-R,0.0]
bottomPlaten_upperRight = [0.04+R/2.0+3*cellsize,0.04+R/2.0-R-0.0025,0.0] 
bottomPlaten = pyck.Cuboid(2,bottomPlaten_lowerLeft,bottomPlaten_upperRight)
### PYCK - Add all shapes to pack
pack.AddShape(solid1)
pack.AddShape(solid2)
pack.AddShape(solid3)
pack.AddShape(solid4)
pack.AddShape(solid5)
pack.AddShape(solid6)
pack.AddShape(solid7)
pack.AddShape(solid8)
pack.AddShape(solid9)
# pack.AddShape(solid10)
# pack.AddShape(solid11)
# pack.AddShape(solid12)
# pack.AddShape(solid13)
# pack.AddShape(solid14)
# pack.AddShape(solid15)
pack.AddShape(topPlaten)
pack.AddShape(bottomPlaten)
pack.Process()
### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)
print("+++++ GRID INFO+++++ No. Nodes:", nodes.numNodes)

################################################################
#--------------------------------------------------------------#
#                        Materials                             #
#--------------------------------------------------------------#
################################################################

### PYCK - Materials
rock1 = graphyt.MaterialModels()
rock1.setBulkMod(3.31e9)
rock1.setShearMod(2.48e9)
rock1.setDensity(1540.0)
rock1.setCohesion(5.20e6)
rock1.setCriticalStrain(0.1)
rock1.setFrictionAngle(55.81)
rock1.setDilatancyAngle(55.81)#(13.95)
rock1.setIsDamage(True)
rock1.setDamageModel(1)
rock1.setCrackSpeed(586.52)
rock1.setDamageM(9)
rock1.setDamageK(7.0e35)#(5.39e36)
rock1.setModel(4)
rock2 = graphyt.MaterialModels()
rock2.setBulkMod(0.7*3.31e9)
rock2.setShearMod(0.7*2.48e9)
rock2.setDensity(1540.0)
rock2.setCohesion(5.20e6)
rock2.setCriticalStrain(0.1)
rock2.setFrictionAngle(55.81)
rock2.setDilatancyAngle(55.81)#(13.95)
rock2.setIsDamage(True)
rock2.setDamageModel(1)
rock2.setCrackSpeed(586.52)
rock2.setDamageM(9)
rock2.setDamageK(7.0e35)#(5.39e36)
rock2.setModel(4)

platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [platen,rock1, rock2]
### PYCK - Create Model for Simulation
model = pyck.Model(pack)
matID = model.CreateIntField("material",1)
objID = model.CreateIntField("objectID",1)
# Assign material IDs to each object
model.SetIntField(matID,1,0)
model.SetIntField(matID,2,0)
model.SetIntField(matID,3,1)
model.SetIntField(matID,4,2)
model.SetIntField(matID,5,2)
model.SetIntField(matID,6,2)
model.SetIntField(matID,7,2)
model.SetIntField(matID,8,2)
model.SetIntField(matID,9,2)
model.SetIntField(matID,10,2)
model.SetIntField(matID,11,2)
# model.SetIntField(matID,12,2)
# model.SetIntField(matID,13,1)
# model.SetIntField(matID,14,2)
# model.SetIntField(matID,15,1)
# model.SetIntField(matID,16,2)
# model.SetIntField(matID,17,1)
matIDfield = model.GetIntField(matID)
matpoints.addPyckIntField(matIDfield)
# Assign object IDs to each object
model.SetIntField(objID,1,1)
model.SetIntField(objID,2,2)
model.SetIntField(objID,3,3)
model.SetIntField(objID,4,3)
model.SetIntField(objID,5,3)
model.SetIntField(objID,6,3)
model.SetIntField(objID,7,3)
model.SetIntField(objID,8,3)
model.SetIntField(objID,9,3)
model.SetIntField(objID,10,3)
model.SetIntField(objID,11,3)
# model.SetIntField(objID,12,3)
# model.SetIntField(objID,13,3)
# model.SetIntField(objID,14,3)
# model.SetIntField(objID,15,3)
# model.SetIntField(objID,16,3)
# model.SetIntField(objID,17,3)
objIDField = model.GetIntField(objID)
matpoints.addPyckIntField(objIDField)
################################################################
#--------------------------------------------------------------#
#                     Boundary Conditions                      #
#--------------------------------------------------------------#
################################################################

### Material Point BCs
for p in range(numParticles):
    obID = matpoints.getObjID(p)
    if(obID == 1):
        matpoints.setVel(p,0.0,-0.2,0)
    if(obID == 2):
        matpoints.setVel(p,0.0,0.2,0)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 8e-3
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
# parameters.setFriction(0.3)
parameters.setCFL(0.1)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")
################################################################
#--------------------------------------------------------------#
#                   Simulation Solver                          #
#--------------------------------------------------------------#
################################################################

## Create the solver object
solver = graphyt.Solver(nodes, matpoints, parameters, materials)
## Define the output arrays
matpntPos = matpoints.getPosArr()
matpntDisp = matpoints.getDispArr()
matpntMat = matpoints.getMatArr()
matpntMass = matpoints.getMassArr()
matpntVel = matpoints.getVelArr()
matpntStress = matpoints.getStressArr()
matpntNormal = matpoints.getNormalsArr()
matpntAccPlasStr = matpoints.getAccPlasticStrainArr()
matpntDamage = matpoints.getDamageArr()
matpntSMax = matpoints.getSMaxArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
## Create a VTP writer and add the arrays
vtp = VTKWriter.VTPWriter(numParticles)
vtp.AddPositions("Position", matpntPos,3, VTKWriter.Order.ijk)
vtp.AddArray("Displacement",matpntDisp,3,VTKWriter.Order.ijk)
vtp.AddArray("Material", matpntMat,1, VTKWriter.Order.ijk)
vtp.AddArray("Mass", matpntMass,1, VTKWriter.Order.ijk)
vtp.AddArray("Velocity",matpntVel,3, VTKWriter.Order.ijk)
vtp.AddArray("Stress",matpntStress,6, VTKWriter.Order.ijk)
vtp.AddArray("Normals",matpntNormal,3,VTKWriter.Order.ijk)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,VTKWriter.Order.ijk)
vtp.AddArray("Damage",matpntDamage,1,VTKWriter.Order.ijk)
vtp.AddArray("Smax",matpntSMax,1,VTKWriter.Order.ijk)
print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    # vtpfname='output/kidneyModelstrat'+str(step)+'.vtp'
    vtpfname='output/kidneyModelmarbled'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 1000 == 0):
        vtp.Write(vtpfname)
print("Total Simulation Duration (s)")
