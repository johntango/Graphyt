# Import the necessary modules
#import sys
#sys.path.insert(0,'/f/sjr/Graphyt/build/')
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
import graphyt
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
#sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
# sys.path.insert(0,'/f/sjr/VTPWriter/build/')
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/VTKWriter/build')
import VTKWriter

import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.0001
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [0.10,0.10,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry

cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0.0])
pack = pyck.StructuredPack(cubic)

# square1_lowerLeft = [0.045,0.02,0.0]
# square1_upperRight = [0.055,0.05,0.0] 
# square1 = pyck.Cuboid(1,square1_lowerLeft,square1_upperRight)
disk1 = pyck.Sphere(1,[0.050,0.03,0.0],0.005)
# square2_lowerLeft = [0.045,square1_upperRight[1]+3*cellsize,0.0]
# square2_upperRight = [0.055,square2_lowerLeft[1]+0.03,0.0] 
# square2 = pyck.Cuboid(2,square2_lowerLeft,square2_upperRight)
disk2 = pyck.Sphere(2,[0.050,0.03+0.01 + 3*cellsize,0.0],0.005)
### PYCK - Add all shapes to pack
pack.AddShape(disk1)
pack.AddShape(disk2)
pack.Process()
### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)
################################################################
#--------------------------------------------------------------#
#                        Materials                             #
#--------------------------------------------------------------#
################################################################

### PYCK - Materials
mat1 = graphyt.MaterialModels()
mat1.setBulkMod(1000)
mat1.setShearMod(500)
mat1.setDensity(10)
mat2 = graphyt.MaterialModels()
mat2.setBulkMod(1000)
mat2.setShearMod(500)
mat2.setDensity(10)
materials = [mat1, mat2]
### PYCK - Create Model for Simulation
model = pyck.Model(pack)
matID = model.CreateIntField("material",1)
objID = model.CreateIntField("objectID",1)
# Assign material IDs to each object
model.SetIntField(matID,1,0)
model.SetIntField(matID,2,1)
matIDfield = model.GetIntField(matID)
matpoints.addPyckIntField(matIDfield)
# Assign object IDs to each object
model.SetIntField(objID,1,1)
model.SetIntField(objID,2,2)
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
        matpoints.setVel(p,0.0,1.50,0)
    if(obID == 2):
        matpoints.setVel(p,0.0,-1.50,0)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 5e-3
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
# parameters.setDampingCoef(1.0)
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

gridMass = nodes.getMassArr()
gridForce = nodes.getForceArr()
gridVel = nodes.getVelArr()
vti = VTKWriter.VTIWriter(Le,[-0.00025,-0.00025,-0.00025],cellsize)
vti.AddArray("Mass",gridMass,1, VTKWriter.Order.ijk)
vti.AddArray("Force",gridForce,3, VTKWriter.Order.ijk)
vti.AddArray("Velocity",gridVel,3,VTKWriter.Order.ijk)
print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    vtpfname='output/disks2D'+str(step)+'.vtp'
    vtifname='output/disks2Dvti'+str(step)+'.vti'
    if (step % 200 == 0):
        vtp.Write(vtpfname)
        vti.Write(vtifname)

print(matpoints.getDispY(int(round(matpoints.numParticles/4))))
print(matpoints.getDispY(int(round(3*matpoints.numParticles/4))))
print(matpoints.getStressYY(int(round(matpoints.numParticles/4))))
print(matpoints.getStressYY(int(round(3*matpoints.numParticles/4))))
print("Total Simulation Duration (s)")
