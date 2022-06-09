
# Import the necessary modules
import sys
sys.path.insert(0,'/f/sjr/Graphyt/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
sys.path.insert(0,'C:/Users/sjr/Desktop/PhD/Graphyt/build-win/Release/')
import graphyt
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
sys.path.insert(0,'C:/Users/sjr/Desktop/PhD\pyck/build-win/swig/Release')
sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
sys.path.insert(0,'/f/sjr/VTPWriter/build/')
sys.path.insert(0,'C:/Users/sjr/Desktop/PhD/VTPWriter/build-win/Release')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
import VTPWriter
import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.05
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [10.0,1.0,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
block = pyck.Cuboid(1,[0.3,0.3,0.0],[0.5,0.4,0.0])
# disk = pyck.Sphere(1,[0.4,0.35,0.0],0.05)
ramp = pyck.Cuboid(2,[0.02,0.3-4*cellsize,0.0],[9.98,0.3,0.0])
### PYCK - Add all shapes to pack
pack.AddShape(block)
# pack.AddShape(disk)
pack.AddShape(ramp)
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
blockMat = graphyt.MaterialModels()
blockMat.setBulkMod(1e6)
blockMat.setShearMod(7e5)
blockMat.setDensity(1000.0)
blockMat.setModel(1)

rampMat = graphyt.MaterialModels()
rampMat.setBulkMod(1e6)
rampMat.setShearMod(7e5)
rampMat.setDensity(1000.0)
rampMat.setModel(1)
# rampMat.setIsRigid(True)
materials = [blockMat, rampMat]
### PYCK - Create Model for Simulation
model = pyck.Model(pack)
matID = model.CreateIntField("material",1)
objID = model.CreateIntField("objectID",1)
# Assign material IDs to each object
model.SetIntField(matID,1,0)
model.SetIntField(matID,2,0)
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
# for p in range(matpoints.numParticles):
#     objID = matpoints.getObjID(p)
#     if(objID==2):
#         matpoints.setVel(p,0,0,0)
## Nodal BCs
for n in range(nodes.numNodes):
    yPos = nodes.getPosY(n)
    if(yPos <0.3 - 2*cellsize):
        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0.0,0.0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 20
is3D = False
mu = 30
parameters = graphyt.Parameters(tmax, is3D)
parameters.setGravity(0.4,-0.6,0)
parameters.setFriction(mu/100.0)
parameters.setCFL(0.001)
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
matpntAcc = matpoints.getAccArr()
matpntStress = matpoints.getStressArr()
matpntNormal = matpoints.getNormalsArr()
matpntAccPlasStr = matpoints.getAccPlasticStrainArr()
matpntDamage = matpoints.getDamageArr()
matpntSMax = matpoints.getSMaxArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
## Create a VTP writer and add the arrays
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Displacement",matpntDisp,3,2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Acceleration",matpntAcc,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("Normals",matpntNormal,3,2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
vtp.AddArray("Smax",matpntSMax,1,2)
print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    # vtpfname='output/inclinedPlane_rolling'+str(step)+'.vtp'
    # vtpfname='output/friction/coef01/inclinedPlanecoef01'+str(step)+'.vtp'
    # vtpfname='output/friction/coef03/inclinedPlanecoef03'+str(step)+'.vtp'
    vtpfname='output/inclinedPlanecoef_'+str(mu)+'_'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 10000 == 0):
        vtp.Write(vtpfname)
print("Total Simulation Duration (s)")
