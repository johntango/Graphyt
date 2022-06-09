###########################################################

# USER INPUTS FOR GRAIN MODEL

# Loading
v = 0.1
tension = False

# Material Properties
# material 1
K1 = 3e9
G1 = 1e9
rho1 = 2000
# material 2
K2 =  K1
G2 = 0.5*G1
rho2 = rho1
###########################################################
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
import numpy
import math
## GEOMETRY ##
# Set the resolution and the pack
cellsize = 0.1
psep = cellsize/2.0
# Background Grid
L = [2,2,0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK
# create materials
mat1 = graphyt.MaterialModels()
mat1.setBulkMod(K1)
mat1.setShearMod(G1)
mat1.setDensity(rho1)
mat1.setModel(1)
mat2 = graphyt.MaterialModels()
mat2.setBulkMod(K2)
mat2.setShearMod(G2)
mat2.setDensity(rho2)
mat2.setModel(1)
materials = [mat1,mat2]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
# create shapes
topsolidLowerLeft = [0.5,1.0+1.25*cellsize,0.0]
topsolidUpperRight = [1.5-0.25*cellsize,1.5-0.25*cellsize,0.0]
topsolid = pyck.Cuboid(1,topsolidLowerLeft,topsolidUpperRight)
pack.AddShape(topsolid)
bottomsolidLowerLeft = [0.5,0.5,0.0]
bottomsolidUpperRight = [1.5-0.25*cellsize,1.0-1.25*cellsize,0.0]
bottomsolid = pyck.Cuboid(2,bottomsolidLowerLeft,bottomsolidUpperRight)
pack.AddShape(bottomsolid)
interfaceLowerLeft = [0.5,1.0-1.25*cellsize,0.0]
interfaceUpperRight = [1.5-0.25*cellsize,1.0+1.25*cellsize,0.0]
interface = pyck.Cuboid(3,interfaceLowerLeft,interfaceUpperRight)
pack.AddShape(interface)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
model.SetIntField(matID,2,0) 
model.SetIntField(matID,3,1) 

# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,1) 
model.SetIntField(objID,2,1) 
model.SetIntField(objID,3,1) 
# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
objIDfield = model.GetIntField(objID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
## Testing material by cell ID
# matpoints.addPyckIntField(matIDfield)
material_array = numpy.arange(0,numParticles,dtype=int)
for p in range(numParticles):
    # if p < 200:
    #     material_array[p]= 0
    # else:
        material_array[p] = 0
matpoints.setMaterialByCell(material_array)

matpoints.addPyckIntField(objIDfield)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
if tension:
    v = -v

for p in range(matpoints.numParticles):  
    y = matpoints.getPosY(p)
    if(y<0.6):
            matpoints.setBC(p,graphyt.BCTypes.mp_vy,[v])
    if(y>1.40):
            matpoints.setBC(p,graphyt.BCTypes.mp_vy,[-v])

# //// NODAL BOUNDARY CONDITIONS /////

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 100.0e-3  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, 0, 0)
parameters.setDampingCoef(0.0)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")

# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
#Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
print("+++++ SOLVER  INFORMATION+++++ Starting Solver")
solver = graphyt.Solver(nodes, matpoints, parameters, materials)

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
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Displacement",matpntDisp,3,2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("Normals",matpntNormal,3,2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
vtp.AddArray("Smax",matpntSMax,1,2)
# Iterate over nsteps number of timesteps
nStep = (tmax/3) / (solver.dt)
vStep = v / nStep
for step in range(int(solver.tmax / solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.
    if t > tmax / 3:
        v = v - vStep
        if v < 0.0:
            v = 0.0
    if t > (2 * tmax / 3):
        v = 0.0
        for p in range(matpoints.numParticles):
            y = matpoints.getPosY(p)
            if y < 0.6 or y > 1.40:
                matpoints.holdParticle(p)
 
    for p in range(matpoints.numParticles):
        y = matpoints.getPosY(p)
        if(y < 0.6):
            matpoints.setBC(p, graphyt.BCTypes.mp_vy,[v])
        if(y > 1.40):
            matpoints.setBC(p, graphyt.BCTypes.mp_vy,[-v])
    vtpfname = 'output/grainTest' + str(step) + '.vtp'
    # vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 1000 == 0:
        vtp.Write(vtpfname)
        # vti.Write(vtifname)
 
print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







