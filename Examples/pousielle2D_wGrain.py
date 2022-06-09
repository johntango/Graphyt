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
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTIWriter/build/')
# sys.path.insert(0,'/f/sjr/VTIWriter/build/')
# import VTIWriter
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.01
psep = cellsize/2.0
# Background Grid
L = [60*cellsize, 40*cellsize, 0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(2e9)
water.setShearMod(0.0)
water.setViscosity(1.0e-3)
water.setDensity(1000)
water.setModel(2)
grain = graphyt.MaterialModels()
# grain.setBulkMod(3e9)
# grain.setShearMod(1e9)
# grain.setDensity(1500)
grain.setModel(1)
grain.setIsRigid(True)
materials = [water,grain]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
# create shapes
liquidLowerLeft = [2*cellsize,L[1]/4,0.0]
liquidUpperRight = [L[0]-2*cellsize,3*L[1]/4,0.0]
liquid = pyck.Cuboid(1,liquidLowerLeft,liquidUpperRight)
R = 3*cellsize
grain = pyck.Sphere(2,[L[0]/2,L[1]/2,0.0],R)
void = pyck.Sphere(0,[L[0]/2,L[1]/2,0.0],R+2*cellsize)
pack.AddShape(liquid)
pack.AddShape(void)
pack.AddShape(grain)
pack.Process()
### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
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

print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
for p in range(matpoints.numParticles):
    objId = matpoints.getObjID(p)
    if(objId == 2):
        matpoints.setVel(p,0.0,0,0)


# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posY = nodes.getPosY(n)
    if ((posY>L[1]-13*cellsize ) or (posY < 12*cellsize )):
        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(-10, 0, 0)
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
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/pouiselle_grainVTP'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 2000 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







