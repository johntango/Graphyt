# Import the necessary modules

import sys
sys.path.insert(0,'/f/sjr/Graphyt/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
import graphyt
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
sys.path.insert(0,'/f/sjr/VTPWriter/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
import VTPWriter
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTIWriter/build/')
sys.path.insert(0,'/f/sjr/VTIWriter/build/')
import VTIWriter
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 1.0
psep = cellsize/2.0
# Background Grid
L = [35*cellsize, 20*cellsize, 4*cellsize]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(1e6)
water.setShearMod(1e5)
water.setDensity(10)
water.setModel(1)
materials = [water]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
liquidLowerLeft = [19*cellsize,7*cellsize,2*cellsize]
liquidUpperRight = [29*cellsize,13*cellsize,2*cellsize]
liquid = pyck.Cuboid(2,liquidLowerLeft,liquidUpperRight)
pack.AddShape(liquid)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) # Set the liquid to mat 1
# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
matpoints.addPyckIntField(matIDfield)

print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
for p in range(matpoints.numParticles):
    matpoints.setVel(p,1.0,0.0,0.0)

# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
   posY = nodes.getPosY(n)
   if ((posY> 15*cellsize ) and (posY < 17*cellsize )):
       nodes.setFixedBC(n)
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 0.1  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, 0, 0)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")

# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
#Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
print("+++++ SOLVER  INFORMATION+++++ Starting Solver")
solver = graphyt.Solver(nodes, matpoints, parameters, materials)

matpntPos = matpoints.getPosArr()
matpntMat = matpoints.getMatArr()
matpntMass = matpoints.getMassArr()
matpntVel = matpoints.getVelArr()
matpntStress = matpoints.getStressArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
gridPos = nodes.getPosArr()
gridMass = nodes.getMassArr()
gridVel = nodes.getVelArr()
gridFor = nodes.getForceArr()
vti = VTIWriter.VTIWriter(Le,cellsize)
vti.AddArray("Position", gridPos,3, 2)
vti.AddArray("Mass", gridMass,1, 2)
vti.AddArray("Velocity", gridVel,3, 2)
vti.AddArray("Forces", gridFor,3, 2)
# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/periodicTestmatpoints'+str(step)+'.vtp'
    vtifname='output/periodicTestnodes'+str(step)+'.vti'
    if step % 500 == 0:
        vtp.Write(vtpfname)
        vti.Write(vtifname)

print("Ratio: MatPoints/Nodes="+str(matpoints.numParticles)+"/"+str(nodes.numNodes)+"\nTotal Simulation Duration (s)")
# // ///////////////////////////////////////////////////







