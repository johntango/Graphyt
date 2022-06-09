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
sys.path.insert(0,'/f/sjr/VTIWriter/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTIWriter/build/')
import VTIWriter
import math
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.0115
psep = cellsize/2.0
# Background Grid
L = [1.60,0.50,0.750]
nodes = graphyt.Nodes(L, cellsize)
print("+++++ Grid INFO+++++ No. Nodes:", nodes.numNodes)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(2e9)
water.setShearMod(0.0)
water.setDensity(1000)
water.setModel(2)
materials = [water]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
liquidLowerLeft = [8.75*cellsize,8.75*cellsize,8.75*cellsize]
liquidUpperRight = [0.4,0.6,0.4]
liquid = pyck.Cuboid(1,liquidLowerLeft,liquidUpperRight)

# Add shapes to pack
pack.AddShape(liquid)
pack.Process()
model = pyck.Model(pack)
# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,0) 
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
objIDfield = model.GetIntField(objID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
matpoints.addPyckIntField(matIDfield)
matpoints.addPyckIntField(objIDfield)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////

# Set Hydrostatic pressure
#for p in range(matpoints.numParticles):
#    y = matpoints.getPosY(p)
#    pressure = 10*(0.75*L[1]-y)*1000
#    matpoints.setStress(p,pressure,pressure,pressure,0,0,0)

# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 8.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posX = nodes.getPosX(n)
    posY = nodes.getPosY(n)
    posZ = nodes.getPosZ(n)
    if (((posX>=L[0]-wallthickness) or (posX <= wallthickness)) or ((posY>=L[1]-wallthickness) or (posY <= wallthickness)) or ((posZ>=L[2]-wallthickness) or (posZ <= wallthickness))):
        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy_vz,[0,0,0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = True
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, -10, 0)
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
vti = VTIWriter.VTIWriter(Le,cellsize)
vti.AddArray("Position", gridPos,3, 2)
vti.AddArray("Mass", gridMass,1, 2)
vti.AddArray("Velocity", gridVel,3, 2)

# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.
    # enforce rotation
  #  matpoints.setRotatingBodyVel(0,wheelX,wheelY,0.0,omega)
    #vtifname='output/dambreak3D'+str(step)+'.vti'
    vtpfname='output/dambreak3D'+str(step)+'.vtp'
    if step % 200 == 0:
    #    vti.Write(vtifname)
       vtp.Write(vtpfname)


print("Total Simulation Duration (s)")
#// ///////////////////////////////////////////////////







