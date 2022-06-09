# Import the necessary modules
import sys
import graphyt
import pyck
import VTKWriter
import math
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.02
psep = cellsize/2.0
# Background Grid
L = [1.0,1.0,0.0]
nodes = graphyt.Nodes(L, cellsize)
print("+++++ Grid INFO+++++ No. Nodes:", nodes.numNodes)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(1e6)
water.setShearMod(0.0)
water.setDensity(1000)
water.setModel(2)
materials = [water]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0.0])
pack = pyck.StructuredPack(cubic)
# create shapes
liquidLowerLeft = [2.75*cellsize,2.75*cellsize,0.0]
liquidUpperRight = [0.5*L[0],0.75*L[1]-2*cellsize,0.0]
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


# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posX = nodes.getPosX(n)
    posY = nodes.getPosY(n)
    posZ = nodes.getPosZ(n)
    if (((posX>=L[0]-wallthickness) or (posX <= wallthickness)) or ((posY>=L[1]-wallthickness) or (posY <= wallthickness))):
        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 5   # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, -10, 0)
parameters.setCFL(0.0001)
parameters.setDampingCoef(0.5)
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


# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    vtpfname='dambreak2D'+str(step)+'.vtp'
    if step % 10000 == 0:
        vtp.Write(vtpfname)
print("Total Simulation Duration (s)")
#// ///////////////////////////////////////////////////







