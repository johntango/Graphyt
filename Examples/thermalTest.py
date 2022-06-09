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
import math

# Set the resolution and the pack
cellsize = 0.01
psep = cellsize/2.0
# Background Grid
L = [1.0,1.0,4*cellsize]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK
# create materials
copper = graphyt.MaterialModels()
copper.setBulkMod(1e7)
copper.setShearMod(1000)
copper.setDensity(1)
copper.setThermalConductivity(391)
copper.setSpecificHeatCapacity(385.2)
copper.setModel(1)
materials = [copper]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
cube = pyck.Cuboid(1,[0.2,0.2,2*cellsize],[0.8,0.8,2*cellsize])
# Add shapes to pack
pack.AddShape(cube)
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
# // Now that nodes and material points have been discretised, boundary conditions are applied
# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
for p in range(matpoints.numParticles):
       # Set Temp for bottom layer
    if(matpoints.getPosY(p)<0.3):
       matpoints.setTemp(p,50)
    if(matpoints.getPosY(p)>0.72):
       matpoints.setTemp(p,100)
# //// NODAL BOUNDARY CONDITIONS /////
#     //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posX = nodes.getPosX(n)
    posY = nodes.getPosY(n)
    posZ = nodes.getPosZ(n)
    if (((posX> (L[0]-wallthickness)) or (posX < wallthickness )) or ((posY> (L[1]-wallthickness)) or (posY < wallthickness)) ):
        #  // Want a format like this: if(isInside(nodes,n,walls))
        nodes.setFixedBC(n)
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
# /////////// Parameters Set /////////////// #
tmax = 100.0  # //(time in seconds)
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
matpntTemp = matpoints.getTemperatureArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Temperature",matpntTemp,1, 2)

#//     // Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters) 
    #//     // After nsteps, output the data using serialize for nodes and material points.
    vtpfname='output/thermalTest'+str(step)+'.vtp'
    if step % 2500 == 0:
        vtp.Write(vtpfname)
    for p in range(matpoints.numParticles):
       # Set Temp for bottom layer
        if(matpoints.getPosY(p)<0.3):
            matpoints.setTemp(p,50)
        if(matpoints.getPosY(p)>0.72):
            matpoints.setTemp(p,100)
print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////
