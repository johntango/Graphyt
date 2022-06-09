# Import the necessary modules

import sys
#sys.path.insert(0,'/f/sjr/Graphyt/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
import graphyt
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
#sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
import VTPWriter
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.075
psep = cellsize/2.0
# Background Grid
L = [1, 1, 1]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
cubeLowerLeft = [cellsize/4,cellsize/4,cellsize/4]
cubeUpperRight = [1-cellsize/2,1-cellsize/2,0.5]
cube = pyck.Cuboid(1,cubeLowerLeft,cubeUpperRight)
sphere = pyck.Sphere(2,[0.5,0.5,0.8],0.2)
#veyron = pyck.StlShape(1,'Veyron.stl',[0.5,0.0,0.5],0.15,[0.0,0.0,0.0],0.0)
#boat = pyck.StlShape(1,'boat.stl',[0.0,0.0,0.0],1.0,[0.0,0.0,0.0],0.0);
# Add shapes to pack
#pack.AddShape(boat)
#pack.AddShape(veyron)
pack.AddShape(cube)
pack.AddShape(sphere)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,1) # Set the cube mat 1
model.SetIntField(matID,2,1) # Set the sphere mat 2

# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
print("adding geo")
matpoints.addPyckGeom(pyckPos,numParticles)
print("adding int fields")
matpoints.addPyckIntField(matIDfield)
matpoints.Serialise(0,"initial")
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
for p in range(matpoints.numParticles):
    if (matpoints.getPosZ(p)>0.6):
        matpoints.setVel(p,0.0,0.0,-2.0)

# //// NODAL BOUNDARY CONDITIONS /////
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
tmax = 2.0  # //(time in seconds)
is3D = True
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, -10, 0)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")

# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
#Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
print("+++++ SOLVER  INFORMATION+++++ Starting Solver")
solver = graphyt.Solver(nodes, matpoints, parameters)

# gridPos = nodes.getPosArr()
# gridMass = nodes.getMassArr()
# gridVel = nodes.getVelArr()
# Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
# vti = VTIWriter.VTIWriter(Le,cellsize,nodes.numNodes)
# vti.AddArray("Position", gridPos,3, 2)
# vti.AddArray("Mass", gridMass,1, 2)
# vti.AddArray("Velocity", gridVel,3, 2)
matpntPos = matpoints.getPosArr()
matpntMass = matpoints.getMassArr()
matpntVel = matpoints.getVelArr()
matpntStress = matpoints.getStressArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
#0//     // Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #//     // After nsteps, output the data using serialize for nodes and material points.
    filename='graphytDamBreak3DTEST'
    #vtifname='grid_graphytPyckVTI_'+str(step)+'_.vtp'
    vtpfname='output/matpnt_graphytPyckVTP_'+str(step)+'_.vtp'
    if step % 500 == 0:
        vtp.Write(vtpfname)
        #matpoints.Serialise(step,filename)
        #nodes.Serialise(step,filename)
#endTimeLoop =omp_get_wtime()#; // Get finishing time of total timeloop
#totSimTime = endTimeLoop - solver.startTimeLoop;
print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







