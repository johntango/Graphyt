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
import math
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.020
psep = cellsize/2.0
# Background Grid
L = [1.0,1.0,4*cellsize]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(1e6)
water.setShearMod(0.0)
#water.setViscosity(1000.0)
water.setDensity(10)
water.setModel(2)
rock = graphyt.MaterialModels()
rock.setBulkMod(1e6)
rock.setShearMod(1e3)
rock.setDensity(500)
rock.setModel(1)
materials = [rock,water]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
box = pyck.Cuboid(1,[2.5*cellsize,2.5*cellsize,2*cellsize],[L[0]-2.5*cellsize,L[1]-2.5*cellsize,2*cellsize])
box_gap = pyck.Cuboid(0,[4*cellsize,4*cellsize,2*cellsize],[L[0]-4*cellsize,L[1]-2.5*cellsize,2*cellsize])
dam = pyck.Cuboid(2,[5.5*cellsize,5.5*cellsize,2*cellsize],[0.5*L[0],0.5*L[1]-4*cellsize,2*cellsize])
# Add shapes to pack
pack.AddShape(box)
pack.AddShape(box_gap)
pack.AddShape(dam)
pack.Process()
model = pyck.Model(pack)
# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,0) 
model.SetIntField(objID,2,1) 
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
model.SetIntField(matID,2,1) 
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
for n in range(nodes.numNodes):
    posX = nodes.getPosX(n)
    posY = nodes.getPosY(n)
    if (((posX>=L[0]-wallthickness) or (posX <=wallthickness)) or ((posY>=L[1]-wallthickness) or (posY <=wallthickness)) ):
         #  // Want a format like this: if(isInside(nodes,n,walls))
        nodes.setFixedBC(n)
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0,-10,0)
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

# ## Velocity of Pelton Wheel
# wheelY=0.15
# wheelX= 0.105
# omega = 20.0 # rad/s

# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.
    # enforce rotation
    matpoints.setRotatingBodyVel(0,0,0,0,0)
    vtpfname='output/fsiDamBreak2D'+str(step)+'.vtp'
    if step % 5000 == 0:
        vtp.Write(vtpfname)


print("Total Simulation Duration (s)")
#// ///////////////////////////////////////////////////







