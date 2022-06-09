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
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.025
psep = cellsize/2.0
# Background Grid
L = [50*cellsize, 50*cellsize,0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(1e9)
water.setDensity(1000)
water.setViscosity(10000)
water.setModel(2)
rock = graphyt.MaterialModels()
rock.setBulkMod(10e9)
rock.setShearMod(2e9)
rock.setDensity(5000)
rock.setModel(1)
materials = [water,rock]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0,0,0])
pack = pyck.StructuredPack(cubic)
# create shapes
solid1LowerLeft = [2*cellsize,4*cellsize,0.0]
solid1UpperRight = [L[0]-2*cellsize,L[1]/4,0.0]
solid1 = pyck.Cuboid(1,solid1LowerLeft,solid1UpperRight)
pack.AddShape(solid1)
liquidLowerLeft = [2*cellsize,L[1]/4,0.0]
liquidUpperRight = [L[0]-2*cellsize,3*L[1]/4,0.0]
liquid = pyck.Cuboid(2,liquidLowerLeft,liquidUpperRight)
pack.AddShape(liquid)
solid2LowerLeft = [2*cellsize,3*L[1]/4,0.0]
solid2UpperRight = [L[0]-2*cellsize,L[1]-4*cellsize,0.0]
solid2 = pyck.Cuboid(3,solid2LowerLeft,solid2UpperRight)
pack.AddShape(solid2)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,1) 
model.SetIntField(matID,2,0) 
model.SetIntField(matID,3,1) 
# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
matpoints.addPyckIntField(matIDfield)

print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////


# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posY = nodes.getPosY(n)
    if ((posY>L[1]-7.5*cellsize ) or (posY < 8*cellsize )):
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

# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/pouiselleCOUPLEDVTP'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 50000 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







