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
cellsize = 0.01
psep = cellsize/2.0
# Background Grid
L = [0.2,0.45,4*cellsize]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials

rock = graphyt.MaterialModels()
rock.setBulkMod(28.79e9)
rock.setShearMod(17.27e9)
rock.setDensity(3786.8)
rock.setYieldStress(2e7)
rock.setModel(3)
materials = [rock]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft = [0.11,0.225,2*cellsize]
solidUpperRight = [0.158,0.3845,2*cellsize]
solid = pyck.Cuboid(1,solidLowerLeft,solidUpperRight)
pack.AddShape(solid)

pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 

# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
# objIDfield = model.GetIntField(objID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
matpoints.addPyckIntField(matIDfield)
# matpoints.addPyckIntField(objIDfield)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)

# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////


# //// NODAL BOUNDARY CONDITIONS /////
for  n in range(nodes.numNodes):
   posY = nodes.getPosY(n)
   if (posY < 0.25):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
   if (posY > 0.37):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-1.015])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 10.0e-3  # //(time in seconds)
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
matpntPlasStr = matpoints.getPlasticStrainArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("PlasticStrain",matpntPlasStr,6,2)
# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/perfectplasticCompression'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 1000 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







