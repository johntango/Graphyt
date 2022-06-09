# Import the necessary modules

import sys
sys.path.insert(0,'/f/sjr/Graphyt/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
import graphyt
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
sys.path.insert(0,'/f/sjr/VTIWriter/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
import VTPWriter
import math
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 5
psep = cellsize/2.0
# Background Grid
L = [1000.0,1000.0, 1000.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials
water = graphyt.MaterialModels()
water.setBulkMod(1000)
water.setShearMod(0.0)
water.setDensity(10)
water.setModel(1)
rock = graphyt.MaterialModels()
rock.setBulkMod(1000)
rock.setShearMod(1000)
rock.setDensity(10)
rock.setModel(1)
materials = [rock, water]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
leg = pyck.StlShape(1,'Football_Helmet.STL',[0.0,0.0,0.0],1.0,[0.0,0.0,0.0],0.0)
# Add shapes to pack
pack.AddShape(leg)
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
# wallthickness = 2.0*nodes.cellsize
# for  n in range(nodes.numNodes):
#     posX = nodes.getPosX(n)
#     posY = nodes.getPosY(n)
#     posZ = nodes.getPosZ(n)
#     if (((posX>=0.2-wallthickness) or (posX <= 0.01+wallthickness)) or ((posY>=0.3-wallthickness) or (posY <=0.01+wallthickness)) or ((posZ>= 0.15-wallthickness) or (posZ <= 0.05+wallthickness )) ):
#         #  // Want a format like this: if(isInside(nodes,n,walls))
#         nodes.setFixedBC(n)
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = True
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
  #  matpoints.setRotatingBodyVel(0,wheelX,wheelY,0.0,omega)
    vtpfname='output/helmet'+str(step)+'.vtp'
    if step % 250 == 0:
        vtp.Write(vtpfname)


print("Total Simulation Duration (s)")
#// ///////////////////////////////////////////////////







