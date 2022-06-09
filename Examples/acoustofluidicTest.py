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
cellsize = 0.001
psep = cellsize/2.0
# Background Grid
L = [50*cellsize, 40*cellsize, 0.0]
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
solid = graphyt.MaterialModels()
solid.setBulkMod(6.48e9)
solid.setShearMod(1.24e9)
solid.setDensity(1050)
solid.setModel(1)
materials = [water]#,solid]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0,0,0])
pack = pyck.StructuredPack(cubic)
# # create shapes
# R = 2*cellsize
# bead = pyck.Cylinder(1,[L[0]/3,L[1]/2.5,0.0],R,[0,0,0.1*cellsize])

liquidLowerLeft = [2*cellsize,L[1]/4,-1]
liquidUpperRight = [L[0]-2*cellsize,3*L[1]/4,1]
liquid = pyck.Cuboid(2,liquidLowerLeft,liquidUpperRight)
pack.AddShape(liquid)
# pack.AddShape(bead)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,1) 
model.SetIntField(matID,2,0) 
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,1) # Set the liquid to mat 1
model.SetIntField(objID,2,2) 
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
for p in range(matpoints.numParticles):
    posY = matpoints.getPosY(p)
    if (posY>L[1]-12*cellsize):
        matpoints.setBC(p,graphyt.BCTypes.mp_vysin,[0.001,3*166376.5/5.5])

# //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
for  n in range(nodes.numNodes):
    posY = nodes.getPosY(n)
    posX = nodes.getPosX(n)
    # if (posY>L[1]-13*cellsize):
    #     nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
        # if((posX>1*L[0]/4) and (posX < 3*L[0]/4)):
        #nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-1])
    if (posY < 13*cellsize ):
        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0e-4  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, 0, 0)
parameters.setCFL(0.001)
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

    vtpfname='output/acoustofluidics_f3'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 1000 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







