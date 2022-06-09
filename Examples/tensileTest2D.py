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
# sys.path.insert(0,'/f/sjr/VTPWriter/build/')
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/VTKWriter/build')
import VTKWriter
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.05
psep = cellsize/2.0
# Background Grid
L = [2,2,0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials

# rock = graphyt.MaterialModels()
# rock.setBulkMod(28.79e6)
# rock.setShearMod(17.27e6)
# rock.setDensity(3786.8)
# rock.setModel(1)
steel = graphyt.MaterialModels()
steel.setBulkMod(215e9)
steel.setShearMod(77e9)
steel.setDensity(7500)
steel.setYieldStress(250e6)
steel.setModel(1)
platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [steel,platen]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft = [0.5,0.5,0.0]
solidUpperRight = [1.5-0.25*cellsize,1.5-0.25*cellsize,0.0]
solid = pyck.Cuboid(1,solidLowerLeft,solidUpperRight)
pack.AddShape(solid)

pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 

# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,1) 

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
    y = matpoints.getPosY(p)
    if(y<0.6):
            matpoints.setBC(p,graphyt.BCTypes.mp_fy,[-1000])
    if(y>1.40):
            matpoints.setBC(p,graphyt.BCTypes.mp_fy,[1000])

# //// NODAL BOUNDARY CONDITIONS /////
#for  n in range(nodes.numNodes):
#   posY = nodes.getPosY(n)
#   if (posY < 0.25):
#       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
#   if (posY > 0.37):
#       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-0.015])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 100.0e-3  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, 0, 0)
#parameters.setDampingCoef(20)
parameters.setCFL(0.95)
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
vtp = VTKWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddPositions("Position", matpntPos,3, VTKWriter.Order.ijk)
vtp.AddArray("Displacement",matpntDisp,3,VTKWriter.Order.ijk)
vtp.AddArray("Material", matpntMat,1, VTKWriter.Order.ijk)
vtp.AddArray("Mass", matpntMass,1, VTKWriter.Order.ijk)
vtp.AddArray("Velocity",matpntVel,3, VTKWriter.Order.ijk)
vtp.AddArray("Stress",matpntStress,6, VTKWriter.Order.ijk)
vtp.AddArray("Normals",matpntNormal,3,VTKWriter.Order.ijk)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,VTKWriter.Order.ijk)
vtp.AddArray("Damage",matpntDamage,1,VTKWriter.Order.ijk)
vtp.AddArray("Smax",matpntSMax,1,VTKWriter.Order.ijk)

gridMass = nodes.getMassArr()
gridForce = nodes.getForceArr()
gridVel = nodes.getVelArr()
vti = VTKWriter.VTIWriter(Le,[-psep,-psep,-psep],cellsize)
vti.AddArray("Mass",gridMass,1, VTKWriter.Order.ijk)
vti.AddArray("Force",gridForce,3, VTKWriter.Order.ijk)
vti.AddArray("Velocity",gridVel,3,VTKWriter.Order.ijk)
# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/tensileTest2D'+str(step)+'.vtp'
    vtifname='output/tensileTest2Dvti'+str(step)+'.vti'
    if step % 1000 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







