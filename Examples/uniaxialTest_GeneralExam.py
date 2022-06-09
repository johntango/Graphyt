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
cellsize = 0.005
psep = cellsize/2.0
# Background Grid
L = [0.2,0.45,0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials

rock = graphyt.MaterialModels()
rock.setBulkMod(28.79e6)
rock.setShearMod(17.27e6)
rock.setDensity(3786.8)
rock.setModel(1)
platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [rock]
# create packers and pack
cubic = pyck.CubicPacker(L,psep,[0,0,0])
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft = [0.11,0.225,0.0]
solidUpperRight = [0.11+4*cellsize,0.225+16*cellsize,0.0]
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
# for p in range(matpoints.numParticles):  
# 	y = matpoints.getPosY(p)
#         if(y<solidLowerLeft[1]+cellsize):
#             matpoints.setStress(p,0,-10,0,0,0,0)
#         if(y<solidUpperRight[1]-cellsize):
#             matpoints.setStress(p,0,-10,0,0,0,0)


# for p in range(matpoints.numParticles):  
# 	idOb = matpoints.getObjID(p)
# 	y = matpoints.getPosY(p)
# 	if(idOb==3):
#             matpoints.setVel(p,0.0,-0.001,0.0)
# 	if(idOb==2):
#             matpoints.setVel(p,0.0,0.001,0.0)

# //// NODAL BOUNDARY CONDITIONS /////
#for  n in range(nodes.numNodes):
#   posY = nodes.getPosY(n)
#   if (posY < 0.25):
#       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
#   if (posY > 0.37):
#       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-0.015])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0, 0, 0)
parameters.setCFL(0.1)
parameters.setDampingCoef(0.002)
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
matpntStrain = matpoints.getStrainArr()
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
vtp.AddArray("Strain",matpntStrain,6, 2)
vtp.AddArray("Normals",matpntNormal,3,2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
vtp.AddArray("Smax",matpntSMax,1,2)
# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    # for p in range(matpoints.numParticles):  
	#     y = matpoints.getPosY(p)
    #         if(y<solidLowerLeft[1]+2*cellsize):
    #             matpoints.setVel(p,0,1/(step+1),0)
    #         elif(y>solidUpperRight[1]-2*cellsize):
    #             matpoints.setVel(p,0,-1/(step+1),0)
    speed = 0.01# - step*0.000001
    if step > 1000:
        speed = 0
    
    for  n in range(nodes.numNodes):
        posY = nodes.getPosY(n)
        if (posY < solidLowerLeft[1]+3*cellsize):
            nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,speed])
        elif (posY > solidUpperRight[1]-2*cellsize):
            nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-speed])
    vtpfname='output/ge2Results'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if step % 2 == 0:
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







