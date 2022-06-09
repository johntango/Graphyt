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
cellsize = 0.005
psep = cellsize/2.0
# Background Grid
L = [0.3,0.3,4*cellsize]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials

rock = graphyt.MaterialModels()
rock.setBulkMod(3.31e9)
rock.setShearMod(2.48e9)
rock.setDensity(1540.0)
rock.setCohesion(5.20e6)
rock.setCriticalStrain(0.01)
rock.setFrictionAngle(55.81)
rock.setDilatancyAngle(13.95)
rock.setIsDamage(True)
rock.setDamageModel(1)
rock.setCrackSpeed(586.52)
rock.setDamageM(9)
rock.setDamageK(5.39e36)
rock.setModel(4)
materials = [rock]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft = [0.05,0.01,2*cellsize]
solidUpperRight = [0.05+0.0762,0.01+0.1524,2*cellsize]
solid = pyck.Cuboid(1,solidLowerLeft,solidUpperRight)
pack.AddShape(solid)
flawLowerLeft = [0.05+0.03175,0.01+0.0699,2*cellsize]
flawUpperRight = [0.05+0.04445,0.01+0.0712,2*cellsize]
flaw = pyck.Cuboid(0,flawLowerLeft,flawUpperRight)
pack.AddShape(flaw)

pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
#model.SetIntField(matID,2,0)
#model.SetIntField(matID,3,0)
# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,0) 
#model.SetIntField(objID,2,1)
#model.SetIntField(objID,3,2)  
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
#for p in range(matpoints.numParticles):  
#	idOb = matpoints.getObjID(p)
#	y = matpoints.getPosY(p)
#	if(y>0.335):
#               matpoints.setBC(p,graphyt.BCTypes.mp_vy,[-0.0200])
#	if(y<0.275):
#	       matpoints.setBC(p,graphyt.BCTypes.mp_vy,[0.0])

# //// NODAL BOUNDARY CONDITIONS /////
for  n in range(nodes.numNodes):
   posY = nodes.getPosY(n)
   if (posY < 0.011):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
   if (posY > 0.15):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0.1])
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
matpntAccPlasStr = matpoints.getAccPlasticStrainArr()
matpntDamage = matpoints.getDamageArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
vtp = VTPWriter.VTPWriter(Le,cellsize,numParticles)
vtp.AddArray("Position", matpntPos,3, 2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
# Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #// After nsteps, output the data using serialize for nodes and material points.

    vtpfname='output/tensile_singleFlaw'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 50 == 0):
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







