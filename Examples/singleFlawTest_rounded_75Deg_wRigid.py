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
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTIWriter/build/')
#sys.path.insert(0,'/f/sjr/VTIWriter/build/')
#import VTIWriter
import math
## GEOMETRY ##

# Set the resolution and the pack
cellsize = 0.00025
psep = cellsize/2.0
# Background Grid
L = [0.3,0.2,0.0]
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
rock.setCriticalStrain(0.1)
rock.setFrictionAngle(55.81)
rock.setDilatancyAngle(13.95)
rock.setIsDamage(True)
rock.setDamageModel(1)
rock.setCrackSpeed(586.52)
rock.setDamageM(9)
rock.setDamageK(5.39e36)#(7.0e35)
rock.setModel(4)

platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [rock, platen]
# create packers and pack
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft = [0.05,0.01,0.0]
solidUpperRight = [0.05+0.0762,0.01+0.1524,0.0]
solid = pyck.Cuboid(1,solidLowerLeft,solidUpperRight)
pack.AddShape(solid)
theta = 75.0
flaw_a= [0.05+(0.5*0.0762 - 0.5*0.0127*math.cos(theta*3.14159/180)),0.01+0.5*0.1524 - 0.5*0.0127*math.sin(theta*3.14159/180),0.0]
flaw_b = [flaw_a[0]+0.0127*math.cos(theta*2.14159/180),flaw_a[1]+0.0127*math.sin(theta*3.14159/180),0.0]
flaw_c= [flaw_b[0],flaw_b[1]+0.00127,0.0]
flaw_d = [flaw_a[0],flaw_a[1]+0.00127,0.0]
flaw = pyck.ConvexHull2D(0,[flaw_a,flaw_b,flaw_c,flaw_d])
roundR = 0.00127/2
x0 = 0.05 + 0.5*0.0762
y0 = 0.01 + 0.5*0.1524 + cellsize
c1x = flaw_a[0]
c2x = flaw_b[0]
c1y = flaw_a[1]+0.5*0.00127
c2y = flaw_b[1]+0.5*0.00127
round1 = pyck.Cylinder(0,[c1x,c1y,2*cellsize],roundR,[0,0,0.1*cellsize])
round2 = pyck.Cylinder(0,[c2x,c2y,2*cellsize],roundR,[0,0,0.1*cellsize])

upperPlatenLowerLeft = [solidLowerLeft[0]-2*cellsize, solidUpperRight[1], solidLowerLeft[2]]
upperPlatenUpperRight = [solidUpperRight[0]+2*cellsize, upperPlatenLowerLeft[1] +2*cellsize, upperPlatenLowerLeft[2]]
upperPlaten = pyck.Cuboid(2,upperPlatenLowerLeft,upperPlatenUpperRight)
lowerPlatenLowerLeft = [solidLowerLeft[0]-2*cellsize, solidLowerLeft[1]-2*cellsize, solidLowerLeft[2]]
lowerPlatenUpperRight = [solidUpperRight[0]+2*cellsize, solidLowerLeft[1],solidLowerLeft[2]]
lowerPlaten = pyck.Cuboid(3,lowerPlatenLowerLeft,lowerPlatenUpperRight)
pack.AddShape(flaw)  # hole
pack.AddShape(round1)# hole
pack.AddShape(round2)# hole
pack.AddShape(upperPlaten)
pack.AddShape(lowerPlaten)
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
model.SetIntField(matID,2,1)
model.SetIntField(matID,3,1)
# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,1) 
model.SetIntField(objID,2,2)
model.SetIntField(objID,3,3)  
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
	idOb = matpoints.getObjID(p)
	if(idOb == 2):
               matpoints.setVel(p,0,-0.20,0)
	if(idOb == 3):
	       matpoints.setVel(p,0,0.0,0)

# //// NODAL BOUNDARY CONDITIONS /////
#for  n in range(nodes.numNodes):
#   posY = nodes.getPosY(n)
#   if (posY < 0.25):
#       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 3.0e-3  # //(time in seconds)
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

    vtpfname='paper_output/singleFlaw_rounded_75Deg_wRigid'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 5000 == 0):
        vtp.Write(vtpfname)
        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







