
# Import the necessary modules
import sys
# sys.path.insert(0,'/f/sjr/Graphyt/build/')
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
sys.path.insert(0,'/Users/samraymond/Desktop/Software_Development/Graphyt/build/')
# sys.path.insert(0,'C:/Users/sjr/Desktop/PhD/Graphyt/build-win/Release/')
import graphyt
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
# sys.path.insert(0,'C:/Users/sjr/Desktop/PhD\pyck/build-win/swig/Release')
sys.path.insert(0,'/Users/samraymond/Desktop/Software_Development/pyck/build/swig')
import pyck
sys.path.insert(0,'/Users/samraymond/Desktop/Software_Development/VTKWriter/build/')
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
# sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/VTKWriter/build')
import VTKWriter
import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.001
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [0.2,0.2,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
R = 0.027
solid = pyck.Sphere(1,[40e-3+R/2.0, 40e-3+R/2.0,0.0],R)
topPlaten_lowerLeft = [0.04+R/2.0-3*cellsize,0.04+R/2.0+R,0.0]
topPlaten_upperRight = [0.04+R/2.0+3*cellsize,0.04+R/2.0+R+0.0025,0.0] 
topPlaten = pyck.Cuboid(2,topPlaten_lowerLeft,topPlaten_upperRight)
bottomPlaten_lowerLeft = [0.04+R/2.0-3*cellsize,0.04+R/2.0-R-0.0025,0.0]
bottomPlaten_upperRight = [0.04+R/2.0+3*cellsize,0.04+R/2.0-R,0.0] 
bottomPlaten = pyck.Cuboid(3,bottomPlaten_lowerLeft,bottomPlaten_upperRight)
### PYCK - Add all shapes to pack
pack.AddShape(solid)
pack.AddShape(topPlaten)
pack.AddShape(bottomPlaten)
pack.Process()
### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos,numParticles)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)
print("+++++ GRID INFO+++++ No. Nodes:", nodes.numNodes)

################################################################
#--------------------------------------------------------------#
#                        Materials                             #
#--------------------------------------------------------------#
################################################################

### PYCK - Materials
rock = graphyt.MaterialModels()
rock.setBulkMod(3.31e9)
rock.setShearMod(2.48e9)
rock.setDensity(1540.0)
rock.setCohesion(5.20e6)
rock.setCriticalStrain(0.1)
rock.setFrictionAngle(55.81)
rock.setDilatancyAngle(55.81)#(13.95)
rock.setIsDamage(True)
rock.setDamageModel(1)
rock.setCrackSpeed(586.52)
rock.setDamageM(9)
rock.setDamageK(7.0e35)#(5.39e36)
rock.setModel(4)
platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [rock, platen]
### PYCK - Create Model for Simulation
model = pyck.Model(pack)
matID = model.CreateIntField("material",1)
objID = model.CreateIntField("objectID",1)
# Assign material IDs to each object
model.SetIntField(matID,1,0)
model.SetIntField(matID,2,1)
model.SetIntField(matID,3,1)
matIDfield = model.GetIntField(matID)
matpoints.addPyckIntField(matIDfield)
# Assign object IDs to each object
model.SetIntField(objID,1,1)
model.SetIntField(objID,2,2)
model.SetIntField(objID,3,3)
objIDField = model.GetIntField(objID)
matpoints.addPyckIntField(objIDField)
################################################################
#--------------------------------------------------------------#
#                     Boundary Conditions                      #
#--------------------------------------------------------------#
################################################################

### Material Point BCs
for p in range(numParticles):
    obID = matpoints.getObjID(p)
    if(obID == 2):
        matpoints.setVel(p,0.0,-0.2,0)
    if(obID == 3):
        matpoints.setVel(p,0.0,0.2,0)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 2e-3
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
parameters.setFriction(0.0)
parameters.setCFL(0.1)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")
################################################################
#--------------------------------------------------------------#
#                   Simulation Solver                          #
#--------------------------------------------------------------#
################################################################

## Create the solver object
solver = graphyt.Solver(nodes, matpoints, parameters, materials)
## Define the output arrays
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
vtp = VTKWriter.VTPWriter(numParticles)
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
print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
# nsteps = 10000
# solver.iterate_CPU_steps(nsteps, parameters)
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    vtpfname='output/brazillian2D'+str(step)+'.vtp'
    vtifname='output/brazillian2Dvti'+str(step)+'.vti'
    if (step % 1000 == 0):
       vtp.Write(vtpfname)
       vti.Write(vtifname)

print("Total Simulation Duration (s)")
