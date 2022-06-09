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
import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.004
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [0.20,0.20,0.20]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
R = 0.05
cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
solid = pyck.Sphere(1,[0.1,0.1,0.1],R)
topPress = pyck.Cylinder(2,[],R/4,[0,0,R/4])
bottomPress = pyck.Cylinder(3,[],R/4,[0,0,-R/4])
### PYCK - Add all shapes to pack
pack.AddShape(solid)
pack.AddShape(topPress)
pack.AddShape(bottomPress)
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
        matpoints.setVel(p,0.0,-0.05,0)
    if(obID == 3):
        matpoints.setVel(p,0.0,0.05,0)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 10e-3
is3D = True
parameters = graphyt.Parameters(tmax, is3D)
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
print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    vtpfname='output/brazillian3D'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 10000 == 0):
        vtp.Write(vtpfname)
print("Total Simulation Duration (s)")
