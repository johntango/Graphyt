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
cellsize = 0.02
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [2,5.5,2]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry

cubic = pyck.CubicPacker(L,psep)
pack = pyck.StructuredPack(cubic)
stone = pyck.StlShape(1,'kidneystone.stl',[-0.25,0.5,-1.75],1,[0.0,0.0,0.0],0.0)
topPlaten_lowerLeft = [1.25,1.50,0.2]
topPlaten_upperRight = [1.4,2.3,1.6] 
topPlaten = pyck.Cylinder(2,[1.425,2.1,0.7],0.2,[0.0,0.0,0.6])
#topPlaten = pyck.Cuboid(2,topPlaten_lowerLeft,topPlaten_upperRight)
bottomPlaten_lowerLeft = [0.1,1.5,0.2]
bottomPlaten_upperRight = [0.25,2.3,1.6] 
bottomPlaten = pyck.Cuboid(3,bottomPlaten_lowerLeft,bottomPlaten_upperRight)
### PYCK - Add all shapes to pack
pack.AddShape(stone)
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
rock.setDamageK(7.0e34)
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
v = 0.35
### Material Point BCs
for p in range(numParticles):
    obID = matpoints.getObjID(p)
    if(obID == 2):
        matpoints.setVel(p,-v,0.0,0)
    if(obID == 3):
        matpoints.setVel(p,0.00,0.0,0)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 0.20
is3D = True
parameters = graphyt.Parameters(tmax, is3D)
parameters.setDampingCoef(1.50)
parameters.setCFL(0.1)
parameters.setGravity(-98.1,0,0)
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
    solver.iterate_CPU(t, step, parameters)
    vtpfname='output/kidneystone'+str(step)+'.vtp'
    if ((step > 5000) and (step < 7000)):
        vnew = v - v*((float(step)-5000)/2000)
        for p in range(numParticles):
            obID = matpoints.getObjID(p)
            if(obID == 2):
                matpoints.setVel(p,-vnew,0.0,0)
            if(obID == 3):
                matpoints.setVel(p,vnew,0.0,0)
    if step == 7000:
        vnew = 0.0
        for p in range(numParticles):
            obID = matpoints.getObjID(p)
            if(obID == 2):
                matpoints.setVel(p,-vnew,0.0,0)
            if(obID == 3):
                matpoints.setVel(p,vnew,0.0,0)
    # if step > 5000:
    #     solver.dt+=1e-9
    #     solver.setDt(solver.dt)

    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 1000 == 0):
        vtp.Write(vtpfname)
    # if t > tmax:
    #     print(t)
    #     print(tmax)
    #     break

print("Total Simulation Duration (s)")
