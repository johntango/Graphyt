
# Import the necessary modules
#import sys
#sys.path.insert(0,'/f/sjr/Graphyt/build/')
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/Graphyt/build/')
#sys.path.insert(0,'C:/Users/sjr/Desktop/PhD/Graphyt/build-win/Release/')
import graphyt
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/pyck/build/swig/')
#sys.path.insert(0,'C:/Users/sjr/Desktop/PhD\pyck/build-win/swig/Release')
#sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
#sys.path.insert(0,'/f/sjr/VTPWriter/build/')
#sys.path.insert(0,'C:/Users/sjr/Desktop/PhD/VTPWriter/build-win/Release')
#sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/PhD/VTPWriter/build/')
import VTKWriter
import math
import numpy as np
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
import genFracture
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
L = [0.10,0.10,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
botBox = pyck.Cuboid(4,[7*cellsize,5*cellsize,0],[L[0]-7*cellsize,L[1]/2 - 2*cellsize,0.0])
topBox = pyck.Cuboid(3,[7*cellsize,L[1]/2+2*cellsize,0],[L[0]-7*cellsize,L[1]-5*cellsize,0.0])
solid_bot = pyck.Cuboid(2,[11*cellsize,9*cellsize,0],[L[0]-11*cellsize,L[1]/2,0.0])
solid_top = pyck.Cuboid(1,[11*cellsize,L[1]/2,0],[L[0]-11*cellsize,L[1]-9*cellsize,0.0])

### PYCK - Add all shapes to pack
pack.AddShape(topBox)
pack.AddShape(botBox)
pack.AddShape(solid_top)
pack.AddShape(solid_bot)

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
rock.setBulkMod(2*3.31e7)
rock.setShearMod(2*2.48e7)
rock.setDensity(1540.0)
rock.setCohesion(5.20e6)
rock.setCriticalStrain(0.1)
rock.setFrictionAngle(55.81)
rock.setDilatancyAngle(55.81)#(13.95)
rock.setIsDamage(True)
rock.setDamageModel(-1)
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
model.SetIntField(matID,2,0)
model.SetIntField(matID,3,1)
model.SetIntField(matID,4,1)

matIDfield = model.GetIntField(matID)
matpoints.addPyckIntField(matIDfield)
# Assign object IDs to each object
model.SetIntField(objID,1,1)
model.SetIntField(objID,2,2)
model.SetIntField(objID,3,3)
model.SetIntField(objID,4,4)

objIDField = model.GetIntField(objID)
matpoints.addPyckIntField(objIDField)
## Refine interface using a sine wave boundary and change the objectID based on the location
# wavelength = L[0] - 21*cellsize
# for p in range(matpoints.numParticles):
#     ypos = matpoints.getPosY(p)
#     xpos = matpoints.getPosX(p) + 8*cellsize
#     yint = 0.5*cellsize*math.sin(4.*2.*3.14159*xpos/wavelength) + L[1]/2. 
#     if((matpoints.getObjID(p) is not 3)  and (matpoints.getObjID(p) is not 4) ):
#         if(ypos > yint):
#             matpoints.setObjectID(p,1)
#         else:
#             matpoints.setObjectID(p,2)
## Redefine the interface using the fracture appeture from Morris
zeta = 0.35
nx = 1  + int(round((L[0] - 14*cellsize)/cellsize) )
appeture = genFracture.drazerKoplik(nx=nx,ny=1,zeta=zeta,seed=2)
target = find_nearest(appeture,0.05)
err = abs(appeture[0] - target) 
while (err> cellsize/2 ):
    appeture = np.roll(appeture,1)
    err = abs(appeture[0] - target) 

# print(appeture)
for p in range(matpoints.numParticles):
    ypos = matpoints.getPosY(p)
    xpos = matpoints.getPosX(p) - 7*cellsize
    cellID_x = round(xpos/cellsize)
    yint =  0.004*appeture[cellID_x] + 0.49*L[1]
    if((matpoints.getObjID(p) is not 3)  and (matpoints.getObjID(p) is not 4) ):
        if(ypos > yint):
            matpoints.setObjectID(p,1)
        else:
            matpoints.setObjectID(p,2)



################################################################
#--------------------------------------------------------------#
#                     Boundary Conditions                      #
#--------------------------------------------------------------#
################################################################

### Material Point BCs
for p in range(numParticles):
    obID = matpoints.getObjID(p)
    if(obID == 4):
        # matpoints.setVel(p,0.05,0.0,0)
        matpoints.setBC(p,graphyt.BCTypes.mp_fx,[0.10])
    if(obID == 3):
        # matpoints.setVel(p,-0.05,0.0,0)
        matpoints.setBC(p,graphyt.BCTypes.mp_fx,[-0.10])
## Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 80e-3
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
parameters.setFriction(0.0)
parameters.setDampingCoef(100.750)
parameters.setCFL(0.05)
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
matpntObjID = matpoints.getObjectIDArr()
Le = [1+int(L[0]/cellsize), 1+int(L[1]/cellsize), 1+int(L[2]/cellsize)]
## Create a VTP writer and add the arrays
vtp = VTKWriter.VTPWriter(numParticles)
vtp.AddPositions("Position", matpntPos,3, 2)
vtp.AddArray("Displacement",matpntDisp,3,2)
vtp.AddArray("Material", matpntMat,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("Normals",matpntNormal,3,2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
vtp.AddArray("Smax",matpntSMax,1,2)
vtp.AddArray("ObjectID",matpntObjID,1,2)

print("+++++ Solver  INFORMATION+++++ Solver Ready")

################################################################
#--------------------------------------------------------------#
#                   Time Integration                           #
#--------------------------------------------------------------#
################################################################
for step in range(int(solver.tmax/solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    vtpfname='output/shearBoxLF'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 1000 == 0):
        vtp.Write(vtpfname)
print("Total Simulation Duration (s)")




