
# Import the necessary modules

import graphyt
import pyck
import VTKWriter
import math
################################################################
#--------------------------------------------------------------#
#                        GEOMETRY                              #
#--------------------------------------------------------------#
################################################################

### GRAPHYT - Initialising nodes and material points
# Grid cellsize
cellsize = 0.05
# Initial particle separation
psep = cellsize/2
# Compuational Domain
L = [5.0,5.0,0.0]
# nodes and matpoint objects
nodes = graphyt.Nodes(L,cellsize)
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
### PYCK - Particle Geometry
cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0])
pack = pyck.StructuredPack(cubic)
R = 0.2
solid_lowerLeft = [10*cellsize,10*cellsize,0.0]
solid_upperRight = [L[0]-10*cellsize,L[1]-10*cellsize,0.0] 
solid = pyck.Cuboid(1,solid_lowerLeft,solid_upperRight)

solidTop_lowerLeft = [10*cellsize,L[1]-15*cellsize,0.0]
solidTop_upperRight = [L[0]-10*cellsize,L[1]-10*cellsize,0.0] 
solidTop = pyck.Cuboid(3,solidTop_lowerLeft,solidTop_upperRight)

solidBot_lowerLeft = [10*cellsize,10*cellsize,0.0]
solidBot_upperRight = [L[0]-10*cellsize,15*cellsize,0.0] 
solidBot = pyck.Cuboid(4,solidBot_lowerLeft,solidBot_upperRight)

solidRight_lowerLeft = [L[0]-15*cellsize,10*cellsize,0.0]
solidRight_upperRight = [L[0]-10*cellsize,L[1]-10*cellsize,0.0] 
solidRight = pyck.Cuboid(5,solidRight_lowerLeft,solidRight_upperRight)

solidLeft_lowerLeft = [10*cellsize,10*cellsize,0.0]
solidLeft_upperRight = [15*cellsize,L[1]-10*cellsize,0.0] 
solidLeft = pyck.Cuboid(6,solidLeft_lowerLeft,solidLeft_upperRight)

water = pyck.Sphere(2,[L[0]/2,L[1]/2,0.0],R)
### PYCK - Add all shapes to pack
pack.AddShape(solid)
pack.AddShape(solidTop)
pack.AddShape(solidBot)
pack.AddShape(solidLeft)
pack.AddShape(solidRight)

pack.AddShape(water)
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
rock.setIsDamage(False)
rock.setDamageModel(1)
rock.setCrackSpeed(586.52)
rock.setDamageM(9)
rock.setDamageK(7.0e35)#(5.39e36)
rock.setModel(1)
rock2 = graphyt.MaterialModels()
rock2.setBulkMod(3.31e9)
rock2.setShearMod(2.48e9)
rock2.setDensity(1540.0)
rock2.setCohesion(5.20e6)
rock2.setCriticalStrain(0.1)
rock2.setFrictionAngle(55.81)
rock2.setDilatancyAngle(55.81)#(13.95)
rock2.setIsDamage(False)
rock2.setDamageModel(1)
rock2.setCrackSpeed(586.52)
rock2.setDamageM(9)
rock2.setDamageK(7.0e35)#(5.39e36)
rock2.setModel(1)
rock3 = graphyt.MaterialModels()
rock3.setBulkMod(3.31e9)
rock3.setShearMod(2.48e9)
rock3.setDensity(1540.0)
rock3.setCohesion(5.20e6)
rock3.setCriticalStrain(0.1)
rock3.setFrictionAngle(55.81)
rock3.setDilatancyAngle(55.81)#(13.95)
rock3.setIsDamage(False)
rock3.setDamageModel(1)
rock3.setCrackSpeed(586.52)
rock3.setDamageM(9)
rock3.setDamageK(7.0e35)#(5.39e36)
rock3.setModel(1)
rock4 = graphyt.MaterialModels()
rock4.setBulkMod(3.31e9)
rock4.setShearMod(2.48e9)
rock4.setDensity(1540.0)
rock4.setCohesion(5.20e6)
rock4.setCriticalStrain(0.1)
rock4.setFrictionAngle(55.81)
rock4.setDilatancyAngle(55.81)#(13.95)
rock4.setIsDamage(False)
rock4.setDamageModel(1)
rock4.setCrackSpeed(586.52)
rock4.setDamageM(9)
rock4.setDamageK(7.0e35)#(5.39e36)
rock4.setModel(1)
rock5 = graphyt.MaterialModels()
rock5.setBulkMod(3.31e9)
rock5.setShearMod(2.48e9)
rock5.setDensity(1540.0)
rock5.setCohesion(5.20e6)
rock5.setCriticalStrain(0.1)
rock5.setFrictionAngle(55.81)
rock5.setDilatancyAngle(55.81)#(13.95)
rock5.setIsDamage(False)
rock5.setDamageModel(1)
rock5.setCrackSpeed(586.52)
rock5.setDamageM(9)
rock5.setDamageK(7.0e35)#(5.39e36)
rock5.setModel(1)
water = graphyt.MaterialModels()
water.setBulkMod(3.31e8)
water.setShearMod(2.48e8)
water.setDensity(1540.0)
water.setModel(1)
materials = [rock, water, rock2, rock3, rock4, rock5]
### PYCK - Create Model for Simulation
model = pyck.Model(pack)
matID = model.CreateIntField("material",1)
objID = model.CreateIntField("objectID",1)
# Assign material IDs to each object
model.SetIntField(matID,1,0)
model.SetIntField(matID,2,1)
model.SetIntField(matID,3,2)
model.SetIntField(matID,4,3)
model.SetIntField(matID,5,4)
model.SetIntField(matID,6,5)
matIDfield = model.GetIntField(matID)
matpoints.addPyckIntField(matIDfield)
# Assign object IDs to each object
model.SetIntField(objID,1,1)
model.SetIntField(objID,2,1)
model.SetIntField(objID,3,1)
model.SetIntField(objID,4,1)
model.SetIntField(objID,5,1)
model.SetIntField(objID,6,1)
objIDField = model.GetIntField(objID)
matpoints.addPyckIntField(objIDField)
################################################################
#--------------------------------------------------------------#
#                     Boundary Conditions                      #
#--------------------------------------------------------------#
################################################################

### Material Point BCs
for p in range(numParticles):
    matID = matpoints.getMatID(p)
    if(matID == 1):
        matpoints.setBC(p,graphyt.BCTypes.mp_stress_xx_yy,[-5e6,-5e6])
        matpoints.holdParticle(p)
    if(matID == 2):
        matpoints.setBC(p,graphyt.BCTypes.mp_stress_xx_yy,[0.0,-2e6])
        matpoints.holdParticle(p)
    if(matID == 3):
        matpoints.setBC(p,graphyt.BCTypes.mp_stress_xx_yy,[0.0,-2e6])
        matpoints.holdParticle(p)
    if(matID == 4):
        matpoints.setBC(p,graphyt.BCTypes.mp_stress_xx_yy,[-1e6,0.0])
        matpoints.holdParticle(p)
    if(matID == 5):
        matpoints.setBC(p,graphyt.BCTypes.mp_stress_xx_yy,[-1e6,0.0])
        matpoints.holdParticle(p)
### Nodal BCs

print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
################################################################
#--------------------------------------------------------------#
#                   Simulation Parameters                      #
#--------------------------------------------------------------#
################################################################

tmax = 4e-2
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
parameters.setFriction(0.0)
parameters.setCFL(0.1)
parameters.setDampingCoef(0.2)
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
vtp.AddArray("Volume",matpntVel,1, VTKWriter.Order.ijk)
vtp.AddArray("Density",matpntVel,1, VTKWriter.Order.ijk)
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

    vtpfname='output/pressBoreHole_steady_part'+str(step)+'.vtp'
    # vtifname='output/pressBoreHole_persistant_grid'+str(step)+'.vti'
    if (step % 500 == 0):
       vtp.Write(vtpfname)
    #    vti.Write(vtifname)

print("Total Simulation Duration (s)")
