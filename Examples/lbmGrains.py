###########################################################
 
# USER INPUTS FOR GRAIN MODEL
 
# Loading
v = 1.0
maxDisplacement = 0.01
tension = False
 
# Material Properties
# material 1
K1 = 3e7
G1 = 1e7
rho1 = 2000
# material 2
K2 = 0.2 * K1
G2 = 0.2 * G1
rho2 = rho1
###########################################################
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
import math
import numpy
import math
import csv
from itertools import product
 
## GEOMETRY ##
# Set the resolution and the pack
cellsize = 0.01
psep = cellsize / 2.0
# Background Grid
L = [1.2, 1.2, 0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK
# create materials
mat1 = graphyt.MaterialModels()
mat1.setBulkMod(K1)
mat1.setShearMod(G1)
mat1.setDensity(rho1)
mat1.setModel(1)
mat2 = graphyt.MaterialModels()
mat2.setBulkMod(K2)
mat2.setShearMod(G2)
mat2.setDensity(rho2)
mat2.setCohesion(5.20e6)
mat2.setCriticalStrain(0.1)
mat2.setFrictionAngle(55.81)
mat2.setDilatancyAngle(55.81)#(13.95)
mat2.setIsDamage(True)
mat2.setDamageModel(1)
mat2.setCrackSpeed(586.52)
mat2.setDamageM(9)
mat2.setDamageK(7.0e35)#(5.39e36)
mat2.setModel(4)
 
platen = graphyt.MaterialModels()
platen.setIsRigid(True)
materials = [mat1, mat2, platen]
# create packers and pack
cubic = pyck.CubicPacker(L, psep, [0.25 * cellsize, 0.25 * cellsize, 0])
pack = pyck.StructuredPack(cubic)
 
# create shapes
sampleLeft = [0.1, 0.1, 0.0]
sampleRight = [1.1, 1.1, 0.0]
sample = pyck.Cuboid(1, sampleLeft, sampleRight)
pack.AddShape(sample)
 
bottomPlatenLeft = [0.1 - 1.5 * cellsize, 0.1 - 1.5 * cellsize, 0.0]
bottomPlatenRight = [1.1 + 1.5 * cellsize, 0.1, 0.0]
bottomPlaten = pyck.Cuboid(2, bottomPlatenLeft, bottomPlatenRight)
pack.AddShape(bottomPlaten)
 
topPlatenLeft = [0.1 - 1.5 * cellsize, 1.1, 0.0]
topPlatenRight = [1.1 + 1.5 * cellsize, 1.1 + 1.5 * cellsize, 0.0]
topPlaten = pyck.Cuboid(3, topPlatenLeft, topPlatenRight)
pack.AddShape(topPlaten)
 
pack.Process()
model = pyck.Model(pack)
# Set the material ID for all objects
matID = model.CreateIntField("material", 1)
model.SetIntField(matID, 1, 0)
model.SetIntField(matID, 2, 2)
model.SetIntField(matID, 3, 2)
 
# Set the objectID's for all objects
objID = model.CreateIntField("objectID", 1)
model.SetIntField(objID, 1, 0)
model.SetIntField(objID, 2, 1)
model.SetIntField(objID, 3, 2)
 
# Add the geometry from pyck to the matpoints
matIDfield = model.GetIntField(matID)
objIDfield = model.GetIntField(objID)
numParticles = pack.GetNumParticles()
pyckPos = pack.GetPositions()
matpoints.addPyckGeom(pyckPos, numParticles)
# Testing material by cell ID
matpoints.addPyckIntField(matIDfield)
 
# material_array = numpy.arange(0,numParticles,dtype=int)
# material_array = [0] * numParticles
# matpoints.setMaterialByCell(material_array)
N = [int(math.ceil(L[0]/cellsize)),int(math.ceil(L[1]/cellsize)),0]
material_grid = [0] * (N[0] * N[1])
 
grain_map = [0]*100*100
with open('6_grain_map.csv', 'r') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\n')
    idx = 0
    for row in spamreader:
        grain_map[idx] = int(row[0])
        idx = idx+1
 
def expandLabel(grainMap, label, L):
    if(L[2] ==0): L[2] = 1
 
    tmpLabels = [0] * (L[0] * L[1] * L[2])
    for i in range(L[0] * L[1] * L[2]):
        tmpLabels[i] = grainMap[i]
 
    kMin = 0
    kMax = L[2]
    if L[2] < 2:
        kMin = 0
        kMax = 1
 
    for k, j, i in product(range(kMin, kMax), range(0, L[1]), range(0, L[0])):
        idx = i + (j + k * L[1]) * L[0]
        thisLabel = grainMap[idx]
        if(thisLabel == label):
            nkMin = k - 1
            nkMax = k + 2
            if L[2] < 2:
                nkMin = k
                nkMax = k + 1
 
            for nk2, nj2, ni2 in product(range(nkMin, nkMax), range(j - 1, j + 2), range(i - 1, i + 2)):
                ni = ni2 % L[0]
                nj = nj2 % L[1]
                nk = nk2 % L[2]
                nidx = ni + (nj + nk * L[1]) * L[0]
                tmpLabels[nidx] = label
    return tmpLabels
grain_map = expandLabel(grain_map,1,[100,100,1])
lower = int(N[1]/2)-1
upper = int(N[1]/2)
lx = int(sampleLeft[0]/cellsize)
ly = int(sampleLeft[1]/cellsize)
ux = int(sampleRight[0]/cellsize)
uy = int(sampleRight[1]/cellsize)
i2 = 0
j2 = 0
for j in range(ly,uy):
    for i in range(lx,ux):
        if(grain_map[i2+j2*100] == 1):
            matpoints.setMaterialByCell(i+j*N[0],1)
        i2 = i2+1
    i2 = 0
    j2 = j2+1
 
# matpoints.setMaterialByCell(material_grid)
 
matpoints.addPyckIntField(objIDfield)
print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)
 
# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////
if tension:
    v = -v
 

def setPlatenVel(v):
    global matpoints
    global numParticles
    for p in range(numParticles):
        obID = matpoints.getObjID(p)
        if(obID == 1):
            matpoints.setVel(p, 0.0, v, 0)
        if(obID == 2):
            matpoints.setVel(p, 0.0, -v, 0)
 
# //// NODAL BOUNDARY CONDITIONS /////
 
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
 
# /////////// Parameters Set /////////////// #
tmax = 1000000.0e-3  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax, is3D)
parameters.setGravity(0, 0, 0)
parameters.setDampingCoef(10.0)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")
 
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# Now Mat points and nodes are set and ready to go, create a solver to run
# an MPM analysis on them
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
Le = [1 + int(L[0] / cellsize), 1 + int(L[1] / cellsize),
      1 + int(L[2] / cellsize)]
# Create a VTP writer and add the arrays
vtp = VTPWriter.VTPWriter(Le, cellsize, numParticles)
vtp.AddArray("Position", matpntPos, 3, 2)
vtp.AddArray("Displacement", matpntDisp, 3, 2)
vtp.AddArray("Material", matpntMat, 1, 2)
vtp.AddArray("Mass", matpntMass, 1, 2)
vtp.AddArray("Velocity", matpntVel, 3, 2)
vtp.AddArray("Stress", matpntStress, 6, 2)
vtp.AddArray("Normals", matpntNormal, 3, 2)
vtp.AddArray("AccPlasticStrain", matpntAccPlasStr, 1, 2)
vtp.AddArray("Damage", matpntDamage, 1, 2)
vtp.AddArray("Smax", matpntSMax, 1, 2)
 

def writeFile(step):
    global vtp
    vtpfname = 'output/grainTestExpanded' + str(step) + '.vtp'
    if step % 100 == 0:
        vtp.Write(vtpfname)
        print("Output Step: "+str(step))
currV = 0.0
VStep = v/1000
 
setPlatenVel(0.0)
currDisplacement = 0.0
step = 0
while(currDisplacement < maxDisplacement):
    currV = currV+VStep
    if(currV>v):
        currV = v
    setPlatenVel(currV)
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    currDisplacement += solver.dt * v
    writeFile(step)
    step = step + 1
 
setPlatenVel(0.0)
print("Platen Stopped")
 
for step in range(step, int(solver.tmax / solver.dt)):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    writeFile(step)
 
print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////