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
cellsize = 0.0075
psep = cellsize/2.0
# Background Grid
L = [1.0,1.0,0.0]
nodes = graphyt.Nodes(L, cellsize)
# Material Points
matpoints = graphyt.MaterialPoints(nodes)
matpoints.psep = psep
# PYCK

# create materials

solid = graphyt.MaterialModels()
solid.setBulkMod(1e6)
solid.setShearMod(10.0)
solid.setDensity(2500.0)
solid.setModel(1)

funnel = graphyt.MaterialModels()
funnel.setIsRigid(True)

materials = [solid,solid,solid,solid,solid,solid,solid,solid,funnel]
# create packers and pack


cubic = pyck.CubicPacker(L,psep,[0,0,0])
pack = pyck.StructuredPack(cubic)
# create shapes
solidLowerLeft1 = [0.55,0.31,0.0]
solidUpperRight1 = [0.67,0.37,0.0]
solid1 = pyck.Cuboid(1,solidLowerLeft1,solidUpperRight1)
pack.AddShape(solid1)
solidLowerLeft2 = [0.55,solidUpperRight1[1],0.0]
solidUpperRight2 = [0.67,0.44,0.0]
solid2 = pyck.Cuboid(2,solidLowerLeft2,solidUpperRight2)
pack.AddShape(solid2)
solidLowerLeft3 = [0.55,solidUpperRight2[1],0.0]
solidUpperRight3 = [0.67,0.51,0.0]
solid3 = pyck.Cuboid(3,solidLowerLeft3,solidUpperRight3)
pack.AddShape(solid3)
solidLowerLeft4 = [0.55,solidUpperRight3[1],0.0]
solidUpperRight4 = [0.67,0.58,0.0]
solid4 = pyck.Cuboid(4,solidLowerLeft4,solidUpperRight4)
pack.AddShape(solid4)
solidLowerLeft5 = [0.55,solidUpperRight4[1],0.0]
solidUpperRight5 = [0.67,0.65,0.0]
solid5 = pyck.Cuboid(5,solidLowerLeft5,solidUpperRight5)
pack.AddShape(solid5)
solidLowerLeft6 = [0.55,solidUpperRight5[1],0.0]
solidUpperRight6 = [0.67,0.72,0.0]
solid6 = pyck.Cuboid(6,solidLowerLeft6,solidUpperRight6)
pack.AddShape(solid6)
solidLowerLeft7 = [0.55,solidUpperRight6[1],0.0]
solidUpperRight7 = [0.67,0.78,0.0]
solid7 = pyck.Cuboid(7,solidLowerLeft7,solidUpperRight7)
pack.AddShape(solid7)
solidLowerLeft8 = [0.55,solidUpperRight7[1],0.0]
solidUpperRight8 = [0.67,0.85,0.0]
solid8 = pyck.Cuboid(8,solidLowerLeft8,solidUpperRight8)
pack.AddShape(solid8)

funnel_top_left = pyck.Cuboid(9,[0.55,0.3,0.0],[0.57,0.9,0.0])
pack.AddShape(funnel_top_left)
#funnel_bot_left = pyck.Cuboid(2,[0.55,0.25,0.0],[0.59,0.3,0.0])
funnel_bot_left = pyck.TriPrism(9,[0.57,0.3,0.0],[0.60,0.3,0.0],[0.57,0.4,0.0],0.0)
pack.AddShape(funnel_bot_left)

funnel_top_right = pyck.Cuboid(9,[0.66,0.3,0.0],[0.69,0.90,0.0])
pack.AddShape(funnel_top_right)
#funnel_bot_right = pyck.Cuboid(2,[0.65,0.25,0.0],[0.69,0.30,0.0])
funnel_bot_right = pyck.TriPrism(9,[0.63,0.3,0.0],[0.67,0.3,0.0],[0.67,0.44,0.0],0.0)
pack.AddShape(funnel_bot_right)



pack.Process()
model = pyck.Model(pack)

# Set the material ID for all objects
matID = model.CreateIntField("material",1)
model.SetIntField(matID,1,0) 
model.SetIntField(matID,2,1)
model.SetIntField(matID,3,2) 
model.SetIntField(matID,4,3)
model.SetIntField(matID,5,4) 
model.SetIntField(matID,6,5)
model.SetIntField(matID,7,6) 
model.SetIntField(matID,8,7)
model.SetIntField(matID,9,8)
# Set the objectID's for all objects
objID = model.CreateIntField("objectID",1)
model.SetIntField(objID,1,0) 
model.SetIntField(objID,2,0) 
model.SetIntField(objID,3,0) 
model.SetIntField(objID,4,0) 
model.SetIntField(objID,5,0) 
model.SetIntField(objID,6,0) 
model.SetIntField(objID,7,0) 
model.SetIntField(objID,8,0) 
model.SetIntField(objID,9,1)

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
# 	idOb = matpoints.getObjID(p)
# 	if(idOb==0):
#             matpoints.setVel(p,0.0,-1.0,0.0)


# //// NODAL BOUNDARY CONDITIONS /////
for  n in range(nodes.numNodes):
   posX = nodes.getPosX(n)
   posY = nodes.getPosY(n)
   if ((posY >0.45) and (posX>(1.5*cellsize+0.57) and (posX<(0.65+2*psep))) ):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0,-0.50])
#    if ((posY >0.30) and ((posX<=0.54) or (posX>=0.68))):    
#        nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0.0,0.0])
   if (posY < 0.28):
       nodes.setBC(n,graphyt.BCTypes.grid_vx,[-0.50])
   if (posY < 0.20):
       nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[-0.50,0.0])
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")

# /////////// Parameters Set /////////////// #
tmax = 1.0  # //(time in seconds)
is3D = False
parameters = graphyt.Parameters(tmax,is3D)
parameters.setGravity(0,-10,0)
parameters.setDampingCoef(0.2)
parameters.setCFL(0.2)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")

# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
#Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
print("+++++ SOLVER  INFORMATION+++++ Starting Solver")
solver = graphyt.Solver(nodes, matpoints, parameters, materials)
print("+++++ SOLVER  INFORMATION+++++ Solver Ready")
matpntPos = matpoints.getPosArr()
matpntDisp = matpoints.getDispArr()
matpntMat = matpoints.getMatArr()
matpntCellID = matpoints.getCellIDArr()
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
vtp.AddArray("CellID", matpntCellID,1, 2)
vtp.AddArray("Mass", matpntMass,1, 2)
vtp.AddArray("Velocity",matpntVel,3, 2)
vtp.AddArray("Stress",matpntStress,6, 2)
vtp.AddArray("Normals",matpntNormal,3,2)
vtp.AddArray("AccPlasticStrain",matpntAccPlasStr,1,2)
vtp.AddArray("Damage",matpntDamage,1,2)
vtp.AddArray("Smax",matpntSMax,1,2)
# Iterate over nsteps number of timesteps
print("Beginning Timestepping")
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate_CPU(t, step, parameters)
    
    # //// NODAL BOUNDARY CONDITIONS /////
    if(step==int(0.5*int(solver.tmax/solver.dt)) ):
        for  n in range(nodes.numNodes):
            posY = nodes.getPosY(n)
            if (posY < 0.28):
                nodes.setBC(n,graphyt.BCTypes.grid_vx,[0.20])
            if (posY < 0.20):
                nodes.setBC(n,graphyt.BCTypes.grid_vx_vy,[0.20,0.0])


    #// After nsteps, output the data using serialize for nodes and material points.
    #matpoints.setRotatingBodyVel(1,0,0,0,0)
    vtpfname='output/extrusionTest'+str(step)+'.vtp'
    #vtifname='output/pouiselleVTI'+str(step)+'.vti'
    if (step % 1000 == 0):
        
        vtp.Write(vtpfname)

        #vti.Write(vtifname)

print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////







