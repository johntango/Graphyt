import sys
#sys.path.insert(0,'/f/sjr/Desktop/Graphyt/build/
sys.path.insert(0,'/mnt/c/Users/sjr/Desktop/Graphyt/build/')
import graphyt

L = [1.0,5.0, 0.3]
cellsize = 0.025
nodes = graphyt.Nodes(L, cellsize)
matpoints = graphyt.MaterialPoints(nodes)
psep = 0.5*cellsize
matpoints.psep = psep

# Define Drop
G = 0.0
K = 1.0e6
rho1=1100
material = 1
xmin = cellsize
xmax = 1.0-cellsize
ymin = 2.50
ymax = 5-2*cellsize
# 2D Version
zmin = 0.1#0.2
zmax = 0.1#0.8
x = xmin
y = ymin
z = zmin
i=0
while z <= zmax:
    y = ymin
    while y <= ymax:
        x = xmin
        while x <= xmax:
            # // Create Material Po with Initialized quantities set to zero.
            # Assuming we now have a state for this particular particle, add it.
            # state - [x,y,z,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,exx,eyy,ezz,exy,eyz,exz,
            # ...eplxx,eplyy,eplzz,eplxy,eplyz,eplxz,material]
            # // 31 property entries
            matpoints.add(i++,x, y, z, material,G,K,rho1) 
            # // adds particle state to the materialpoints array
            x += psep
        y += psep
    z += psep




# Define Water
G = 0.0
K = 1.0e6
rho2=1000.0
material = 1
xmin = cellsize
xmax = 1.0-cellsize
ymin = cellsize
ymax = 2.50
# 2D Version
zmin = 0.1#0.2
zmax = 0.1#0.8
x = xmin
y = ymin
z = zmin
while z <= zmax:
    y = ymin
    while y <= ymax:
        x = xmin
        while x <= xmax:
            # // Create Material Po with Initialized quantities set to zero.
            # Assuming we now have a state for this particular particle, add it.
            # state - [x,y,z,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,exx,eyy,ezz,exy,eyz,exz,
            # ...eplxx,eplyy,eplzz,eplxy,eplyz,eplxz,material]
            # // 31 property entries
            matpoints.add(i++,x, y, z, material,G,K,rho2) 
            # // adds particle state to the materialpoints array
            x += psep
        y += psep
    z += psep



print("+++++ MAT POINT INFO+++++ No. Particles:", matpoints.numParticles)
# // Now that nodes and material points have been discretised, boundary conditions are applied
# /////////-----------BOUNDARY CONDITIONS--------------///////
# //// MATERIAL POINT BOUNDARY CONDITIONS /////


# //// NODAL BOUNDARY CONDITIONS /////
#     //// NODAL BOUNDARY CONDITIONS /////
#// Wall boundaries for Grid Nodes
wallthickness = 2.0*nodes.cellsize
for  n in range(nodes.numNodes):
    posX = nodes.getPosX(n)
    posY = nodes.getPosY(n)
    posZ = nodes.getPosZ(n)
    if (((posX> (L[0]-wallthickness)) or (posX < wallthickness )) or ((posY> (L[1]-wallthickness)) or (posY < wallthickness)) ):
        #  // Want a format like this: if(isInside(nodes,n,walls))
        nodes.setFixedBC(n)
print("+++++ MODEL  INFORMATION+++++ Boundary Conditions Applied")
# /////////// Parameters Set /////////////// #
tmax = 10.0  # //(time in seconds)
int_scheme = 1
materialModel = 0  # // 0 - Elastic 1 - EOS Pure Fluid
G = 0.0  # //48.0e5;//0.0;//0.3e9;//0.0;//0.0;//48e9;//0.0;//48e9;//0.0;//48e9;;
K = 1.0e6  # //140.0e5;//22e4;//1.6e9;//22e9;//140e9;//22e4;//140e9;//22e5;//140e9;
is3D = False
parameters = graphyt.Parameters(tmax, int_scheme, materialModel,is3D, G, K)
parameters.setGravity(0, -10, 0)
print("+++++ MODEL  INFORMATION+++++ Parameters Set")

# // //+++++++++++++++++++++++++++++++++++++++++++++++//
# // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
# // //+++++++++++++++++++++++++++++++++++++++++++++++//
#Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
print("+++++ SOLVER  INFORMATION+++++ Starting Solver")
solver = graphyt.Solver(nodes, matpoints, parameters)

#0//     // Iterate over nsteps number of timesteps
for step in range(int(solver.tmax/solver.dt)):  # ; t<solver.tmax; t+=solver.dt):
    t = solver.dt * step
    solver.iterate(t, step, parameters)
    #//     // After nsteps, output the data using serialize for nodes and material points.
    filename='taylorinst2D'
    if step % 2500 == 0:
        matpoints.Serialise(step,filename)
       # nodes.Serialise(step,filename)
#endTimeLoop =omp_get_wtime()#; // Get finishing time of total timeloop
#totSimTime = endTimeLoop - solver.startTimeLoop;
print("Total Simulation Duration (s)")
# // ///////////////////////////////////////////////////
