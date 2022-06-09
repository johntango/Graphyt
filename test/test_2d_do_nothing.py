import unittest
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

class DoNothing2D(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Set the resolution and the pack
        cellsize = 0.0025
        psep = cellsize/2.0
        # Background Grid
        L = [0.2,0.2,0.0]
        nodes = graphyt.Nodes(L, cellsize)
        # Material Points
        matpoints = graphyt.MaterialPoints(nodes)
        matpoints.psep = psep
        # PYCK

        # create material
        rock = graphyt.MaterialModels()
        rock.setBulkMod(3.31e9)
        rock.setShearMod(2.48e9)
        rock.setDensity(1540.0)
        rock.setModel(1)

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
        pack.Process()
        model = pyck.Model(pack)
        # Set the material ID for all objects
        matID = model.CreateIntField("material",1)
        model.SetIntField(matID,1,0) 
        # Set the objectID's for all objects
        objID = model.CreateIntField("objectID",1)
        model.SetIntField(objID,1,1) 
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

        # /////////// Parameters Set /////////////// #
        tmax = 3.0e-5  # //(time in seconds)
        is3D = False
        parameters = graphyt.Parameters(tmax,is3D)
        parameters.setGravity(0, 0, 0)
        # // //+++++++++++++++++++++++++++++++++++++++++++++++//
        # // //+++++++++++++ MODEL RUN +++++++++++++++++++++++//
        # // //+++++++++++++++++++++++++++++++++++++++++++++++//
        #Now Mat points and nodes are set and ready to go, create a solver to run an MPM analysis on them
        solver = graphyt.Solver(nodes, matpoints, parameters, materials)
        # Iterate over nsteps number of timesteps
        for step in range(0,1500):  # ; t<solver.tmax; t+=solver.dt):
            t = solver.dt * step
            solver.iterate_CPU(t, step, parameters)
            cls.solver = solver
            cls.matpoints = matpoints
            cls.parameters = parameters
            cls.nodes = nodes
            super(DoNothing2D, cls).setUpClass()

    def test_Displacement(self):
            displacementX = self.matpoints.getDispX(5)
            testValue = 0.0
            self.assertAlmostEqual(displacementX,testValue)
    def test_Stress(self):
            stressXX = self.matpoints.getStressXX(5)
            testValue = 0.0
            self.assertAlmostEqual(stressXX,testValue)
    
    
