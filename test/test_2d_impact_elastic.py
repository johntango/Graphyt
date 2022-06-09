import unittest
# Import the necessary modules
import sys
# sys.path.insert(0,'/f/sjr/Graphyt/build/')
import graphyt
sys.path.insert(0,'/f/sjr/pyck/build/swig/')
import pyck
sys.path.insert(0,'/f/sjr/VTPWriter/build/')
import math

class DiskImpact2D(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Import the necessary modules

        ################################################################
        #--------------------------------------------------------------#
        #                        GEOMETRY                              #
        #--------------------------------------------------------------#
        ################################################################

        ### GRAPHYT - Initialising nodes and material points
        # Grid cellsize
        cellsize = 0.0005
        # Initial particle separation
        psep = cellsize/2
        # Compuational Domain
        L = [0.10,0.10,0.0]
        # nodes and matpoint objects
        nodes = graphyt.Nodes(L,cellsize)
        matpoints = graphyt.MaterialPoints(nodes)
        matpoints.psep = psep
        ### PYCK - Particle Geometry

        cubic = pyck.CubicPacker(L,psep,[0.25*cellsize,0.25*cellsize,0.0])
        pack = pyck.StructuredPack(cubic)

        # square1_lowerLeft = [0.045,0.02,0.0]
        # square1_upperRight = [0.055,0.05,0.0] 
        # square1 = pyck.Cuboid(1,square1_lowerLeft,square1_upperRight)
        disk1 = pyck.Sphere(1,[0.050,0.03,0.0],0.005)
        # square2_lowerLeft = [0.045,square1_upperRight[1]+3*cellsize,0.0]
        # square2_upperRight = [0.055,square2_lowerLeft[1]+0.03,0.0] 
        # square2 = pyck.Cuboid(2,square2_lowerLeft,square2_upperRight)
        disk2 = pyck.Sphere(2,[0.050,0.03+0.01 + 3*cellsize,0.0],0.005)
        ### PYCK - Add all shapes to pack
        pack.AddShape(disk1)
        pack.AddShape(disk2)
        pack.Process()
        ### PYCK -> GRAPHYT - add pyck geometry to graphyt matpoints
        numParticles = pack.GetNumParticles()
        pyckPos = pack.GetPositions()
        matpoints.addPyckGeom(pyckPos,numParticles)
        ################################################################
        #--------------------------------------------------------------#
        #                        Materials                             #
        #--------------------------------------------------------------#
        ################################################################

        ### PYCK - Materials
        mat1 = graphyt.MaterialModels()
        mat1.setBulkMod(1000)
        mat1.setShearMod(500)
        mat1.setDensity(10)
        mat2 = graphyt.MaterialModels()
        mat2.setBulkMod(1000)
        mat2.setShearMod(500)
        mat2.setDensity(10)
        materials = [mat1, mat2]
        ### PYCK - Create Model for Simulation
        model = pyck.Model(pack)
        matID = model.CreateIntField("material",1)
        objID = model.CreateIntField("objectID",1)
        # Assign material IDs to each object
        model.SetIntField(matID,1,0)
        model.SetIntField(matID,2,1)
        matIDfield = model.GetIntField(matID)
        matpoints.addPyckIntField(matIDfield)
        # Assign object IDs to each object
        model.SetIntField(objID,1,1)
        model.SetIntField(objID,2,2)
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
            if(obID == 1):
                matpoints.setVel(p,0.0,1.50,0)
            if(obID == 2):
                matpoints.setVel(p,0.0,-1.50,0)
        ### Nodal BCs

        ################################################################
        #--------------------------------------------------------------#
        #                   Simulation Parameters                      #
        #--------------------------------------------------------------#
        ################################################################

        tmax = 5e-3
        is3D = False
        parameters = graphyt.Parameters(tmax, is3D)
        # parameters.setDampingCoef(1.0)
        parameters.setCFL(0.1)
        ################################################################
        #--------------------------------------------------------------#
        #                   Simulation Solver                          #
        #--------------------------------------------------------------#
        ################################################################

        ## Create the solver object
        solver = graphyt.Solver(nodes, matpoints, parameters, materials)

        ################################################################
        #--------------------------------------------------------------#
        #                   Time Integration                           #
        #--------------------------------------------------------------#
        ################################################################
        for step in range(int(solver.tmax/solver.dt)):
            t = solver.dt * step
            solver.iterate_CPU(t, step, parameters)
            cls.solver = solver
            cls.matpoints = matpoints
            cls.parameters = parameters
            cls.nodes = nodes
            super(DiskImpact2D, cls).setUpClass()

    def test_Displacement(self):
            displacementY = self.matpoints.getDispY(640)
            testValue = -0.0024606676268147457
            self.assertAlmostEqual(displacementY,testValue)
    def test_Stress(self):
            stressYY = self.matpoints.getStressYY(640)
            testValue = 4.943548848378589
            self.assertAlmostEqual(stressYY,testValue)
    
    
