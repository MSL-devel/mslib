CCDEFAULT = g++ -O3 -msse2 -fopenmp
CCDEBUG = g++ -Wall -Wno-sign-compare -g 

CC = ${CCDEFAULT}

VPATH = src


SOURCE  = Atom AtomAngleRelationship AtomContainer AtomDihedralRelationship AtomDistanceRelationship \
          AtomGeometricRelationship AtomGroup AtomicPairwiseEnergy AtomSelection AtomVector CartesianGeometry \
          BBQTable BBQTableReader BBQTableWriter CartesianPoint CCD\
          Chain CharmmAngleInteraction CharmmBondInteraction CharmmDihedralInteraction \
          CharmmElectrostaticInteraction CharmmEnergy CharmmImproperInteraction CharmmParameterReader \
          CharmmSystemBuilder CharmmTopologyReader CharmmTopologyResidue CharmmUreyBradleyInteraction \
          CharmmVdwInteraction ChiStatistics CoiledCoils CrystalLattice DeadEndElimination EnergySet Enumerator EnvironmentDatabase \
          EnvironmentDescriptor File FourBodyInteraction Frame Helanal HelixFusion IcEntry IcTable Interaction \
          InterfaceResidueDescriptor Line LogicalParser MIDReader Matrix Minimizer MoleculeInterfaceDatabase \
          MslTools OptionParser PairwiseEnergyCalculator PDBFormat PDBReader PDBWriter PhiPsiReader PhiPsiStatistics PolymerSequence PSFReader \
          Position PotentialTable Predicate PrincipleComponentAnalysis PyMolVisualization Quaternion Quench Reader Residue ResiduePairTable \
          ResiduePairTableReader ResidueSubstitutionTable ResidueSubstitutionTableReader RotamerLibrary \
          RotamerLibraryReader SelfPairManager SphericalPoint SurfaceAreaAndVolume Symmetry System SystemRotamerLoader TBDReader ThreeBodyInteraction \
          Timer Transforms Tree TwoBodyDistanceDependentPotentialTable TwoBodyInteraction Writer UserDefinedInteraction  UserDefinedEnergy UserDefinedEnergySetBuilder 


HEADER = Hash.h MslExceptions.h Real.h Selectable.h Tree.h release.h 

TESTS   = testAtomGroup testAtomSelection testAtomVector testBBQ testBBQ2 testCCD testCharmmBuild testCharmmEnergies \
          testCharmmTopologyReader testCoiledCoils testDerivatives testEnergySet testEnvironmentDatabase \
          testEnvironmentDescriptor testFrame testGenerateCrystalLattice testHelixFusion testIcBuilding testLinkedPositions testLoopOverResidues \
          testMolecularInterfaceDatabase testMslToolsFunctions testPDBIO testPhiPsi testPolymerSequence testPSFReader testQuench \
          testResiduePairTable testResidueSubstitutionTable testSurfaceAreaAndVolume testSymmetry testSystemCopy \
          testSystemIcBuilding testTransforms testTree 


PROGRAMS = getSphericalCoordinates fillInSideChains generateCrystalLattice

GSL=T
GLPK=T
BOOST=T
32BIT=F


EXTERNAL_LIB_DIR=/usr/lib
#EXTERNAL_LIB_DIR=/library/sharedlibs64

ifeq ($(32BIT),T)
    EXTERNAL_LIB_DIR=/library/sharedlibs
endif

# GLPK Libraries
ifeq ($(GLPK),T)
    FLAGS          += -D__GLPK__
    SOURCE         += LinearProgrammingOptimization
    TESTS          += testRotamerOptimization
    PROGRAMS       += optimizeLP
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libglpk.a
endif


# GSL Libraries 
ifeq ($(GSL),T)

    FLAGS          += -D__GSL__
    SOURCE         += RandomNumberGenerator GSLMinimizer MonteCarloOptimization
    PROGRAMS       += optimizeMC
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libgsl.a ${EXTERNAL_LIB_DIR}/libgslcblas.a

endif

# BOOST Libraries
ifeq ($(BOOST),T)

    FLAGS          += -D__BOOST__ -DBOOST_DISABLE_THREADS
#    TESTS          += testBoost
#    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_serialization-gcc42-1_34_1.a
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_serialization-gcc43-mt-1_37.a


endif

# Include local Makefile
-include myProgs/myProgs.mk

# Generic Includes,Flags.  Static compile.
INCLUDE  = src -I/library/sharedincludes
FLAGS   += -static -DUSE_REAL_EQ_DOUBLE 


# Add proper suffix
OBJECTS       = $(patsubst %,objs/%.o, $(SOURCE)) 
MYOBJS        = $(patsubst %,objs/%.o, $(MYSOURCE)) 
BINARIES      = $(patsubst %,bin/%, $(PROGRAMS)) 
MYBINS        = $(patsubst %,bin/%, $(MYPROGS)) 
TESTBINS      = $(patsubst %,bin/%, $(TESTS)) 
MYHEADERFILES = $(patsubst, %,myProgs/%, $(MYHEADERS))
PHEADERS      = $(patsubst %,programs/%.h, $(PROGRAMS_HEADERS))


# Compile/Link commands
all: ${BINARIES} ${MYBINS} ${TESTBINS}

${OBJECTS}: objs/%.o : src/%.cpp src/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@  

${TESTBINS}: bin/% : tests/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}

${BINARIES}: bin/% : programs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}

${MYOBJS}: objs/%.o : myProgs/%.cpp myProgs/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@  

${MYBINS}: bin/% : myProgs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
	${CC} ${FLAGS} -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}

.PHONY : clean
clean :
	-rm -f ${OBJECTS} ${BINARIES} ${TESTBINS} ${MYOBJS} ${MYBINS}


python:
	gcc ${FLAGS} -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  
	g++ ${FLAGS} -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
	cp PythonMSL.so /usr/share/python-support/pymol/pymol/
