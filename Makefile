CCDEFAULT = g++ -O3 -msse2 -fopenmp
CCDEBUG = g++ -Wall -Wno-sign-compare -g 

GSLDEFAULT = T
GLPKDEFAULT = T
BOOSTDEFAULT = T
ARCH32BITDEFAULT = F

EXTERNAL_LIB_DIR_DEFAULT=/usr/lib

CC = ${CCDEFAULT}

VPATH = src


SOURCE  = Atom Atom3DGrid AtomAngleRelationship AtomContainer AtomDihedralRelationship AtomDistanceRelationship \
          AtomGeometricRelationship AtomGroup AtomicPairwiseEnergy AtomSelection AtomVector CartesianGeometry \
          BBQTable BBQTableReader BBQTableWriter CartesianPoint\
          Chain CharmmAngleInteraction CharmmBondInteraction CharmmDihedralInteraction \
          CharmmElectrostaticInteraction CharmmEnergy CharmmImproperInteraction CharmmParameterReader \
          CharmmSystemBuilder CharmmTopologyReader CharmmTopologyResidue CharmmUreyBradleyInteraction \
          CharmmVdwInteraction ChiStatistics CoiledCoils CrystalLattice DeadEndElimination EnergySet EnergeticAnalysis Enumerator EnvironmentDatabase \
          EnvironmentDescriptor File FourBodyInteraction Frame Helanal HelixFusion IcEntry IcTable Interaction \
          InterfaceResidueDescriptor Line LogicalParser MIDReader Matrix Minimizer MoleculeInterfaceDatabase \
          MslTools OptionParser PairwiseEnergyCalculator PDBFormat PDBFragments PDBReader PDBWriter PhiPsiReader PhiPsiStatistics PolymerSequence PSFReader \
          Position PotentialTable Predicate PrincipleComponentAnalysis PyMolVisualization Quaternion Reader Residue ResiduePairTable \
          ResiduePairTableReader ResidueSelection ResidueSubstitutionTable ResidueSubstitutionTableReader RotamerLibrary \
          RotamerLibraryReader SelfPairManager SasaAtom SasaCalculator SphericalPoint SurfaceSphere Symmetry System SystemRotamerLoader TBDReader \
          ThreeBodyInteraction Timer Transforms Tree TwoBodyDistanceDependentPotentialTable TwoBodyInteraction Writer UserDefinedInteraction  UserDefinedEnergy \
          UserDefinedEnergySetBuilder HelixGenerator


HEADER = Hash.h MslExceptions.h Real.h Selectable.h Tree.h release.h 

TESTS   = testAtomGroup testAtomSelection testAtomVector testBackRub testBBQ testBBQ2 testCCD testCharmmBuild testCharmmEnergies \
          testCharmmTopologyReader testCoiledCoils testDerivatives testEnergySet testEnergeticAnalysis testEnvironmentDatabase \
          testEnvironmentDescriptor testFrame testGenerateCrystalLattice testHelixFusion testIcBuilding testLinkedPositions testLoopOverResidues \
          testMolecularInterfaceDatabase testMslToolsFunctions testPDBIO testPDBFragments testPhiPsi testPolymerSequence testPSFReader testQuench \
          testRegEx testResiduePairTable testResidueSubstitutionTable testSurfaceAreaAndVolume testSymmetry testSystemCopy \
          testSystemIcBuilding testTransforms testTree testHelixGenerator



PROGRAMS = getSphericalCoordinates fillInSideChains generateCrystalLattice createFragmentDatabase getDihedrals energyTable analEnergy grepSequence \
           getSelection alignMolecules calculateSasa runQuench runKBQuench searchFragmentDatabase tableEnergies optimizeLP.cpp optimizeMC.cpp


# To ever-ride the defaults, set the GSL GLPK BOOST
ifndef GSL
   GSL=${GSLDEFAULT}
endif
ifndef GLPK
   GLPK=${GLPKDEFAULT}
endif
ifndef BOOST
   BOOST=${BOOSTDEFAULT}
endif
ifndef ARCH32BIT
   ARCH32BIT=${ARCH32BITDEFAULT}
endif

ifndef EXTERNAL_LIB_DIR
   EXTERNAL_LIB_DIR=${EXTERNAL_LIB_DIR_DEFAULT}
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
    SOURCE         += RandomNumberGenerator GSLMinimizer MonteCarloOptimization CCD BackRub Quench SurfaceAreaAndVolume
    PROGRAMS       += optimizeMC
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libgsl.a ${EXTERNAL_LIB_DIR}/libgslcblas.a

endif

# BOOST Libraries
ifeq ($(BOOST),T)

    FLAGS          += -D__BOOST__ -DBOOST_DISABLE_THREADS
    SOURCE         +=  RegEx
#    TESTS          += testBoost
#   NOTE   CHANGE!!! THE FOLLOWING SHOULD NOT HAVE A VERSION, IF USING SPECIAL LOCATIONS FOR LIBRARIES USE SYMLINKS TO POINT TO THE DESIRED ONE
#    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_serialization-gcc43-mt-1_37.a ${EXTERNAL_LIB_DIR}/libboost_regex-mt.a
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_serialization.a ${EXTERNAL_LIB_DIR}/libboost_regex-mt.a
endif



# Generic Includes,Flags.  Static compile.  
# NOTE IS THE FOLLOWING STILL NECESSARY?
INCLUDE  = src
ifdef CUSTOMINCLUDES
   INCLUDE += -I${CUSTOMINCLUDES}
endif
FLAGS   += -static -DUSE_REAL_EQ_DOUBLE 


# Include local Makefile
-include myProgs/myProgs.mk
# -include Makefile.local

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
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${BINARIES}: bin/% : programs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${MYOBJS}: objs/%.o : myProgs/%.cpp myProgs/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@  

${MYBINS}: bin/% : myProgs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
	${CC} ${FLAGS} -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}  -lpthread

.PHONY : clean
clean :
	-rm -f ${OBJECTS} ${BINARIES} ${TESTBINS} ${MYOBJS} ${MYBINS}


python:
	gcc ${FLAGS} -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  
	g++ ${FLAGS} -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
	cp PythonMSL.so /usr/share/python-support/pymol/pymol/
