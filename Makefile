CCOPTIM = g++ -O3 -msse3 -mfpmath=sse -funroll-loops  
CCDEBUG = g++ -Wall -msse3 -mfpmath=sse -funroll-loops -Wno-sign-compare -g


GSLDEFAULT = F
GLPKDEFAULT = F
BOOSTDEFAULT = F
ARCH32BITDEFAULT = F
FFTWDEFAULT = F
MACOSDEFAULT = F
DEBUGDEFAULT = F

EXTERNAL_LIB_DIR_DEFAULT=/usr/lib


VPATH = src


SOURCE  = ALNReader Atom Atom3DGrid AtomAngleRelationship AtomContainer AtomDihedralRelationship AtomDistanceRelationship \
          AtomGeometricRelationship AtomGroup AtomicPairwiseEnergy AtomSelection AtomPointerVector CartesianGeometry \
          BBQTable BBQTableReader BBQTableWriter CartesianPoint\
          Chain CharmmAngleInteraction CharmmBondInteraction CharmmDihedralInteraction \
          CharmmElectrostaticInteraction CharmmEnergy CharmmImproperInteraction CharmmParameterReader CharmmEEF1ParameterReader \
          CharmmSystemBuilder CharmmTopologyReader CharmmTopologyResidue CharmmUreyBradleyInteraction \
          CharmmVdwInteraction CharmmEEF1Interaction CharmmEEF1RefInteraction ChiStatistics CoiledCoils CrystalLattice DeadEndElimination EnergySet EnergeticAnalysis Enumerator EnvironmentDatabase \
          EnvironmentDescriptor File FourBodyInteraction Frame Helanal HelixFusion IcEntry IcTable Interaction \
          InterfaceResidueDescriptor Line LogicalParser MIDReader Matrix Minimizer MoleculeInterfaceDatabase \
          MslTools OptionParser PairwiseEnergyCalculator PDBFormat PDBFragments PDBReader PDBWriter PhiPsiReader PhiPsiStatistics PolymerSequence PSFReader \
          Position PotentialTable Predicate PrincipleComponentAnalysis PyMolVisualization Quaternion Reader Residue ResiduePairTable \
          ResiduePairTableReader ResidueSelection ResidueSubstitutionTable ResidueSubstitutionTableReader RotamerLibrary \
          RotamerLibraryReader SelfPairManager SasaAtom SasaCalculator SphericalPoint SurfaceSphere Symmetry System SystemRotamerLoader TBDReader \
          ThreeBodyInteraction Timer Transforms Tree TwoBodyDistanceDependentPotentialTable OneBodyInteraction TwoBodyInteraction Writer UserDefinedInteraction  UserDefinedEnergy \
          UserDefinedEnergySetBuilder HelixGenerator RotamerLibraryBuilder RotamerLibraryWriter AtomBondBuilder


HEADER = Hash.h MslExceptions.h Real.h Selectable.h Tree.h release.h 

TESTS   = testAtomGroup testAtomSelection testAtomPointerVector testBackRub testBBQ testBBQ2 testCCD testCharmmBuild testCharmmEnergies \
          testCharmmTopologyReader testCoiledCoils testDerivatives testEnergySet testEnergeticAnalysis testEnvironmentDatabase \
          testEnvironmentDescriptor testFrame testGenerateCrystalLattice testHelixFusion testIcBuilding testLinkedPositions testLoopOverResidues \
          testMolecularInterfaceDatabase testMslToolsFunctions testPDBIO testPDBFragments testPhiPsi testPolymerSequence testPSFReader testQuench \
          testResiduePairTable testResidueSubstitutionTable testSasaCalculator testSurfaceAreaAndVolume testSymmetry testSystemCopy \
          testSystemIcBuilding testTransforms testTree testHelixGenerator testRandomSeqGenerator testRotamerLibraryWriter testNonBondedCutoff  testALNReader \
	  testAtomAndResidueId testAtomBondBuilder testTransformBondAngleDiheEdits testAtomContainer testCharmmEEF1ParameterReader testEEF1 testEEF1_2





PROGRAMS = getSphericalCoordinates fillInSideChains generateCrystalLattice createFragmentDatabase getDihedrals energyTable analEnergy grepSequence \
	   getSelection alignMolecules calculateSasa runQuench runKBQuench searchFragmentDatabase tableEnergies printSequence generateCoiledCoils getSurroundingResidues



# To ever-ride the defaults, set the GSL GLPK BOOST MSLDEBUG environmental variables
ifndef MSL_DEBUG
   MSL_DEBUG = DEBUGDEFAULT
endif
ifndef MSL_GSL
   MSL_GSL=${GSLDEFAULT}
endif
ifndef MSL_GLPK
   MSL_GLPK=${GLPKDEFAULT}
endif
ifndef MSL_BOOST
   MSL_BOOST=${BOOSTDEFAULT}
endif
ifndef ARCH32BIT
   ARCH32BIT=${ARCH32BITDEFAULT}
endif
ifndef FFTW
   FFTW=${FFTWDEFAULT}
endif

ifndef EXTERNAL_LIB_DIR
   EXTERNAL_LIB_DIR=${EXTERNAL_LIB_DIR_DEFAULT}
endif

ifndef MACOS
    MACOS=${MACOSDEFAULT}
endif

ifeq ($(MSL_DEBUG),T)
   CC = ${CCDEBUG}
else
   CC= ${CCOPTIM}
endif

# GLPK Libraries
ifeq ($(MSL_GLPK),T)
    FLAGS          += -D__GLPK__
    SOURCE         += LinearProgrammingOptimization
    TESTS          += testRotamerOptimization
    PROGRAMS       += optimizeLP
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libglpk.a
endif


# GSL Libraries 
ifeq ($(MSL_GSL),T)

    FLAGS          += -D__GSL__
    SOURCE         += RandomNumberGenerator GSLMinimizer MonteCarloOptimization CCD BackRub Quench SurfaceAreaAndVolume
    PROGRAMS       += optimizeMC
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libgsl.a ${EXTERNAL_LIB_DIR}/libgslcblas.a

endif

# BOOST Libraries
ifeq ($(MSL_BOOST),T)

    FLAGS          += -D__BOOST__ -DBOOST_DISABLE_THREADS
    SOURCE         +=  RegEx RandomSeqGenerator
    TESTS          += testRegEx
#    TESTS          += testBoost
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_serialization.a 

# For MAC I only compiled the non-multithreaded library, sometime I'll figure it out, but we do not use multi-threading so for now this is ok.
ifeq ($(MACOS),T)
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libboost_regex.a
else
    STATIC_LIBS    +=  ${EXTERNAL_LIB_DIR}/libboost_regex-mt.a
endif


endif

ifeq ($(FFTW),T)
    STATIC_LIBS    += ${EXTERNAL_LIB_DIR}/libfftw3.a
endif



# Generic Includes,Flags.  Static compile.  
# NOTE IS THE FOLLOWING STILL NECESSARY?
INCLUDE  = src
ifdef CUSTOMINCLUDES
   INCLUDE += -I${CUSTOMINCLUDES}
endif


# Add a MACOS flag for certain code breaks (see bottom of Tree.h, Selectable.h ... templated classes don't need pre-instantiations?)
ifeq ($(MACOS),T)
    FLAGS += -D__MACOS__ -DUSE_REAL_EQ_DOUBLE 
else
    FLAGS   += -static -DUSE_REAL_EQ_DOUBLE 
endif

# Include local Makefile
-include myProgs/myProgs.mk
# -include Makefile.local

# Include local Makefile
-include examples/examples.mk
# -include Makefile.local

# Add proper suffix
OBJECTS       = $(patsubst %,objs/%.o, $(SOURCE)) 
MYOBJS        = $(patsubst %,objs/%.o, $(MYSOURCE)) 
BINARIES      = $(patsubst %,bin/%, $(PROGRAMS)) 
EXAMPLEBINS   = $(patsubst %,bin/%, $(EXAMPLES)) 
MYBINS        = $(patsubst %,bin/%, $(MYPROGS)) 
TESTBINS      = $(patsubst %,bin/%, $(TESTS)) 
MYHEADERFILES = $(patsubst, %,myProgs/%, $(MYHEADERS))
PHEADERS      = $(patsubst %,programs/%.h, $(PROGRAMS_HEADERS))




# Compile/Link commands
all: ${BINARIES} ${MYBINS} ${TESTBINS} ${EXAMPLEBINS}
programs: ${BINARIES}
tests: ${TESTBINS}
examples: ${EXAMPLEBINS}
mybins: ${MYBINS}

${OBJECTS}: objs/%.o : src/%.cpp src/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@  

${TESTBINS}: bin/% : tests/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${BINARIES}: bin/% : programs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${EXAMPLEBINS}: bin/% : examples/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${MYOBJS}: objs/%.o : myProgs/%.cpp myProgs/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@  

${MYBINS}: bin/% : myProgs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
	${CC} ${FLAGS} -I${INCLUDE} -o $@  ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}  -lpthread

.PHONY : clean
clean :
	-rm -f ${OBJECTS} ${BINARIES} ${EXAMPLEBINS} ${TESTBINS} ${MYOBJS} ${MYBINS}


pythonLin:
	gcc ${FLAGS} -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/Library/Frameworks/Python.framework/Versions/2.6/Headers/
	g++ ${FLAGS} -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
	cp PythonMSL.so /usr/share/python-support/pymol/pymol/


pythonMac:
	gcc ${FLAGS} -fPIC -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/System/Library/Frameworks/Python.framework/Versions/2.5/Headers/
	g++ ${FLAGS} -bundle -undefined dynamic_lookup objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
#       sudo cp PythonMSL.so /Applications/PyMOLX11Hybrid.app/pymol/modules
