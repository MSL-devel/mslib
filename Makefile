############################################################################################
#
#  Enviromental variables that control for the presence of external libraries and
#  debug level
#
#  $MSL_GSL    set to "T" if GSL is installed or else to "F" (default)
#  $MSL_BOOST  set to "T" if BOOST is installed or else to "F" (default)
#  $MSL_GLPK   set to "T" if GLPK is installed or else to "F" (default)
#  $MSL_DEBUG  set to "T" to compile in "debug" mode, or else to "F" (default)
# 
#  Set in your .cshrc
#    setenv MSL_GSL T
#    setenv MSL_GLPK T
#    setenv MSL_BOOST T
#    setenv MSL_DEBUG F
#
#  Set in your .bash
#    export MSL_GSL=T
#    export MSL_GLPK=T
#    export MSL_BOOST=T
#    export MSL_DEBUG=F
############################################################################################

# Define compiler command
#CC  = g++ 
CCOPTIM = g++ -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp
CCDEBUG = g++ -Wall -Wno-sign-compare -msse3 -mfpmath=sse -funroll-loops -g


GSLDEFAULT = F
GSLOLDDEFAULT = F
GLPKDEFAULT = F
BOOSTDEFAULT = F
ARCH32BITDEFAULT = F
FFTWDEFAULT = F
RDEFAULT = F
MACOSDEFAULT = F
DEBUGDEFAULT = F
TESTINGDEFAULT = F
MSLOUT_DEBUG_OFFDEFAULT = T
STATICDEFAULT = T
EXTERNAL_LIB_DIR_DEFAULT=/usr/lib
EXTERNAL_INCLUDE_DIR_DEFAULT=/usr/include

VPATH = src


SOURCE  = ALNReader Atom Atom3DGrid AtomAngleRelationship AtomContainer AtomDihedralRelationship AtomDistanceRelationship \
          AtomGeometricRelationship AtomGroup AtomicPairwiseEnergy AtomSelection AtomPointerVector CartesianGeometry \
          BaselineEnergyBuilder BaselineInteraction BBQTable BBQTableReader BBQTableWriter CartesianPoint\
          Chain CharmmAngleInteraction CharmmBondInteraction CharmmDihedralInteraction \
          CharmmElectrostaticInteraction CharmmEnergy CharmmImproperInteraction CharmmParameterReader CharmmEEF1ParameterReader \
          CharmmSystemBuilder CharmmTopologyReader CharmmTopologyResidue CharmmUreyBradleyInteraction \
          CharmmVdwInteraction CharmmEEF1Interaction CharmmEEF1RefInteraction ChiStatistics CoiledCoils CrystalLattice DeadEndElimination EnergySet EnergeticAnalysis Enumerator EnvironmentDatabase \
          EnvironmentDescriptor File FourBodyInteraction Frame FuseChains Helanal HelixFusion HydrogenBondBuilder IcEntry IcTable Interaction \
          InterfaceResidueDescriptor Line LogicalParser MIDReader Matrix Minimizer MoleculeInterfaceDatabase \
          MslOut MslTools OptionParser PairwiseEnergyCalculator CRDFormat PDBFormat PDBFragments PDBReader PDBWriter PDBTopology CRDReader CRDWriter PolymerSequence PSFReader \
          Position PotentialTable Predicate PrincipleComponentAnalysis PyMolVisualization Quaternion Reader Residue ResiduePairTable \
          ResiduePairTableReader ResidueSelection ResidueSubstitutionTable ResidueSubstitutionTableReader RotamerLibrary \
          RotamerLibraryReader SelfPairManager SasaAtom SasaCalculator Scwrl4HBondInteraction SphericalPoint SurfaceSphere Symmetry System SystemRotamerLoader TBDReader \
          ThreeBodyInteraction Timer Transforms Tree TwoBodyDistanceDependentPotentialTable OneBodyInteraction TwoBodyInteraction Writer UserDefinedInteraction  UserDefinedEnergy \
          UserDefinedEnergySetBuilder HelixGenerator RotamerLibraryBuilder RotamerLibraryWriter AtomBondBuilder LogicalCondition MonteCarloManager \
	  SelfConsistentMeanField PhiPsiReader PhiPsiStatistics RandomNumberGenerator \
	  BackRub CCD MonteCarloOptimization Quench SpringConstraintInteraction SurfaceAreaAndVolume VectorPair VectorHashing PDBTopologyBuilder SysEnv


HEADER = Hash.h MslExceptions.h Real.h Selectable.h Tree.h release.h 

TESTS   = testAtomGroup testAtomSelection testAtomPointerVector testBBQ testBBQ2 testCharmmBuild testCharmmEnergies \
          testCharmmTopologyReader testCoiledCoils testEnergySet testEnergeticAnalysis testEnvironmentDatabase \
          testEnvironmentDescriptor testFrame testGenerateCrystalLattice testHelixFusion testIcBuilding testLinkedPositions testLoopOverResidues \
          testMolecularInterfaceDatabase testMslToolsFunctions testCRDIO testPDBIO testPDBFragments testPhiPsi testPolymerSequence testPSFReader \
          testResiduePairTable testResidueSubstitutionTable testSasaCalculator testSymmetry testSystemCopy \
          testSystemIcBuilding testTransforms testTree testHelixGenerator testRotamerLibraryWriter testNonBondedCutoff  testALNReader \
	  testAtomAndResidueId testAtomBondBuilder testTransformBondAngleDiheEdits testAtomContainer testCharmmEEF1ParameterReader testEEF1 testEEF1_2 \
	  testResidueSelection testAddCharmmIdentity testMslOut testMslOut2 testRandomNumberGenerator \
	  testPDBTopology testVectorPair testSharedPointers2 testTokenize testMinimization testSaveAtomAltCoor testPDBTopologyBuild testSysEnv testBoost

PROGRAMS = getSphericalCoordinates fillInSideChains generateCrystalLattice createFragmentDatabase getDihedrals energyTable analEnergy \
	   getSelection alignMolecules calculateSasa searchFragmentDatabase printSequence getSurroundingResidues \
           insertLoopIntoTemplate setConformation coiledCoilBuilder findClashes mutate calculateDistanceOrAngle minimize repackSideChains backrubPdb

# PROGRAMS/TESTS_THAT_DO_NOT_COMPLILE =  generateCoiledCoils testBoost testRInterface 
# PROGRAMS/TESTS_WITHOUT_SOURCE_FILES =  testBoostSpriit testBoostSpirit2 testLogicalCondition testCoiledCoil testDistanceHashing 


# To ever-ride the defaults, set the $MSL_GSL $MSL_GLPK $MSL_BOOST $MSL_DEBUG and $MSL_TESTING environmental variables
ifndef MSL_DEBUG
   MSL_DEBUG = DEBUGDEFAULT
endif
ifndef MSL_TESTING
   MSL_TESTING = TESTINGDEFAULT
endif
ifndef MSL_GSL
   MSL_GSL=${GSLDEFAULT}
endif
ifndef MSL_GSL_OLD
   MSL_GSL_OLD=${GSLOLDDEFAULT}
endif
ifndef MSL_GLPK
   MSL_GLPK=${GLPKDEFAULT}
endif
ifndef MSL_BOOST
   MSL_BOOST=${BOOSTDEFAULT}
endif
ifndef MSL_STATIC
   MSL_STATIC=${STATICDEFAULT}
endif

ifndef ARCH32BIT
   ARCH32BIT=${ARCH32BITDEFAULT}
endif
ifndef FFTW
   FFTW=${FFTWDEFAULT}
endif

ifndef MSL_R
    MSL_R=${RDEFAULT}
endif

ifndef MSL_EXTERNAL_LIB_DIR
   MSL_EXTERNAL_LIB_DIR=${EXTERNAL_LIB_DIR_DEFAULT}
endif

ifndef MSL_EXTERNAL_INCLUDE_DIR
   MSL_EXTERNAL_INCLUDE_DIR=${EXTERNAL_INCLUDE_DIR_DEFAULT}
endif

ifndef MSL_MACOS
    MSL_MACOS=${MACOSDEFAULT}
endif

ifndef MSL_MSLOUT_DEBUG_OFF
   MSL_MSLOUT_DEBUG_OFF=${MSLOUT_DEBUG_OFFDEFAULT}
endif

ifeq ($(MSL_DEBUG),T)
    CC = ${CCDEBUG}
#   FLAGS =   -Wall -Wno-sign-compare -msse3 -mfpmath=sse -funroll-loops -g 
#   LINKFLAGS =
else
    CC= ${CCOPTIM}
#   FLAGS =  -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp
#   LINKFLAGS =
endif


# compilation with alternative testing code
ifeq ($(MSL_TESTING),T)
    FLAGS          += -D__TESTING__
    SOURCE         += 
    TESTS          += 
    PROGRAMS       += 
    STATIC_LIBS    += 
endif

# GLPK Libraries
ifeq ($(MSL_GLPK),T)
    FLAGS          += -D__GLPK__
    SOURCE         += LinearProgrammingOptimization
    TESTS          += testRotamerOptimization
    PROGRAMS       += optimizeLP
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libglpk.a
endif


# GSL Libraries 
ifeq ($(MSL_GSL),T)
    FLAGS          += -D__GSL__
    SOURCE         += GSLMinimizer
    TESTS          += testQuench testDerivatives testCCD testBackRub testSurfaceAreaAndVolume
    PROGRAMS       += tableEnergies runQuench runKBQuench optimizeMC
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libgsl.a ${MSL_EXTERNAL_LIB_DIR}/libgslcblas.a
endif
# GSL OLD remove some functions in the Minimizer when the GSL version available is pre-v1.14
ifeq ($(MSL_GSL_OLD),T)
    FLAGS          += -D__GSL_OLD__
endif

# BOOST Libraries
ifeq ($(MSL_BOOST),T)
    FLAGS          += -D__BOOST__ -DBOOST_DISABLE_THREADS
    SOURCE         +=  RegEx RandomSeqGenerator
    TESTS          += testRegEx testRandomSeqGenerator
    PROGRAMS       +=  grepSequence
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_serialization.a 
    # For MAC I only compiled the non-multithreaded library, sometime I'll figure it out, but we do not use multi-threading so for now this is ok.
    ifeq ($(MSL_MACOS),T)
        STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_regex.a
    else
        STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libboost_regex.a
    endif
endif

ifeq ($(FFTW),T)
    STATIC_LIBS    += ${MSL_EXTERNAL_LIB_DIR}/libfftw3.a
endif

# R Libraries
ifeq ($(MSL_R),T)
     R_HOME         = /opt/applications/R/2.12.1/gnu/
     FLAGS         += -D__R__ 

ifeq ($(MSL_STATIC),T)     
     STATIC_LIBS   += ${MSL_EXTERNAL_LIB_DIR}/libRcpp.a ${MSL_EXTERNAL_LIB_DIR}/libRInside.a
else

ifeq ($(MSL_MACOS),T)
     LINKFLAGS     += -framework R
else
     RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags | grep -v WARN)
     RLDFLAGS  := 		$(shell $(R_HOME)/bin/R CMD config --ldflags | grep -v WARN)
     RBLAS     :=	        $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS | grep -v WARN)
     RLAPACK   := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS | grep -v WARN)
     #
     ### if you need to set an rpath to R itself, also uncomment
     RRPATH :=		-Wl,-rpath,$(R_HOME)/lib
     #
     ### include headers and libraries for Rcpp interface classes
     RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave | grep -v WARN)
     RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave | grep -v WARN)
     #
     #
     ### include headers and libraries for RInside embedding classes
     RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave | grep -v WARN)
     RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave | grep -v WARN)
     #
     FLAGS         += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS | grep -v WARN)
     LINKFLAGS     += $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

endif  # IF MACOS


endif # IF MSL_STATIC

endif # IF MSL_R


ifeq ($(MSL_MSLOUT_DEBUG_OFF),T)
    FLAGS += -D__MSL_MSLOUT_DEBUG_OFF__
endif
# Generic Includes,Flags.  Static compile.  
# NOTE IS THE FOLLOWING STILL NECESSARY?
INCLUDE  = src -I${MSL_EXTERNAL_INCLUDE_DIR} 



# Add a MACOS flag for certain code breaks (see bottom of Tree.h, Selectable.h ... templated classes don't need pre-instantiations?)
ifeq ($(MSL_MACOS),T)
    FLAGS += -D__MACOS__ 
    MSL_STATIC=F
endif

ifeq ($(MSL_STATIC),T)
    FLAGS   += -static
endif

FLAGS   +=  -DUSE_REAL_EQ_DOUBLE 

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
objects: ${OBJECTS}
tests: ${TESTBINS}
examples: ${EXAMPLEBINS}
mybins: ${MYBINS}

## Create flags file
#objs/.flags::
#	echo "${FLAGS} -I${INCLUDE} ${SYMBOLS}" > objs/.flags
#
#${OBJECTS}: objs/%.o : src/%.cpp src/%.h 
#	${CC} @objs/.flags -c $< -o $@  
#
#${TESTBINS}: bin/% : tests/%.cpp objs/.flags ${OBJECTS} ${MYOBJS} ${HEADERS}
#	${CC} @objs/.flags ${LINKFLAGS} -Lobjs/ -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread
#
#${BINARIES}: bin/% : programs/%.cpp objs/.flags ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
#	${CC} @objs/.flags ${LINKFLAGS} -Lobjs/ -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread
#
#${EXAMPLEBINS}: bin/% : examples/%.cpp objs/.flags ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
#	${CC} @objs/.flags ${LINKFLAGS} -Lobjs/ -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread
#
#${MYOBJS}: objs/%.o : myProgs/%.cpp myProgs/%.h 
#	${CC} @objs/.flags -c $< -o $@  
#
#${MYBINS}: bin/% : myProgs/%.cpp objs/.flags ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
#	${CC} @objs/.flags ${LINKFLAGS} -o $@  ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}  -lpthread

${OBJECTS}: objs/%.o : src/%.cpp src/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@ 

${TESTBINS}: bin/% : tests/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${BINARIES}: bin/% : programs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${EXAMPLEBINS}: bin/% : examples/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${PHEADERS}
	${CC} ${FLAGS} ${LINKFLAGS} -Lobjs/ -I${INCLUDE} -o $@ ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS} -lpthread

${MYOBJS}: objs/%.o : myProgs/%.cpp myProgs/%.h 
	${CC} ${FLAGS} -I${INCLUDE} ${SYMBOLS} -c $< -o $@

${MYBINS}: bin/% : myProgs/%.cpp ${OBJECTS} ${MYOBJS} ${HEADERS} ${MYHEADERFILES}
	${CC} ${FLAGS} ${LINKFLAGS} -I${INCLUDE} -o $@  ${OBJECTS} ${MYOBJS} $< ${STATIC_LIBS}  -lpthread

.PHONY : clean
clean :
	-rm -f ${OBJECTS} ${BINARIES} ${EXAMPLEBINS} ${TESTBINS} ${MYOBJS} ${MYBINS}
#	-rm -f ${OBJECTS} ${BINARIES} ${EXAMPLEBINS} ${TESTBINS} ${MYOBJS} ${MYBINS} objs/.flags


#pythonLin: objs/.flags
#	gcc @objs/.flags -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/Library/Frameworks/Python.framework/Versions/2.6/Headers/
#	g++ @objs/.flags -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
#	cp PythonMSL.so /usr/share/python-support/pymol/pymol/
pythonLin:
	gcc -fpic -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/Library/Frameworks/Python.framework/Versions/2.6/Headers/
	g++ -lm -shared objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
	cp PythonMSL.so /usr/share/python-support/pymol/pymol/


#pythonMac: objs/.flags
#	gcc @objs/.flags -fPIC -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/System/Library/Frameworks/Python.framework/Versions/2.5/Headers/
#	g++ @objs/.flags -bundle -undefined dynamic_lookup objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
pythonMac:
	gcc -fPIC -c src/PythonMSL.cpp -o objs/PythonMSL.o -Wall -I${INCLUDE} -I/usr/include/python2.6 -I/usr/include  -I/System/Library/Frameworks/Python.framework/Versions/2.5/Headers/
	g++ -bundle -undefined dynamic_lookup objs/PythonMSL.o ${OBJECTS} ${STATIC_LIBS} -o PythonMSL.so 
#       sudo cp PythonMSL.so /Applications/PyMOLX11Hybrid.app/pymol/modules
