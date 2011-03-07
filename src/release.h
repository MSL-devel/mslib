/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
  Sabareesh Subramaniam, Ben Mueller

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#ifndef RELEASE_H
#define RELEASE_H

#define MSLVERSION "0.16.1.2"
#define MSLDATE "March 07, 2011"

/*
HISTORY:
0.16.1.2    March 07, 2011    jedonald
                'src/CharmmParameterReader.cpp' -Fix an uninitialized variable warning for BlockTypes block variable
0.16.1.1    March 07, 2011    jedonald
                'src/CharmmParameterReader.cpp' -Fix an uninitialized variable warning for BlockTypes block variable
0.16.1.0    March 05, 2011    asenes
                'src/Position.h' -changed setActiveIdentity to return bool
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added addIdentity functions in which the bb atoms are
                 passed as a string, space separated, of atoms
                'src/SystemRotamerLoader.h', 'src/SystemRotamerLoader.cpp' -CHANGED API: in loadRotamers and addRotamers the rotlib
                 is passed as a last argument and defaulted to blank, which means use default library. Old functions are still
                 there but declared DEPRECATED
                'src/System.h' -Added setActiveIdentity functions
                'src/Quench.cpp', 'programs/grepSequence.cpp', 'programs/fillInSideChains.cpp', 'programs/coiledCoilBuilder.cpp',
                 'programs/energyOptimizations.h', 'tests/testCharmmEnergies.cpp' -Fixed loadRotamers functions to reflect API
                 change in SystemRotamerLoader
                'Makefile', 'examples/examples.mk' -Commented out programs that do not compile
                'tests/testTokenize.cpp' -Added test
0.16.0.0    March 03, 2011    asenes
                'src/AtomContainer.h', 'src/HelixGenerator.cpp', 'src/Position.h', 'src/FuseChains.cpp', 'src/Position.cpp', 'src/Residue.h',
                 'src/ChiStatistics.h', 'src/Chain.cpp', 'src/Chain.h', 'src/IcTable.h', 'src/BBQTable.cpp', 'src/PDBTopology.cpp',
                 'src/IcTable.cpp', 'src/EnvironmentDescriptor.cpp', 'src/PDBTopology.h', 'src/CharmmSystemBuilder.cpp', 'src/System.h',
                 'src/EnvironmentDatabase.cpp', 'src/SelfPairManager.cpp', 'src/ChiStatistics.cpp', 'src/System.cpp', 'src/UserDefinedEnergySetBuilder.cpp',
                 'src/HelixFusion.cpp', 'src/BackRub.cpp', 'src/PolymerSequence.cpp', 'programs/insertLoopIntoTemplate.cpp', 'programs/printSequence.cpp',
                 'programs/calculateSasa.cpp', 'programs/getSphericalCoordinates.cpp', 'programs/fillInSideChains.cpp', 'programs/getSurroundingResidues.cpp',
                 'tests/testResiduePairTable.cpp', 'tests/testAddCharmmIdentity.cpp', 'tests/testCharmmEnergies.cpp', 'tests/testHelixFusion.cpp',
                 'tests/testPDBFragments.cpp', 'tests/testLinkedPositions.cpp', 'tests/testPhiPsi.cpp', 'tests/testResidueSubstitutionTable.cpp',
                 'examples/example_looping_over_Chain_Residues_Atoms.cpp' -API CHANGE! Removed the size() function from the molecular
                 objects, replaced with chainSize(), positionSize(), residueSize(). Also removed deprecated exists() function
0.15.2.0    March 02, 2011    asenes
                'src/PDBFormat.h', 'src/PDBFormat.cpp' -Added support for the MODEL and ENDMDL tags
                'src/PDBReader.h', 'src/PDBReader.cpp' -Now it counts the number of models (using the MODEL tag)
                'src/System.h', 'src/System.cpp' -Now it gets from the PDBReader the number of models and can switch to a given
                 model (using internally the alternative coordinates of the Atoms)
                'examples/example_multiple_coordinates_from_NMR_multiModel_PDB.cpp', 'exampleFiles/example0007.pdb', 'examples/examples.mk'
                 -Added example for reading and accesing multi model PDB files with the System
                'examples/example_SasaCalculator_usage.cpp' -Fixed text inaccuracy
0.15.1.0    February 20, 2011    asenes
                'src/System.h' -added setActiveRotamer(std::string _identityOrPositionId, unsigned int _n) to set the position
                 in its n
                'src/Position.h' -added void setActiveRotamer(std::string _identity, unsigned int _n) to set the position to the
                 n
0.15.0.1    February 14, 2011    brettth
                'src/SurfaceAreaAndVolume.cpp' -Fixed bug that was causing volume to be incorrectly calculated.
                'src/SurfaceAreaAndVolume.cpp' -Fixed bug that was causing volume to be incorrectly calculated.
0.15.0.0    February 12, 2011    asenes
                'examples/example_add_atoms_to_System_and_AtomContainer.cpp' -Example on how to add atoms to the AtomContainer
                 and the System with the addAtom function
                'examples/examples.mk' -Added example_add_atoms_to_System_and_AtomContainer.cpp
                'src/AtomContainer.h', 'src/AtomContainer.cpp', 'src/Residue.h', 'src/Residue.cpp' -added element to addAtom function
                
                'src/System.h', 'src/System.cpp' -added addAtom function
                'src/Position.cpp' -fixed bug, it was not updating the atom list correctly if atoms were added multiple times
                'src/CharmmTopologyReader.cpp' -fixed small bug
                'src/CartesianPoint.h', 'src/CartesianPoint.cpp' -added angleRadians and dihedralRadians functions
                'src/BBQTable.cpp' -initialized pointers to NULL to remove warning
                'src/File.cpp' -initialized openmode variable to remove warning
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -removed deprecated functions. added fail() function
                 to check for errors after reading input files
                'src/Transforms.cpp' -removed verbose MSLOUT stream
                'Makefile' -removed external flag system because it was not always working well
0.14.1.3    October 21, 2010    brettth
                'src/RandomNumberGenerator.cpp' -Getting rid of copying vector of doubles into double * buffer.
0.14.1.2    September 08, 2010    sabs
                'src/RotamerLibraryBuilder.h' -Added comments on usage.
                'src/RotamerLibraryBuilder.cpp' -Removed unnecessary check for null string in libName
0.14.1.1    September 01, 2010    bkmueller
                'programs/coiledCoilBuilder.cpp' -Created hack so parameter loop goes from start to finish correctly
0.14.1.0    September 01, 2010    bkmueller
                'tests/testRandomNumberGenerator.cpp', 'tests/testDerivatives.cpp' -Forgotten commits for version 0.14.0.0 when
                 the API of RandomNumberGenerator was changed
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Now it takes the variable positions from the System. The System
                 can be told what are the variable positions or do it automatically
                'src/System.h', 'src/System.cpp' -The system can now be told what are the variable positions. Added setVariablePositions
                 and isPositionVariable. Added writePdb with remarks.
                'src/MslTools.h', 'src/MslTools.cpp' -Added comparison functions for atomId, identityId, positionId, etc. compareAtomIds,
                 comparePositionIds, compareIdentityIds, compareAtomOfIdentityIds
                'src/OptionParser.h', 'src/OptionParser.cpp' -Added getUnsignedIntVector, getUnsignedIntVector and getUnsignedMultiIntVector
                 functions
                'programs/coiledCoilBuilder.cpp' -Completed first version of program, supports generation of coils of single parameter
                 with side chain repacking, or grid search over parameter space
0.14.0.3    August 27, 2010    sabs
                'src/SelfPairManager.cpp' -Another bugfix
0.14.0.2    August 24, 2010    sabs
                'src/SelfPairManager.cpp' -Bug fix
0.14.0.1    August 16, 2010    dwkulp
                'Makefile' -define flags once for compilation/linking put into objs/.flags file, results in cleaner compilation/linking
                 output
0.14.0.0    August 12, 2010    bkmueller
                'programs/coiledCoilBuilder.cpp' -Takes in a polymer sequence and CC parameters and uses CoiledCoils and SelfPairManager
                 to create the CC and optimize its residues
                'src/MonteCarloManager.h', 'src/MonteCarloManager.cpp' -Ported from Cub, runs any Monte Carlo Optimization and
                 accepts and rejects postions based on the Metropolis Criterion (takes energy as input, keeps track of temperature
                 and cycles)
                'src/SelfConsistentMeanField.cpp', 'src/SelfConsistentMeanField.h' -Ported from Cub, runs a Self Consistent Mean
                 Field algorithm on side chains and finds the probabilities of each rotamer at each position, can accept the masks
                 created by Dead End Elmination and does not consider these rotamers
                'tests/testCoiledCoils.cpp' -A limited testing program for the CoiledCoil object
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added setRandomNumberGenerator, and getRandomNumberGenerator.
                 Also seed and getSeed functions. Also functions getFixedEnergy, getSelfEnergy and getPairEnergy. Added the runOptimizer
                 function which optimizes the rotamers of a side chain based on Dead End Elimination, Enumeration, Self Consistent
                 Mean Field and Monte Carlo Optimization, all can be turned off or on by the functions: setRunDEE, setRunMC, setRunSCMF
                 and setRunEnum. A setVerbose function can be toggled on and off to display comments from the runOptimizer function.
                 The results of MCO, SCMF and DEE can be returned by getDEEAliveRotamers, getDEEAliveMask, getSCMFstate, and getMCOstate
                
                'src/Enumerator.h', 'src/Enumerator.cpp' -Added new constructor to take a vector<vector<unsigned int> > and calcEnuerationValues
                 and setValues functions
                'src/CoiledCoils.h', 'src/CoiledCoils.cpp' -removed offersCoiledCoils and sotoCoiledCoils. Merged northCoiledCoils
                 and gevorgCoiledCoils into northCoiledCoils which can now take both Gevorg (renamed Crick) parameters or Norths
                 parameters. Added getCoiledCoilBundle and getCoiledCoilBundleCricks which create CN and DN symmetry bundles of
                 Coiled Coils. Also added primarySequenceToCoiledCoils which takes a system and applies CoiledCoil geometry to
                 it based on input starting parameters. All functions parameters are taken in degrees
                'src/MonteCarloOptimization.cpp' -Adjusted functions to work with new RandomNumberGenerator API
                'src/DeadEndElimination.h', 'src/DeadEndElimination.cpp' -Changed getAliveStates to return a vector<vector<unsigned
                 int> > instead of a vector<vector<int> >
                'src/RandomNumberGenerator.h', 'src/RandomNumberGenerator.cpp' -Changed API from setRNGSeed to setSeed, getRNGSeed
                 to getSeed and setRNGTimeBasedSeed to setTimeBasedSeed
                'src/MslTools.h', 'src/MslTools.cpp' -Added normalizeVector and normalizeCumulativeVector functions
                'src/Symmetry.h', 'src/Symmetry.cpp' -Modified applyCN and applyDN so now the newly created symmetry atoms can
                 be pushed into the original AtomPointerVector (not default behavior)
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp' -added subdivideByChainAndPosition and subdivideByChainPositionAndIdentity
                 which divides the AtomPointerVector into a vector<vector<map<string, Atom*> > > or a vector<vector<map<string<map,
                 Atom*> > > >
                'tests/testCCD.cpp' -Adjusted to work with new RandomNumberGenerator API
                'Makefile' -Added new programs and objects
0.13.0.0    August 10, 2010    asenes
                'src/RandomNumberGenerator.h', 'src/RandomNumberGenerator.cpp' -Removed dependency on GSL (usage with GSL still
                 preferred) and revised the API
                'tests/testRandomNumberGenerator.cpp' -Revised test for RandomNumberGenerator
                'src/SurfaceAreaAndVolume.cpp', 'src/CCD.cpp', 'src/PhiPsiStatistics.h', 'src/PhiPsiStatistics.cpp', 'src/MonteCarloOptimization.cpp',
                 'src/Quench.cpp', 'src/CrystalLattice.h', 'src/CrystalLattice.cpp', 'src/BackRub.cpp' -Removed GLS dependency
                 for RandomNumberGenerator
                'src/MslTools.h', 'src/MslTools.cpp' -Removed GLS dependency for RandomNumberGenerator and removed unnecessary
                 function getRandomInt (use RandomNumberGenerator instead)
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp', 'src/AtomicPairwiseEnergy.h', 'src/AtomicPairwiseEnergy.cpp'
                 -Removed call for MslTools getRandomInt, shifted to RandomNumberGenerator
                'Makefile' -Removed GLS dependency for RandomNumberGenerator, PhiPsiStatistics, BackRub, CCD, MonteCarloOptimization,
                 Quench, SurfaceAreaAndVolume
                'src/Atom.h' -Fixed bug in removeAllAltConformations
0.12.0.0    August 02, 2010    sabs
                'src/Scrwl4HBondInteraction.cpp', 'src/Scrwl4HBondInteraction.h', 'src/HydrogenBondBuilder.h', 'src/HydrogenBondBuilder.cpp'
                 -Implementation of the SCRWL Hydrogen Bond term
                'src/EnergySet.cpp', 'src/EnergySet.h' -Added _activeOnly flag to the calcEnergyOfSubset Function
                'src/Reader.cpp', 'src/Reader.h' -Added getAllLines to the Reader
                'Makefile' -Added the hydrogenBondBuilder and Scrwl classes
0.11.0.0    July 14, 2010    dwkulp
                'src/SystemRotamerLoader.cpp', 'src/SystemRotamerLoader.h' -stores name of rotamer library file
                'src/Position.h', 'src/Position.cpp' -New constructor taking a positionId string
                'src/Atom.h', 'src/Atom.cpp' -functions to set HasCoordinates flag and a minimization Index variable
                'src/ChiStatistics.h', 'src/ChiStatistics.cpp' -reorganize chi definitions into PDBTopology
                'src/AtomPointerVector.cpp' -Turn off MslOut in constructor
                'src/CharmmTopologyReader.h', 'src/CharmmTopologyReader.cpp' -store charmmFileName
                'src/PDBTopology.h', 'src/PDBTopology.cpp', 'tests/testPDBTopology.cpp', 'Makefile', 'examples/example_mutation_rotamers.cpp',
                 'examples/examples.mk' -Added PDBToplogy object
0.10.0.0    May 19, 2010    dwkulp
                'src/RandomNumberGenerator.h', 'src/RandomNumberGenerator.cpp', 'src/PhiPsiStatistics.h', 'src/PhiPsiStatistics.cpp',
                 'tests/testRandomNumberGenerator.cpp', 'tests/testPhiPsi.cpp', 'Makefile' -Added Non
0.9.0.1    May 17, 2010    dwkulp
                'examples/example_coiled_coils_and_symmetric_bundles.cpp', 'examples/examples.mk' -coiled coil and symmetric bundle
                 example
0.9.0.0    May 17, 2010    dwkulp
                'Makefile', 'src/OptionParser.h', 'src/OptionParser.cpp', 'src/MslOut.h', 'src/MslOut.cpp', 'tests/testMslOut.cpp',
                 'tests/testMslOut2.cpp', 'src/Atom.cpp', 'src/Transforms.cpp', 'src/AtomPointerVector.cpp' -Added MSLOUT functionality
                 to some objects for testing, plus some tests
0.8.9.0    April 30, 2010    asenes
                'src/CRDWriter.h', 'src/CRDWriter.cpp', 'src/CRDFormat.h', 'src/CRDFormat.cpp', 'src/CRDFormat.h', 'src/CRDReader.h',
                 'src/CRDReader.cpp', 'tests/testCRDIO.cpp', 'tests/testData.h' -CRD IO, added CRD writer, removed a few bugs from
                 CRDReader and CRDFormat
                'src/PDBWriter.cpp' -Removed annoying REMARK credit line
                'programs/alignMolecules.cpp' -Added option not to output a transformed pdb file (just calculate the alignment
                 and RMSD
                'programs/setConformation.cpp' -Fixed help
                'Makefile' -Added CRDWriter
0.8.8.2    April 29, 2010    asenes
                'programs/alignMolecules.cpp' -Fixed bug, not longer using PDBReader and PDBWriter
0.8.8.1    April 29, 2010    sabs
                'tests/testPolymerSequence.cpp' -Added the API addPositionIdentity
                'tests/testData.h' -Added data for CRDReader
                'src/SystemRotamerLoader.cpp' -Added check to keep old Rotamers when the keepOldRotamers flag is on
                'src/RotamerLibrary.cpp' -Increased the number of digits in Rotamer Library
                'src/CharmmTopologyReader.h' -Added API getElement(atomName)
                'Makefile' -Added CRDReader and its testFiles
                'src/CRDReader.cpp', 'src/CRDReader.h', 'src/CRDFormat.cpp', 'src/CRDFormat.h', 'tests/testCRDIO.cpp' -CRDReader,CRDFormat
                 objects and test files
                'src/Atom.h' -getCharge returns double instead of double &
0.8.7.1    April 29, 2010    dwkulp
                'src/CharmmTopologyReader.h' -helper function getResidues() will return a vector of CharmmTopologyResidues
0.8.7.0    April 28, 2010    dwkulp
                'src/CharmmElectrostaticInteraction.h', 'src/CharmmEnergy.h', 'src/CharmmEnergy.cpp', 'src/AtomicPairwiseEnergy.cpp',
                 'src/Quench.cpp', 'src/AtomicPairwiseEnergy.h', 'src/CharmmVdwInteraction.h' -Switchable cutoffs modified, push
                 utility functions into CharmmEnergy such that AtomicPairwiseEnergy could also use them
                'src/SurfaceAreaAndVolume.cpp', 'src/SurfaceAreaAndVolume.h', 'tests/testSurfaceAreaAndVolume.cpp' -Bug fixed in
                 SurfaceArea, now testSurfaceAreaAndVolume computes proper surface area.... still working on volume
                'src/Atom.cpp' -setUnbound.. function commented out in destructor, memory leak
                'programs/runQuench.cpp', 'programs/runQuench.h' -autoFindPositions option added, non bonded cutoffs added
                'programs/insertLoopIntoTemplate.h', 'programs/insertLoopIntoTemplate.cpp' -numClashes arguement, rename chains
                 such that there isn't a duplicate chain in output file
                'programs/getSelection.h', 'programs/getSelection.cpp' -select only residues found in a CHARMM toplogy file, rename
                 HIS
                'programs/analEnergy.h', 'programs/analEnergy.cpp' -compute energy between two selections
                'programs/getDihedrals.cpp' -added filename to R plot output
0.8.6.1    April 27, 2010    asenes
                'src/SystemRotamerLoader.cpp' -Fixed bug, it was not updating the variable position table in the System
0.8.6.0    April 16, 2010    asenes
                'src/Transforms.h', 'src/Transforms.cpp' -Added setImproper function, fixed bug in setDihedral function
                'programs/setConformation.cpp' -A program to edit degrees of freedom (bond distances, angles, dihedrals, impropers)
                 from a PDB file
		 'Makefile' -Updated with setConformation program
0.8.5.2    April 09, 2010    dwkulp
                'examples/examples.mk', 'examples/example_measurements.cpp' -Measurement example
                'src/MslTools.h' -comment line added for getFileName
0.8.5.1    April 09, 2010    jedonald
                'tests/testTransforms.cpp', 'tests/testGenerateCrystalLattice.cpp', 'tests/testSymmetry.cpp', 'tests/testDerivatives.cpp',
                 'tests/testBBQ2.cpp', 'src/EnvironmentDescriptor.cpp', 'src/RotamerLibraryBuilder.cpp', 'src/HelixGenerator.cpp',
                 'src/LinearProgrammingOptimization.cpp', 'src/CartesianGeometry.cpp', 'src/SurfaceAreaAndVolume.cpp', 'src/PythonMSL.cpp',
                 'src/Matrix.cpp', 'src/BBQTableReader.cpp', 'src/EnvironmentDatabase.cpp', 'src/Atom.h', 'src/PhiPsiStatistics.cpp',
                 'src/ChiStatistics.cpp', 'src/AtomPointerVector.cpp', 'src/Line.cpp', 'src/Frame.cpp', 'src/CartesianGeometry.h',
                 'src/System.cpp', 'src/Transforms.cpp', 'src/CartesianPoint.h', 'src/Helanal.cpp', 'src/Matrix.h', 'src/IcEntry.cpp',
                 'src/MslTools.h', 'src/CartesianPoint.cpp', 'src/CrystalLattice.cpp', 'src/HelixFusion.cpp', 'src/BackRub.cpp',
                 'src/BBQTable.cpp', 'src/CoiledCoils.cpp', 'src/Symmetry.cpp', 'programs/generateCoiledCoils.cpp', 'programs/optimizeMC.cpp',
                 'programs/printSequence.cpp', 'programs/getSphericalCoordinates.cpp', 'programs/optimizeLP.cpp' -Make MslTools
                 a sub
0.8.5.0    April 07, 2010    dwkulp
                'src/Chain.h', 'src/Chain.cpp', 'src/Position.h', 'src/Position.cpp', 'src/Residue.cpp', 'Makefile', 'src/FuseChains.h',
                 'src/FuseChains.cpp', 'programs/insertLoopIntoTemplate.cpp', 'programs/insertLoopIntoTemplate.h', 'src/AtomContainer.cpp'
                 -Updated position/chain indexing functions, modified findNeighbors in Residue, removeAllAtoms in AtomContainer,
                 Makefile environment variables, new object FuseChains and program insertLoopIntoTemplate
0.8.4.0    April 07, 2010    asenes
                'src/AtomSelection.cpp' -Fixed bug, it was not clearning the AtomPointerVector when repeating a selection with
                 the same name
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp', 'tests/testCharmmBuild.cpp' -Added buildSystemFromPDB
                 function
0.8.3.2    April 06, 2010    dwkulp
                'programs/getDihedrals.cpp' -Fancy Rama plots get printed in getDihedrals through R
0.8.3.1    April 06, 2010    dwkulp
                'Makefile', 'programs/getDihedrals.cpp', 'tests/testRInterface.cpp' -R interface for MSL
0.8.3.0    April 05, 2010    asenes
                'src/Atom.h', 'src/Atom.cpp' -Now bonds can be delete from atoms. The Atom calls its bonded atoms to delete their
                 bond reciprocally. Called by the distruptor
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added addIdentity function to add new identity post
                 built. Also change the API, the System is passed upon construction
                'src/Symmetry.h', 'src/Symmetry.cpp' -Minor fixes, added copy contructor and = operator
                'src/System.h', 'src/System.cpp' -Now reset is a public function.
                'src/AtomSelection.h', 'src/AtomSelection.cpp' -Made the new selection the default (the old is still available
                 compiling under F = T.
                'src/SystemRotamerLoader.h', 'src/SystemRotamerLoader.cpp' -Now using the positionId in loadRotamers. Removed references
                 to the deprecated exists function
                'src/Residue.h', 'src/BBQTable.h', 'src/BBQTable.cpp' -Removed references to the deprecated exists function
                'src/LogicalCondition.h', 'src/LogicalCondition.cpp' -Second version, fixes and expanded
                'src/Quench.h', 'src/Quench.cpp', 'tests/testCharmmBuild.cpp', 'tests/testCharmmEnergies.cpp', 'tests/testEEF1.cpp',
                 'tests/testEEF1_2.cpp', 'tests/testEnergeticAnalysis.cpp', 'tests/testLinkedPositions.cpp', 'tests/testNonBondedCutoff.cpp',
                 'tests/testSurfaceAreaAndVolume.cpp', 'programs/analEnergy.cpp', 'programs/energyOptimizations.h', 'programs/fillInSideChains.cpp'
                 -Adjusted for change of API of CharmmSystemBuilder
                'src/CharmmTopologyResidue.h', 'src/CharmmTopologyResidue.cpp' -Only updated file header
                'tests/testAddCharmmIdentity.cpp' -A test for the addIdentity in the CharmmSystemBuilder
                'tests/testAtomSelection.cpp' -Revised with more complex logic examples
                'tests/testCoiledCoils.cpp' -Now it writes the output file to tmp
                'tests/testAtomBondBuilder.cpp' -Added tests for the removal of bonds from atoms
                'examples/example_selecting_atoms.cpp' -Replaces the previous example_selecting_atoms_and_residues.cpp which wasn't
                 showing residues (need another program for that)
                'examples/examples.mk' -
                'Makefile' -added testAddCharmmIdentity
0.8.2.6    April 01, 2010    dwkulp
                'examples/examples.mk' -UPdated examples makefile
0.8.2.5    April 01, 2010    dwkulp
                'exampleFiles/example0004.pdb', 'examples/example_regular_expressions.cpp' -Molecular Alignment Example
0.8.2.4    April 01, 2010    dwkulp
                'exampleFiles/example0005.pdb', 'exampleFiles/example0006.pdb', 'examples/example_molecular_alignment.cpp' -Molecular
                 Alignment Example
0.8.2.3    March 31, 2010    jedonald
                'src/Symmetry.h', 'src/Symmetry.cpp', 'programs/generateCoiledCoils.cpp' -Fix symmetry object to correctly do D_N.
                 Also fix hardcoded values to match Gevorgs data
0.8.2.2    March 30, 2010    dwkulp
                'src/CoiledCoils.cpp' -sotoCoiledCoil method added
0.8.2.1    March 26, 2010    asenes
                'src/LogicalCondition.cpp', 'src/LogicalCondition.h' -Fixed a number of bugs
                'src/AtomSelection.cpp' -Added support for keyword HASCOOR in addition to HASCRD
                'tests/testAtomSelection.cpp' -Added more complex tests, fixed bug
0.8.2.0    March 25, 2010    asenes
                'src/MslTools.h', 'src/MslTools.cpp' -Added toBool function
                'src/System.h', 'src/System.cpp', 'src/AtomContainer.h', 'src/AtomContainer.cpp' -Fixed small bug, applySavedCoor
                 was not returning a bool
                'src/LogicalCondition.h', 'src/LogicalCondition.cpp', 'src/AtomSelection.cpp', 'src/AtomSelection.h', 'tests/testAtomSelection.cpp'
                 -A new logical parser and tree for selection logic (LogicalCondition): currently only being tested. one needs
                 to setenv MSL_TESTING T to turn on compiling of this feature
                'src/Atom.h' -Minimal changes, do not affect the code
                'Makefile' -Added LogicalCondition and support for another conditional environmental varialbe T, which turn on
                 compilation of testing code that is off by default
0.8.1.1    March 24, 2010    dwkulp
                'exampleFiles/example0003.pdb', 'examples/examples.mk', 'examples/example_selecting_atoms_and_residues.cpp' -Example
                 Selecting Atoms And Residues
0.8.1.0    March 19, 2010    asenes
                'src/HelixGenerator.cpp', 'programs/grepSequence.cpp', 'programs/alignMolecules.cpp', 'programs/getSurroundingResidues.cpp'
                 -Change of API: change the Transforms function align to rmsdAlignment
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp' -Removed function updateGeometricCenter. Now the cartesian
                 center is calculated directly by the getGeometricCenter function, unless and integer stamp identical to the previous
                 one is given, in which case a chached center is returned (this is good if we need to get the center multiple times
                 and the atoms do not move
                'src/Quaternion.cpp', 'src/HelixFusion.h', 'src/HelixFusion.cpp', 'src/Frame.cpp', 'src/PDBFragments.h', 'src/PDBFragments.cpp',
                 'src/PrincipleComponentAnalysis.cpp', 'src/CrystalLattice.cpp' -Change of API: removed the AtomPointerVector's
                 updateGeometricCenter function
                'src/AtomGroup.h', 'src/AtomGroup.cpp' -Now it uses the parent's (AtomPointerVector) getGeometricCenter function
                
                'src/Transforms.h', 'src/Transforms.cpp' -Changed the RMSD alignment function align (which had the same name of
                 the function that aligns one vector with another) to the more clear rmsdAlignment name
                'examples/example_AtomPointerVector.cpp' -Refined the example program
0.8.0.0    March 19, 2010    asenes
                'examples/example_looping_over_Chain_Residues_Atoms.cpp', 'exampleFiles/example0002.pdb' -An example program for
                 looping over chains residues and atoms in the System
                'examples/example_read_write_PDBs_with_the_AtomContainer.cpp', 'examples/example_read_write_PDBs_with_the_System.cpp'
                 -Minor changes
                'examples/example_AtomPointerVector.cpp' -An example program for using the AtomPointerVector
                'programs/generateCoiledCoils.cpp', 'tests/testGenerateCrystalLattice.cpp', 'tests/testNonBondedCutoff.cpp', 'tests/testSymmetry.cpp',
                 'tests/testCoiledCoils.cpp', 'tests/testCharmmEnergies.cpp' -API changes: removed translate and rotate from AtomPointerVector;
                 also getAtoms changed to getAtomPointers in many objects and getAllAtoms in the System changed to getAllAtomPointers
                
                'programs/createFragmentDatabase.cpp', 'programs/runKBQuench.cpp', 'programs/getSelection.cpp', 'programs/analEnergy.cpp',
                 'programs/grepSequence.cpp', 'programs/alignMolecules.cpp', 'programs/calculateSasa.cpp', 'programs/getSphericalCoordinates.cpp',
                 'programs/fillInSideChains.cpp', 'programs/runQuench.cpp', 'programs/getSurroundingResidues.cpp', 'programs/energyOptimizations.h',
                 'tests/testTransforms.cpp', 'tests/testTransformBondAngleDiheEdits.cpp', 'tests/testPDBFragments.cpp', 'tests/testLinkedPositions.cpp',
                 'tests/testSystemIcBuilding.cpp', 'tests/testPolymerSequence.cpp', 'tests/testQuench.cpp', 'tests/testResiduePairTable.cpp',
                 'tests/testCharmmBuild.cpp', 'tests/testFrame.cpp', 'tests/testPhiPsi.cpp', 'tests/testPDBIO.cpp', 'tests/testSasaCalculator.cpp',
                 'tests/testEEF1_2.cpp', 'tests/testEnvironmentDescriptor.cpp', 'tests/testBBQ.cpp', 'tests/testEnergySet.cpp',
                 'tests/testResidueSubstitutionTable.cpp', 'tests/testLoopOverResidues.cpp', 'tests/testEnergeticAnalysis.cpp',
                 'tests/testAtomBondBuilder.cpp', 'tests/testSystemCopy.cpp', 'tests/testEEF1.cpp', 'tests/testBBQ2.cpp', 'tests/testEnvironmentDatabase.cpp',
                 'tests/testCCD.cpp', 'tests/testSurfaceAreaAndVolume.cpp', 'tests/testHelixFusion.cpp' -API changed, getAtoms
                 changed to getAtomPointers in many objects and getAllAtoms in the System changed to getAllAtomPointers
                'src/Chain.h', 'src/Chain.cpp' -Added toString and << operator; also API changed, getAtoms changed to getAtomPointers
                 in many objects and getAllAtoms in the System changed to getAllAtomPointers
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp' -API changed, : removed translate and rotate from AtomPointerVector,
                 use the Transforms object instead
                'src/Residue.h' -toString now spits the identityId; also API changed, getAtoms changed to getAtomPointers in many
                 objects and getAllAtoms in the System changed to getAllAtomPointers
                'src/PolymerSequence.h' -Added addPositionIdentity function; also API changed, getAtoms changed to getAtomPointers
                 in many objects and getAllAtoms in the System changed to getAllAtomPointers
                'src/Position.h', 'src/Position.cpp' -toString now prints the positionId followed by the list of all identities;
                 also API changed, getAtoms changed to getAtomPointers in many objects and getAllAtoms in the System changed to
                 getAllAtomPointers
                'src/SurfaceAreaAndVolume.h', 'src/SurfaceAreaAndVolume.cpp', 'src/CrystalLattice.h', 'src/CrystalLattice.cpp'
                 -API changes: removed translate and rotate from AtomPointerVector
                'src/Symmetry.h', 'src/Symmetry.cpp' -API changes: removed translate and rotate from AtomPointerVector; also removed
                 translate and rotate from AtomPointerVector
                'src/EnergeticAnalysis.cpp', 'src/EnvironmentDescriptor.cpp', 'src/PDBFragments.cpp', 'src/TwoBodyDistanceDependentPotentialTable.cpp',
                 'src/HelixFusion.cpp', 'src/SasaCalculator.cpp', 'src/PSFReader.h', 'src/SystemRotamerLoader.cpp', 'src/AtomContainer.h',
                 'src/CharmmSystemBuilder.cpp', 'src/BBQTable.cpp', 'src/Quench.cpp', 'src/EnvironmentDatabase.cpp', 'src/PythonMSL.cpp',
                 'src/BackRub.cpp', 'src/SasaCalculator.h', 'src/AtomGeometricRelationship.h', 'src/CoiledCoils.h', 'src/Transforms.cpp',
                 'src/EnergySet.cpp', 'src/RegEx.cpp', 'src/SelfPairManager.cpp', 'src/PDBReader.h', 'src/System.h', 'src/System.cpp',
                 'src/AtomicPairwiseEnergy.cpp', 'src/CCD.cpp', 'src/Interaction.h' -API changed, getAtoms changed to getAtomPointers
                 in many objects and getAllAtoms in the System changed to getAllAtomPointers
                'examples/example_SasaCalculator_usage.cpp' -API changed, getAtoms changed to getAtomPointers in many objects and
                 getAllAtoms in the System changed to getAllAtomPointers
                'examples/examples.mk' -Added example_looping_over_Chain_Residues_Atoms example_AtomPointerVector
0.7.1.3    March 16, 2010    asenes
                'src/System.h' -Just edited a comment
                'src/AtomContainer.h', 'src/AtomContainer.cpp' -Added the ability to save coordinates to buffer; added operator()(unsigned
                 int) and operator[](string); replaced all size_t with unsigned int;
                'examples/example_read_write_PDBs_with_the_AtomContainer.cpp', 'examples/example_read_write_PDBs_with_the_System.cpp'
                 -Example programs for the tutorial: read and write PDB files
                'examples/examples.mk' -Added example_read_write_PDBs_with_the_AtomContainer and example_read_write_PDBs_with_the_System
                
                'Makefile' -added testResidueSelection (forgotten in previous submit)
0.7.1.2    March 16, 2010    sabs
                'src/System.h' -Support to save atom coordinates
                'src/CharmmSystemBuilder.h' -Added Support to read Solvation Parameters
                'src/CharmmSystemBuilder.cpp' -Added Support to read Solvation Parameters
0.7.1.1    March 15, 2010    jedonald
                'programs/getSurroundingResidues.cpp' -Fix output filenames for non
0.7.1.0    March 13, 2010    dwkulp
                'src/Selectable.h', 'src/LogicalParser.cpp', 'src/Residue.cpp', 'tests/testResidueSelection.cpp' -Added new functionality
                 to Selectable objects so you can now selecton on functions that return a bool but also pass in a single string
                 arguement
0.7.0.2    March 09, 2010    sabs
                'Makefile' -Moved programs and tests to appropriate sections of the Makefile based on libraries required
                'programs/analEnergy.cpp' -Removed redundant using namespace MSL
0.7.0.1    March 09, 2010    asenes
                'src/Atom.cpp' -Forgot to upload at rel 0.7.0.0. Added getter for atomId and atomOfIdentityId
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Set R
0.7.0.0    March 09, 2010    asenes
                'tests/testEEF1.cpp', 'tests/testEEF1_2.cpp' -A test for the implementation of Lazaridis's solvation
                'tests/testCharmmEEF1ParameterReader.cpp' -A test for the implementation of the parameter reader for Lazaridis's
                 solvation
                'tests/testAtomContainer.cpp' -A test for the AtomContainer
                'tests/testAtomAndResidueId.cpp' -A test for the newly implemented atom, position and identity ids in MslTools
                
                'src/CharmmEEF1Interaction.h', 'src/CharmmEEF1Interaction.cpp' -Interaction object for the Lazaridis EEF1 solvation,
                 two
                'src/CharmmEEF1ParameterReader.h', 'src/CharmmEEF1ParameterReader.cpp' -The parameter file reader for the Lazaridis
                 EEF1 solvation
                'src/OneBodyInteraction.h', 'src/OneBodyInteraction.cpp' -A one
                'src/CharmmEEF1RefInteraction.h', 'src/CharmmEEF1RefInteraction.cpp' -Interaction object for the Lazaridis EEF1
                 solvation, one
                'tests/testResidueSubstitutionTable.cpp' -Fixed after atom/res/pos/chain
                'tests/testEnergySet.cpp' -Minor. Added now needed includes
                'tests/testPhiPsi.cpp', 'tests/testLoopOverResidues.cpp', 'tests/testRegEx.cpp', 'tests/testBBQ2.cpp', 'tests/testCharmmEnergies.cpp',
                 'tests/testTransformBondAngleDiheEdits.cpp', 'tests/testNonBondedCutoff.cpp', 'tests/testLinkedPositions.cpp',
                 'tests/testResiduePairTable.cpp', 'src/EnergeticAnalysis.cpp', 'src/PDBFragments.cpp', 'src/HelixFusion.cpp',
                 'src/BBQTable.cpp', 'src/EnvironmentDatabase.cpp', 'src/PythonMSL.cpp', 'src/BackRub.cpp', 'src/ResidueSelection.cpp',
                 'src/HelixGenerator.cpp', 'programs/runKBQuench.cpp', 'programs/analEnergy.cpp', 'programs/grepSequence.cpp',
                 'programs/getDihedrals.cpp', 'programs/getSphericalCoordinates.cpp', 'programs/fillInSideChains.cpp', 'programs/energyOptimizations.h',
                 'programs/runQuench.cpp', 'examples/example_AtomContainer_usage.cpp', 'examples/example_multipleResidueIdentities.cpp',
                 'programs/printSequence.cpp', 'programs/getSurroundingResidues.cpp' -Fixedafteratom/res/pos/chain
                'tests/testCharmmBuild.cpp' -Added calculation of the energies
                'src/PolymerSequence.h', 'src/PolymerSequence.cpp' -Added a function to get a sequence directly from reading a
                 PDB file, simplified the implementation of a constructor. This object needs some rethinking.
                'src/Residue.h', 'src/Residue.cpp', 'src/Position.h', 'src/Position.cpp', 'src/Chain.h', 'src/Chain.cpp', 'src/System.h',
                 'src/System.cpp' -Changed the atom/res/pos/chain
                'src/CharmmElectrostaticInteraction.cpp' -Fixed bug, added missing initialization of non
                'src/CharmmParameterReader.h', 'src/CharmmParameterReader.cpp' -Changed all the functions to get param to return
                 a bool, removed functions for EEF1 solvation (they were not implemented), UB now added and get
                'src/Atom.h' -Added getter for atomId and atomOfIdentityId
                'src/Selectable.h' -Added clearFlag and clearAllFlags functions
                'src/Matrix.h', 'src/Matrix.cpp' -???? Added gsl includes
                'src/ALNReader.cpp', 'src/PSFReader.cpp', 'src/TBDReader.cpp', 'src/RotamerLibraryReader.cpp', 'src/PDBReader.cpp',
                 'src/PhiPsiReader.cpp', 'src/CharmmTopologyReader.cpp', 'src/ResiduePairTableReader.cpp', 'src/BBQTableReader.cpp',
                 'src/MIDReader.cpp' -The read() function of the reader returns false if the file wasn't open
                'src/Reader.cpp' -Fixed memory leak in read(string &_inputString)
                'tests/testPDBIO.cpp' -Everything is commented out ????
                'src/CharmmEnergy.h', 'src/CharmmEnergy.cpp' -Added Lazaridis EEF1 implicit solvation function EEF1Ener
                'src/SasaCalculator.h', 'src/SasaCalculator.cpp' -Minor changes
                'src/AtomContainer.h', 'src/AtomContainer.cpp' -Added support for the atomId
                'src/File.h', 'src/File.cpp' -Moved init function to private
                'src/SurfaceAreaAndVolume.cpp' -Removed cout; changed call to vdwParam of the charmm parameter object
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added EEF1 Lazaridis solvation energy
                'src/AtomSelection.h', 'src/AtomSelection.cpp' -Fixed bug, it was returning pointers to a local function variable.
                 Added an empty selection. Now also atoms are called to remove a selection flag when clearStoredSelection is called
                
                'src/EnergySet.h' -Removed unnecessary includes for the specific Interactions (they are all used as Interaction
                 objects anyway)
                'src/MslTools.h', 'src/MslTools.cpp' -Added atomId, identityId, positionId supporting functions. Now splitIntAndString
                 returns a bool
                'src/AtomPointerVector.h' -Nothing major, just a little cleanup
                'src/AtomicPairwiseEnergy.cpp', 'tests/testSurfaceAreaAndVolume.cpp' -Changed call to vdwParam of the charmm parameter
                 object
                'examples/example_SasaCalculator_usage.cpp' -An example file for calculating SASA
                'examples/examples.mk' -Added example_SasaCalculator_usage.cpp
                'Makefile' -Change env variables to MSL_BOOST MSL_GLS etc. Added various objects and tests
0.6.1.0    March 04, 2010    dwkulp
                'tests/testSasaCalculator.cpp' -NO printResidueSasaTable in sasa object
                'src/PythonMSL.cpp' -capping off of residue SASA such that normalized SASA is only from 0 to 1
                'src/MslTools.h', 'src/MslTools.cpp' -fileExists function added
                'src/EnergeticAnalysis.h', 'src/EnergeticAnalysis.cpp' -factored parameter file out, no more hardcoding of path
                 names
                'src/Residue.h', 'src/Residue.cpp' -updated comments and added a bit of code to findNeighbors ... uses a default
                 value and actually uses _atomInOtherResidue now
                'src/PDBReader.cpp' -will store SEG_ID, thought this had been already fixed, but I guess not
                'programs/analEnergy.cpp' -Uses a parameter file for EnergeticAnalysis
                'programs/getSurroundingResidues.h', 'programs/getSurroundingResidues.cpp' -new program to get surrounding residues
                 and align them
                'Makefile' -added new programs
0.6.0.1    February 24, 2010    brettth
                'examples/example_AtomContainer_usage.cpp', 'examples/example_multipleAtomsCoordinates.cpp', 'examples/example_multipleResidueIdentities.cpp'
                 -Adding using namespace MSL; to cpp files so that they will know about our new namespace.
0.6.0.0    February 23, 2010    dwkulp
                'src/Tree.h' -fix namespace bracket placement to compile on mac
                'programs/printSequence.cpp', 'programs/printSequence.h' -new program
0.5.4.0    February 22, 2010    brettth
                I have added namespace MSL to all of the MSL objects in /src.
                Also removed 'using namespace std;' from header files.
0.5.3.0    February 17, 2010    dwkulp
                'src/PythonMSL.cpp' -added getSasa function to python interface
                'exampleFiles/pymolrc.py' -PyMOL init file that includes helper functions for MSL python interface
                'Makefile' -Include testSasaCalculator test
0.5.2.0    February 16, 2010    dwkulp
                'tests/testPDBFragments.cpp' -API change constructor needs second string for BBQTable, empty string in this test
                 just so it compiles
                'tests/testRegEx.cpp' -Added generic use of Regular Expressions in MslTools
                'tests/testALNReader.cpp', 'tests/testData.h', 'src/ALNReader.h', 'src/ALNReader.cpp' -Added ClustalW aln file
                 readers and some test data
                'tests/testDerivatives.cpp' -pow function was ambigous and compile complained on MAC
                'src/MslTools.h', 'src/MslTools.cpp' -some pow functions were ambigous, added regex function
                'src/PDBWriter.cpp' -removed pymol if statment bug in writeREMARKs function
                'programs/energyTable.cpp' -bug effected by CharmmEnergy::setDielectric being removed
                'Makefile' -New testALNReader and ALNReader classes
0.5.1.0    February 14, 2010    asenes
                'src/AtomBondBuilder.h', 'src/AtomBondBuilder.cpp' -New object: it adds bonded information to the atoms based on
                 their distance
                'src/Transforms.h', 'src/Transforms.cpp' -Added functions for direct edits of protein degrees of freedom setBondDistance,
                 setBondAngle and setDihedral, requires that there is bonding information
                'src/Atom.h', 'src/Atom.cpp' -Added function to find all the atoms that are bonded (even through intermediary atoms)
                 to the atom. It takes exclusions, this way one can get all the atoms bonded in a certain branch
                'tests/testSasaCalculator.cpp' -Missing test for the SASA calculator that was not added previously
                'tests/testTransformBondAngleDiheEdits.cpp' -Test for the direct edit functions of protein degrees of freedom setBondDistance,
                 setBondAngle and setDihedral in Transforms
                'tests/testAtomBondBuilder.cpp' -Test for AtomBondBuilder, adds bonded information to the atoms based on their
                 distance
0.5.0.0    February 13, 2010    asenes
                'tests/testBBQ2.cpp', 'tests/testNonBondedCutoff.cpp', 'tests/testIcBuilding.cpp', 'tests/testCCD.cpp', 'tests/testSystemIcBuilding.cpp',
                 'tests/testCharmmBuild.cpp', 'tests/testBoost.cpp', 'tests/testEnvironmentDatabase.cpp', 'tests/testEnergySet.cpp',
                 'tests/testPDBIO.cpp', 'tests/testFrame.cpp', 'tests/testPDBFragments.cpp', 'tests/testPolymerSequence.cpp', 'tests/testBBQ.cpp',
                 'tests/testCharmmEnergies.cpp', 'tests/testCoiledCoils.cpp', 'tests/testSurfaceAreaAndVolume.cpp', 'tests/testSystemCopy.cpp',
                 'tests/testHelixGenerator.cpp', 'tests/testSymmetry.cpp', 'tests/testEnvironmentDescriptor.cpp', 'tests/testGenerateCrystalLattice.cpp',
                 'tests/testAtomSelection.cpp', 'tests/testTransforms.cpp', 'tests/testLoopOverResidues.cpp',
                 'programs/alignMolecules.cpp', 'programs/createFragmentDatabase.cpp', 'programs/grepSequence.cpp', 'programs/calculateSasa.cpp',
                 'programs/searchFragmentDatabase.cpp', 'programs/getSelection.cpp', 'programs/generateCoiledCoils.cpp', 'programs/getSphericalCoordinates.cpp',
                 'Makefile', 'src/PDBReader.cpp', 'src/PDBFragments.h', 'src/CrystalLattice.cpp', 'src/AtomicPairwiseEnergy.h',
                 'src/EnergeticAnalysis.cpp', 'src/CharmmSystemBuilder.cpp', 'src/System.h', 'src/Symmetry.cpp', 'src/Frame.cpp',
                 'src/PDBFragments.cpp', 'src/SystemRotamerLoader.cpp', 'src/Atom3DGrid.cpp', 'src/SasaCalculator.h', 'src/PrincipleComponentAnalysis.h',
                 'src/SasaCalculator.cpp', 'src/TwoBodyDistanceDependentPotentialTable.h', 'src/Position.h', 'src/AtomDihedralRelationship.h',
                 'src/AtomContainer.cpp', 'src/CoiledCoils.h', 'src/Chain.h', 'src/Transforms.h', 'src/AtomAngleRelationship.cpp',
                 'src/Minimizer.cpp', 'src/TwoBodyInteraction.h', 'src/HelixGenerator.h', 'src/FourBodyInteraction.h', 'src/AtomDistanceRelationship.cpp',
                 'src/BBQTable.cpp', 'src/CCD.cpp', 'src/SurfaceAreaAndVolume.h', 'src/AtomDistanceRelationship.h', 'src/BackRub.h',
                 'src/EnergySet.h', 'src/AtomGeometricRelationship.h', 'src/SurfaceAreaAndVolume.cpp', 'src/PDBReader.h', 'src/Atom3DGrid.h',
                 'src/BBQTable.h', 'src/AtomicPairwiseEnergy.cpp', 'src/TwoBodyDistanceDependentPotentialTable.cpp',
                 'src/CCD.h', 'src/CoiledCoils.cpp', 'src/AtomSelection.cpp', 'src/Frame.h', 'src/AtomGroup.h', 'src/BBQTableWriter.cpp',
                 'src/PrincipleComponentAnalysis.cpp', 'src/CrystalLattice.h', 'src/AtomAngleRelationship.h', 'src/Transforms.cpp',
                 'src/Residue.cpp', 'src/AtomContainer.h', 'src/PDBWriter.h', 'src/PolymerSequence.cpp', 'src/HelixGenerator.cpp',
                 'src/Minimizer.h', 'src/Residue.h', 'src/ThreeBodyInteraction.h', 'src/AtomSelection.h', 'src/EnvironmentDescriptor.h',
                 'src/Position.cpp', 'src/PyMolVisualization.h', 'src/PolymerSequence.h', 'src/PythonMSL.cpp', 'src/Quaternion.cpp',
                 'src/BBQTableReader.cpp', 'src/EnvironmentDescriptor.cpp', 'src/Selectable.h', 'src/System.cpp',
                 'src/Quaternion.h', 'src/Symmetry.h', 'src/Chain.cpp', 'src/EnvironmentDatabase.cpp', 'src/AtomDihedralRelationship.cpp',
                 'src/PSFReader.h', 'src/PDBWriter.cpp', 'src/HelixFusion.cpp' -Changed object name from AtomVector to AtomPointerVector
                 'src/AtomVector.h', 'src/AtomVector.cpp', 'tests/testAtomVector.cpp', -Removed
                 'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp', 'tests/testAtomPointerVector.cpp', -Added
0.4.1.2    February 10, 2010    sabs
                'src/Atom.cpp' -Added the call setSelectionFlag(all,true) to Atom::copy()
0.4.1.1    February 10, 2010    jedonald
                'src/Frame.cpp' -Add frames for Asp and Glu oxygens
                'programs/tableEnergies.cpp' -Change default potential (to a non
                'programs/getSphericalCoordinates.cpp' -Print Lys dihedral angles now instead of spherical coordinates
                'programs/generateCoiledCoils.cpp', 'programs/generateCoiledCoils.h' -Program for making coiled coils
0.4.1.0    February 05, 2010    sabs
                'src/SystemRotamerLoader.cpp', 'src/RotamerLibraryReader.cpp', 'src/RotamerLibrary.cpp', 'src/RotamerLibrary.h',
                 'src/RotamerLibraryReader.h', 'src/RotamerLibraryWriter.h', 'src/RotamerLibraryWriter.cpp', 'tests/testRotamerLibraryWriter.cpp'
                 -Added readFile and writeFile to RotamerLibrary to read and write a library file
0.4.0.3    February 04, 2010    dwkulp
                'tests/testAtomSelection.cpp', 'src/Atom.cpp' -all selection flag in Atom
0.4.0.2    February 03, 2010    dwkulp
                'Makefile' -PythonMSL compile for linux and mac
0.4.0.1    February 01, 2010    asenes
                'examples/example_AtomContainer_usage.cpp', 'examples/example_multipleAtomsCoordinates.cpp', 'examples/example_multipleResidueIdentities.cpp',
                 'exampleFiles/example0000.pdb', 'exampleFiles/example0001.pdb' -Files for the example programs
                'tests/testNonBondedCutoff.cpp' -Added function (buildNonBonded) to recalculate the non
0.4.0.0    February 01, 2010    asenes
                'src/EnergySet.h', 'src/EnergySet.cpp' -Added function to delete an energy term (for example, all VDW) resetTerm(string
                 _term)
                'src/MslTools.h', 'src/MslTools.cpp' -Added the getMSLversion function to print MSL version and date
                'src/AtomGroup.h', 'src/AtomGroup.cpp' -Added getGeometricCenter(unsigned int _stamp=0), needed to implement the
                 energy cutoffs
                'src/CharmmEnergy.h', 'src/CharmmEnergy.cpp' -The object is no longer in charge of the dielectric, e14factor, use
                
                'src/CharmmElectrostaticInteraction.h', 'src/CharmmElectrostaticInteraction.cpp', 'src/CharmmVdwInteraction.h',
                 'src/CharmmVdwInteraction.cpp' -Added support for distance cutoffs and switching function
                'src/Atom.h', 'src/Atom.cpp' -Added getGroupGeometricCenter groupDistance getBonds (needed to implement the energy
                 cutoffs). Also, replaced unecessary size_t with unsigned int (small change).
                'src/AtomicPairwiseEnergy.h', 'src/AtomicPairwiseEnergy.cpp' -Added variables, getters and setters for supporting
                 storage of the dielectric constant, usa
                'src/Position.h' -Small change in the toString(), an asterisk indentifies the active identity
                'src/AtomContainer.h', 'src/AtomContainer.cpp' -Many changes. Added readPdb and writePdb functions. Allowed option
                 A 7 CA in addition to the comma separated A, 7, CA for the exists(string) function and ( ) operator. Added addAtom
                 ( string _name, double _x, double _y, double _z ) . Changed getFoundAtom to getLastFoundAtom for consistency
                'src/System.h', 'src/System.cpp' -Added functions getEnergySummary and printEnergySummary. Also, the toString()
                 now print as a PolymerSequence
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp', 'bin/testNonBondedCutoff' -Added function (buildNonBonded)
                 to recalculate the non
                'programs/energyTable.cpp' -Fixed change of API in AtomicPairwiseEnergy setVdwScale function
                'examples/example_multipleAtomsCoordinates.cpp', 'examples/example_multipleResidueIdentities.cpp', 'examples/example_AtomContainer_usage.cpp',
                 'examples/examples.mk' -New collection of examples for the documentation. The examples.mk file is sourced by the
                 main Makefile
                'Makefile' -Added new class of programs EXAMPLES, in the examples directory. They are listed in a local make file
                 examples.mk that is source by the main Makefile
0.3.2.1    January 31, 2010    dwkulp
                'src/PairwiseEnergyCalculator.cpp', 'src/Selectable.h', 'src/System.h', 'src/Atom3DGrid.h', 'src/Tree.h', 'src/AtomVector.h',
                 'src/CharmmTopologyResidue.h', 'src/Matrix.h', 'src/MslTools.h', 'src/CharmmEnergy.cpp', 'src/CharmmEnergy.h',
                 'Makefile' -modifications for comilation under MacOS, consists of sys/types.h header file inclusion, template
                 class type instantation,static const definitions in cpp not h files
0.3.2.0    January 29, 2010    asenes
                'src/RotamerLibrary.h', 'src/RotamerLibrary.cpp' -addeed a bunch of new getter functions
                'src/RotamerLibraryBuilder.h', 'src/RotamerLibraryBuilder.cpp' -An object to create a RotamerLibrary starting from
                 the conformation of molecules
                'src/RotamerLibraryWriter.h', 'src/RotamerLibraryWriter.cpp' -An object to write a rotamer library
                'tests/testRotamerLibraryWriter.cpp' -A test for the RotamerLibraryWriter
                'Makefile' -Added RotamerLibraryWriter and RotamerLibraryBuilder and testRotamerLibraryWriter
0.3.1.1    January 29, 2010    asenes
                'tests/testPolymerSequence.cpp' -Changed to reflect the changes in PolymerSequence uploaded with the previous commit
                 (I forgot this file)
0.3.1.0    January 29, 2010    asenes
                'Makefile' -The debug compile option is now under the control of an environmental variable T. Set it to T if you
                 want to compile in debug mode
                'src/PolymerSequence.h', 'src/PolymerSequence.cpp' -The toString() option now outputs in the same format of the
                 input sequence; added a void setSequence(System &_sys) function; added void setSequence(const AtomVector &_atoms)
                 and relative constructor PolymerSequence(const AtomVector &_atoms); the input sequence no longer require to have
                 chains on multiple lines, one line for everything works; removed a hard
                'src/MslTools.h', 'src/MslTools.cpp' -Added functions to check if a string is only digits, alphanums or spaces
                 (bool isDigitChars(string _input) and isAlphaNumericChars, isAlphaChars, isWhiteSpaces); also added function to
                 get a random int between 0 and _max unsigned int getRandomInt(unsigned int _max)
0.3.0.2    January 27, 2010    brettth
                'Makefile' -Changing FFTW to use library in EXTERNAL_LIB_DIR.
0.3.0.1    January 23, 2010    asenes
                'Makefile' -Removed duplicated entries for optimizeLP and optimizeMC, programs that belong to special conditional
                 sections
0.3.0.0    January 23, 2010    asenes
                'src/PolymerSequence.cpp' -Added flexibility to sequence string, it now accepts both {34} LEU and {34}LEU, as well
                 as {34}[ILE LEU VAL] and [{34}ILE LEU VAL]
                'tests/testPolymerSequence.cpp' -updated for the last changes of PolymerSequence
                'src/EnergySet.h', 'src/EnergySet.cpp' -TWO CHANGES: 1) API CHANGE to fix bug, now calcEnergy(false) is not allowed
                 to calculated all energy terms, including those that refer to inactive atoms: use instead calcEnergyAllAtoms();
                 2) added ability to save list of interactions with saveEnergySubset(string _subsetName, string _selection1, string
                 _selection2), and calcEnergyOfSubset(string _subsetName)
                'src/Chain.cpp' -added some comments
                'src/MslTools.cpp' -added some comments
                'src/System.h', 'src/System.cpp' -made changes to reflect the API change in the EnergySet
                'Makefile' -removed unneded extensions from program names in PROGRAMS section
0.2.1.0    January 19, 2010    brettth
                'src/CartesianGeometry.cpp' -Fixing bug in getXRotationMatrix and getZRotationMatrix. They rotated CW instead of
                 CCW.
0.2.0.4    January 18, 2010    brettth
                'src/RandomSeqGenerator.cpp', 'src/RandomSeqGenerator.h' -RandomSeqGenerator object now holds the random number
                 generator as a member variable instead of creating a new one each time getRandomSeq is called. This ensures better
                 randomization.
                'src/PhiPsiStatistics.h', 'src/PhiPsiStatistics.cpp' -PhiPsiStatistics calls now pass in a const Residue reference
                 instead of a Reference. This way no copy operation is needed for a call. This was a huge overhead when calling
                 PhiPsi repeatedly. Also, new change to add functions where you ask for stats given a given phi/psi combo (i.e.,
                 don't need to pass in residue objects if you already know name, and phi/psi angles.)
0.2.0.3    January 05, 2010    jedonald
                'src/TwoBodyDistanceDependentPotentialTable.cpp' -Never count contacts within a side chain
                'src/Frame.cpp' -Updates to residue frames, especially LYS and HIS, and allow option for multiply defined atoms
                 to be moved together
                'src/Chain.cpp', 'src/Position.cpp', 'src/System.cpp' -Garbage collection updates (from Alessandro Senes, email)
                
                'src/Frame.h' -Updates to residue frames, especially LYS and HIS, and allow option for multiply defined atoms to
                 be moved together
                'src/Residue.h' -Allow option for multiply defined atoms to be found in residue based search
                'src/Quench.cpp' -Assign coordinates (had been accidentally removed)
                'src/Residue.cpp' -Allow option for multiply defined atoms to be found in residue based search
                'programs/getSphericalCoordinates.cpp' -Allow negative residues, always print plane angles
                'programs/getSphericalCoordinates.h' -Allow negative residues, always print plane angles
                'programs/runQuench.cpp' -Use largeRotNum correctly
0.2.0.2    January 04, 2010    brettth
                'src/RandomSeqGenerator.cpp', 'src/RandomSeqGenerator.h' -A simple class which will generate random sequences with
                 a given probability for each character. I use it to generate random peptide sequences.
                'tests/testRandomSeqGenerator.cpp' -A test of the RandomSeqGenerator code.
                'Makefile' -Updating Makefile to include RandomSeqGenerator code.
0.2.0.1    December 30, 2009    brettth
                'src/HelixGenerator.cpp', 'src/HelixGenerator.h' -Adding ability to center or not center the helix. Also adding
                 code to find the helical axis.
                'Makefile' -Adding ability to use fftw. I currently use this in my SinWaveFitter program, but other people may
                 want to use it for other purposes. By default, it still assumes that you don't have it installed.
0.2.0.0    December 17, 2009    dwkulp
                'programs/getSelection.cpp', 'programs/getSelection.h' -add an outPdb option
                'src/Frame.h', 'src/Frame.cpp' -compute frame from 2 orthogonal lines
                'src/Hash.h' -google hash map added
                'src/Line.cpp' -bug in projection function fixed
                'src/Matrix.cpp', 'src/Matrix.h' -added BOOST serialization code, so we can write binary matrix objects out to
                 file
                'src/PDBFormat.cpp' -rearrange reading fields
                'src/PDBReader.h' -needed to initialize scaleTranslation, scaleRotation matrices
                'src/Transforms.h' -moved align and orient using CartesianPoints to public members
                'src/Residue.h', 'src/Residue.cpp' -added additional findNeighbors, to find a neighboring residue by any specified
                 atom type or all atoms
0.1.2.0    November 04, 2009    asenes
                'Makefile' -Using environmental variables (with defaults) to set a number of preferences
                'src/TwoBodyInteraction.h', 'src/FourBodyInteraction.h', 'src/ThreeBodyInteraction.h' -Fixed bug that was creating
                 memory leak, the distructor must be virtual for the inheritance to work properly
                'src/CharmmTopologyResidue.cpp', 'src/FourBodyInteraction.h', 'src/ThreeBodyInteraction.h' -Fixed bug that was
                 creating memory leak, the distructor was not calling the pointer deletion
                'src/SystemRotamerLoader.h', 'src/SystemRotamerLoader.cpp' -Added functions addRotamers(...), they add rotamers,
                 load rotamers delete the old ones by default. Added also a flag to loadRotamers to preserve the old rotamers (addRotamers
                 is actually a wrapper for loadRotamers with the flag on)
                'src/EnergySet.cpp' -Fixed bug that was creating a memory leak. An unnecessary break was cutting pointer deletion
                 short in deletePointers
                'src/LogicalParser.cpp', 'src/Tree.h' -Created temporary fix for memory leak, now the Tree distroyes something
                 created by the LogicalParser, need to go back and revise
                'programs/calculateSasa.cpp' -Given the option
0.1.1.7    November 02, 2009    brettth
                'src/RotamerLibrary.cpp', 'src/CharmmTopologyResidue.cpp' -These files need to include <stdio.h> in order to use
                 sprintf. Fixing this.
0.1.1.6    October 16, 2009    dwkulp
                -moved OptionParser calls, such that readArgv happens after the object is set up
0.1.1.5    October 14, 2009    jedonald
                -Program to calculate the energy of a protein using an contact potential
                -Fix a typo in header file
                -Change default file path
                -Modify initial setup for KB potentials to only use KB energy, also pass by reference now the tables and vectors
                
                -Add function to calculate KB total energy
                -Add functions to calculate total energy and to set up for quench using KB energies. Also add bool to allow local
                 sidechain backbone contacts or not as an option. Defaults to only non
                -Add checked in programs to current Makefile
0.1.1.4    October 13, 2009    brettth
                -Adding code that will generate idealized helices.
                -A test of the Helix Generator code.
0.1.1.3    October 12, 2009    jedonald
                -Additional code to allow potential tables to be used in quenching.
0.1.1.2    October 12, 2009    jedonald
                -Fix compilation warnings
                -Allow potential tables to be used in quenching. Involved moving some include statements
0.1.1.1    October 09, 2009    dwkulp
                -toString(), when printing pymol lines have axis length be a parameter
                -toString(), variable axis length
                -Fixed get frame for LYS, center on NZ , no translation
0.1.0.0    October 09, 2009    dwkulp
                -added MYHEADERS, seems we could handle this in Makefile from MYPROGS
                -added BBQTable arguement to localSamplingPDB, therefore it can be passed in from python
                -output change to add total number of rotamers
                -ifdef for compile was SELECTION, not ATOMSELECTION
                -Residue version of AtomSelection
                -ASN side chain oxygen atom is OD1 not OD2; print error out when atoms don't exist
                -Constructing frames from atoms generated left
                -Made Residue's selectable so that we can use ResidueSelection on them
                -added BBQTable arguement
                -added BBQTable arguement to cpp file as well!
                -include Atom.h
                -printOptions, will print out all possible options
                -phiPsi Statistics is now optional and not the default, only getDihedrals for amino acids
                -moved parameters around sine createSystem sets parameter defaults
                -when GMEC found after DEE, don't run MC, just print out structure.
                -fixed usage statement to include positions
                -moved readArgv calls to after setRequired/setAllowed, this enables the 'printOptions' system to work inside OptionParser
                
                -create a winning PDB structure if option structureConfig is set
                -new programs, some reason we now need -lpthread to link our programs?
                -new program to evaluate either an atom or residue selection and print out what it finds
0.0.2.2    October 08, 2009    brettth
                -Taking out a blank line, just to test submit procedure.
0.0.2.1    October 08, 2009    brettth
                -Modified PhiPsiStatistics so that it no longer assumes a 5 degree increment between Phi/Psi bins. Instead, it
                 looks at the increment of the Psi for the first and second entry in the table, and uses that.
                -There was a slight bug when filling in missing atoms for a chain. The new backbone atoms for the last 2 residues
                 were not being updated correctly.
0.0.2.0   10/01/2009 Alessandro Senes, merged after month local repository, many changes (no API incompatible changes, IIRC)
          ADDED: 
                   src/Atom3DGrid.h                             Object that creates a 3D grid (based on a given size) and distributes atoms in the grid cells.  Useful for finding lists of neighbors within some distance cutoff
                   src/Atom3DGrid.cpp                           Object that creates a 3D grid (based on a given size) and distributes atoms in the grid cells.  Useful for finding lists of neighbors within some distance cutoff
                   src/SasaAtom.h                               A container for calculating the SASA of an atom (uses SurfaceSphere)
                   src/SasaAtom.cpp                             A container for calculating the SASA of an atom (uses SurfaceSphere)
                   src/SasaCalculator.h                         An object that calculates the SASA of a AtomVector (uses SasaAtom)
                   src/SasaCalculator.cpp                       An object that calculates the SASA of a AtomVector (uses SasaAtom)
                   src/SurfaceSphere.h                          An almost homogeneously spaced sphere of points around a center (an atom), used to calculate the SASA
                   src/SurfaceSphere.cpp                        An almost homogeneously spaced sphere of points around a center (an atom), used to calculate the SASA
                   programs/alignMolecules.cpp                  RMSD alignment of two molecules based on selections
                   programs/calculateSasa.cpp                   Calculates the SASA of a PDB
		   myProgs/myProgs.mk.RENAME_ME                 The make file for programs/objects that do not belong to the repository.  Rename necessary so that people do not commit their own modified version
          CHANGED:
                   src/Atom.h                                   Added radius and sasa as properties and relative functions, fixed a bug with the group number
                   src/Atom.cpp                                 Added radius and sasa as properties and relative functions, fixed a bug with the group number
                   src/AtomGroup.cpp                            Fixed a bug with the group number
                   src/AtomVector.h                             Added new constructor for: AtomVector a(22, &atom) add 22 atoms with the same pointer value
                   src/AtomVector.cpp                           Added new constructor for: AtomVector a(22, &atom) add 22 atoms with the same pointer value
                   src/CartesianGeometry.cpp                    Just removed some commented out code
                   src/CCD.h                                    Added #include<algorithm> and moved some includes from .cpp to .h file
                   src/CCD.cpp                                  Added #include<algorithm> and moved some includes from .cpp to .h file
                   src/Chain.h                                  Added renumberChain function, removed commented out code, fixed bug in updatePositionMap
                   src/Chain.cpp                                Added renumberChain function, removed commented out code, fixed bug in updatePositionMap
                   src/CharmmAngleInteraction.cpp               Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmBondInteraction.cpp                Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmDihedralInteraction.h              Reverted small typo in toString. Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmDihedralInteraction.cpp            Reverted small typo in toString. Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmElectrostaticInteraction.cpp       Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmImproperInteraction.cpp            Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmParameterReader.h                  Added operator=, removed commented out code, moved include <math.h> from .cpp to .h, fixed bug in copy constructor and reset (vdwParamPairMap not assigned)
                   src/CharmmParameterReader.cpp                Added operator=, removed commented out code, moved include <math.h> from .cpp to .h, fixed bug in copy constructor and reset (vdwParamPairMap not assigned)
                   src/CharmmSystemBuilder.h                    Fixed copy constructor and added operator=
                   src/CharmmSystemBuilder.cpp                  Fixed copy constructor and added operator=
                   src/CharmmTopologyReader.h                   Fixed copy constructor and added operator=
                   src/CharmmTopologyReader.cpp                 Fixed copy constructor and added operator=
                   src/CharmmTopologyResidue.h                  Fixed copy constructor and added operator=, added reset function
                   src/CharmmTopologyResidue.cpp                Fixed copy constructor and added operator=, added reset function
                   src/CharmmUreyBradleyInteraction.cpp         Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CharmmVdwInteraction.cpp                 Fixed bug for compilation with 32bit: added (Atom*) declaration for NULL pointer vector<Atom*>(<num>, (Atom*)NULL);
                   src/CoiledCoils.h                            Moved include from .cpp to .h and added comment warning for not having a proper copy constructor
                   src/CoiledCoils.cpp                          Moved include from .cpp to .h and added comment warning for not having a proper copy constructor
                   src/EnergySet.h                              Commented out underfined getTotalNumberOfInteractions(unsigned int _type) function, added setAllTermsInactive and setAllTermsActive functions, added check for atom having coordinates in calculateEnergy, fixed bug in getTotalNumberOfInteractions
                   src/EnergySet.cpp                            Commented out underfined getTotalNumberOfInteractions(unsigned int _type) function, added setAllTermsInactive and setAllTermsActive functions, added check for atom having coordinates in calculateEnergy, fixed bug in getTotalNumberOfInteractions
                   src/File.h                                   Moved include from .cpp to .h, added comment ("NOTE: SHOULD THE ENUMS BE IN UPPERCASE?")
                   src/File.cpp                                 Moved include from .cpp to .h, added comment ("NOTE: SHOULD THE ENUMS BE IN UPPERCASE?")
                   src/Interaction.h                            Added function atomsHaveCoordinates()
                   Makefile                                     Moved source CCD BackRub Quench SurfaceAreaAndVolume to GSL section, moved RegEx source to BOOST section, added new programs and libraries
                   src/Matrix.h                                 Added << operator, change toString() to constant function
                   src/Matrix.cpp                               Added << operator, change toString() to constant function
                   src/MonteCarloOptimization.h                 Commented out unused linkedPositions variable and linkPositions function, moved includes to .h from .cpp
                   src/MonteCarloOptimization.cpp               Commented out unused linkedPositions variable and linkPositions function, moved includes to .h from .cpp
                   src/MslTools.h                               Added quickSort(vector<double>&) and quickSortWithIndex(vector<double>&, vector<unsigned int>&) functions, added check for __GSL__ and alternate code for random number generator in getRandomAlphaNumString, moved includes from .cpp to .h
                   src/MslTools.cpp                             Added quickSort(vector<double>&) and quickSortWithIndex(vector<double>&, vector<unsigned int>&) functions, added check for __GSL__ and alternate code for random number generator in getRandomAlphaNumString, moved includes from .cpp to .h
                   src/OptionParser.cpp                         Changed hardcoded option "version" to "mslVersion"
                   src/Position.h                               Added setResidueName function, added toString, renumberNoUpdate, getSasa functions, removed commented out code, added break in loop in updateResidueMap (trivial speedup)
                   src/Position.cpp                             Added setResidueName function, added toString, renumberNoUpdate, getSasa functions, removed commented out code, added break in loop in updateResidueMap (trivial speedup)
                   src/Residue.h                                Added getAtomMap function and << operator, converted toString to constant function, added getSasa function
                   src/Residue.cpp                              Added getAtomMap function and << operator, converted toString to constant function, added getSasa function                             
                   src/RotamerLibrary.cpp                       Removed commented out code
                   src/SelfPairManager.cpp                      Fixed bug for compilation with 32bit: added (Residue*) declaration for NULL pointer vector<Residue*>(<num>, (Residue*)NULL);
                   src/PDBFormat.cpp                            Silenced warning for non-numerical in charge field of PDB (often corrupted but not so important).  Should we use a verbose flag here?
                   src/PDBReader.h                              Added getMissingAtoms and getMissingResidues functions.  In read(), added check string length for substr calls, and added code for getting missing atom and residues from the remarks (not yet using PDBFormat! change!)
                   src/PDBReader.cpp                            Added getMissingAtoms and getMissingResidues functions.  In read(), added check string length for substr calls, and added code for getting missing atom and residues from the remarks (not yet using PDBFormat! change!)
                   src/PDBWriter.cpp                            Added warnings for failed writeln call in the write function.  Changed slighly default REMARK line (should we remove this?), commented on hard-coded if (singleRemarks[i].find("from pymol") == string::npos){ (what is it?)
                   src/System.cpp                               Rationalized code in assignCoordinates
                   src/Transforms.h                             Added lastTranslation CartesianPoint.  Added getLastRotationMatrix and getLastTranslation functions.  Updated many transformation functions.
                   src/Transforms.cpp                           Added lastTranslation CartesianPoint.  Added getLastRotationMatrix and getLastTranslation functions.  Updated many transformation functions.
*/


#endif
