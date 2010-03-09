/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#define MSLVERSION "0.7.0.0"
#define MSLDATE "March 09, 2010"

/*
HISTORY:
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
