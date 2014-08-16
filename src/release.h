/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

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

#define MSLVERSION "1.2.2.5"
#define MSLDATE "August 16, 2014"

/*
HISTORY:
1.2.2.5    August 16, 2014    grigoryan
                'src/Chain.cpp', 'src/Chain.h', 'src/PDBReader.cpp', 'src/System.cpp', 'src/System.h' -Added the ability to preserve
                 residue order (even if residue ID's are out of order) in building a System. Also, changed the behavior of PDBReader,
                 so that if the chain ID is empty, it uses the segment ID as the chain ID for the purpose of reading/parsing and
                 deciding which residues go into which chain. Upon writing, of course, it is just the first character of segment
                 ID that becomes the chain ID in the output file. This behavior is very useful for reading in very large PDB files
                 where the nuber of unique chains exceeds the number of chars.
                'src/Chain.cpp', 'src/Chain.h', 'src/PDBReader.cpp', 'src/System.cpp', 'src/System.h' -Added the ability to preserve
                 residue order (even if residue ID's are out of order) in building a System. Also, changed the behavior of PDBReader,
                 so that if the chain ID is empty, it uses the segment ID as the chain ID for the purpose of reading/parsing and
                 deciding which residues go into which chain. Upon writing, of course, it is just the first character of segment
                 ID that becomes the chain ID in the output file. This behavior is very useful for reading in very large PDB files
                 where the nuber of unique chains exceeds the number of chars.
1.2.2.4    May 29, 2014    bkmueller
                'src/SelfPairManager.cpp', 'src/SelfPairManager.h' -Added recalculateNonSavedEnergies which allows the user to
                 supply a mask which will not recalculate the given pair energies, useful for skipping intramolecular pair energies
                
1.2.2.3    May 29, 2014    sgfc
                'src/SelfPairManager.cpp' -Added
1.2.2.2    May 26, 2014    sgfc
                './src/SelfPairManager.h' -Added a mask to runGreedyOptimizer to exclude particular rotamers/identities from analysis
                
1.2.2.1    February 18, 2014    asenes
               WARNING!  Files submitted without building tree or running tests.
                'RELEASE_NOTES.txt', 'var/header.txt' -Updated year in copyright statement to 2014
1.2.2.0    February 18, 2014    scraven
                'src/RandomNumberGenerator.cpp', 'src/RandomNumberGenerator.h' -added getRandomOrder function to return vector
                 of randomly ordered unsigned integers
1.2.1.8    February 11, 2014    asenes
               WARNING!  Files submitted without building tree or running tests.
                'README.txt' -Updated contributor list and year for copyright
1.2.1.7    February 10, 2014    sgfc
                'src/SystemRotamerLoader.cpp' -Fixed bug in function loadRotamers(string _positionId, string _resName, string _levelName,
                 string _rotLib, bool _keepOldRotamers) where _resName and _rotLib were swapped
1.2.1.6    October 25, 2013    asenes
                'programs/alignMolecules.cpp' -Added option to write all models for a multi-model PDB input file (the default writes
                 only 1 model)
1.2.1.5    September 27, 2013    asenes
               WARNING!  Files submitted without building tree or running tests.
                'Makefile' -The default make, i.e. make all, now does not compile the sandbox and the LEAD programs
1.2.1.4    September 27, 2013    bhannigan
               WARNING!  Files submitted without building tree or running tests.
                'scripts/submit.py' -Fixing a bug which was causing us to miss submitting some files that were meant to be added.
                
1.2.1.3    September 27, 2013    bhannigan
               WARNING!  Files submitted without building tree or running tests.
                'scripts/TEST_WHITELIST' -Adding a file that we can use to exempt certain files from running tests during submit.
                
1.2.1.2    September 27, 2013    bhannigan
               WARNING!  Files submitted without building tree or running tests.
                'scripts/TEST_WHITELIST' -Adding a file that we can use to exempt certain files from running tests during submit.
                
1.2.1.1    September 27, 2013    bhannigan
               WARNING!  Files submitted without building tree or running tests.
                'scripts/TEST_WHITELIST' -Adding a file that we can use to exempt certain files from running tests during submit.
                
                'scripts/submit.py' -Adding flag to allow users to skip build when submitting files. This can be dangerous, so
                 it should be used sparingly.
                'scripts/submit.py' -Adding flag to allow users to skip build when submitting files. This can be dangerous, so
                 it should be used sparingly.
                'scripts/TEST_WHITELIST' -Adding a file that we can use to exempt certain files from running tests during submit.
                
1.2.1.0    September 25, 2013    asenes
                'RELEASE_NOTES.txt', 'README.txt' -Moving the trunk to development toward v.1.3.
1.1.2.16    September 25, 2013    asenes
                'src/SelfPairManager.cpp', 'src/SelfPairManager.h' -Fixed bug in runLP, which was causing a compilation error.
                 It was returning the incorrect type vector<int> which no longer matched the underlying LinearProgramingOptimization
                 function return variable
                'Makefile' -Updated header and removed unused OpenNMP compiling option
1.1.2.15    September 18, 2013    asenes
                'README.txt' -Updated readme file (did not go in previous submit)
1.1.2.14    September 18, 2013    asenes
                'src/Atom.cpp' -Added check to setUnboundFromAll, if the bond list is emtpy returns immediately for speed
                'RELEASE_NOTES.txt' -Updated readme files
1.1.2.13    September 17, 2013    james
                'src/System.cpp' -Fix System::deletePointers() by moving the chain deletion to the end
1.1.2.12    September 17, 2013    james
                'src/System.cpp' -Fix System::deletePointers() by moving the chain deletion to the end
1.1.2.11    September 11, 2013    james
                'src/CharmmSystemBuilder.cpp' -Fixed bug that was passing incorrect cutoff information in the updateNonBonded function
                
1.1.2.10    August 06, 2013    jedonald
                'src/StrideReader.h', 'src/StrideReader.cpp', 'src/DSSPReader.h', 'src/DSSPReader.cpp' -Secondary structure file
                 readers
                'src/Quaternion.cpp' -Remove unused variables for compilation warning
                'src/SurfaceAreaAndVolume.cpp' -Comment out unused firstAngle, lastAngle
                'src/Matrix.cpp' -Remove unused x,y,z,a,eval_i
                'src/OnTheFlyManager.cpp' -Remove or comment unused tid variables
                'src/CharmmTopologyResidue.cpp' -Remove unused iterator pos
                'src/CCD.cpp' -Remove unused convergedDist variable
                'src/MslTools.cpp' -Remove unused absolutePath variable
                'src/ALNReader.cpp' -Modify regular expression searches to match more MSA output formats
                'src/LogicalParser.cpp' -Remove unused numResults variable
                'src/GSLMinimizer.cpp' -Comment out converged variable
                'src/PDBFragments.cpp' -Comment out unused validTriplet variable
                'myProgs/jedonald/betaBetaBakerChirality.cpp', 'myProgs/jedonald/jedonald.mk' -Program for determining beta hairpin
                 chirality (Koga et al, Nature 2012)
                'Makefile' -Have -fopenmp depend on MSL_OPENMP
                'scripts/mslBuildTools.py' -Add quotes to failure tag when there is an error
1.1.2.9    June 15, 2013    bhannigan
                'scripts/submit.py' -Slight bug that placed spaces inbetween letters in the submit comments. Oops.
1.1.2.8    June 15, 2013    bhannigan
                'scripts/submit.py' -U p d a t i n g s u b m i t s c r i p t t o h a n d l e n e w S o u r c e f o r g e r e p
                 o s i t o r y i n f o . A l s o , m e s s a g e s c a n n o w h a v e d a s h e s i n c o m m e n t s . S e e
                 : - t h i s - n o w - w o r k s .
1.1.2.7    June 11, 2013    sabs
                'myProgs/sabs/tmRulesCreator.cpp' -Takes a geometryFile instead of one geometry string.
                'myProgs/sabs/sabs.mk' -Added New programs
                'myProgs/sabs/designViaSequenceMC.cpp' -Print out accepted PDBs with sequence number
                'src/EnergySet.h', 'src/EnergySet.cpp' -Specify interactions to delete by type.
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Greedy optimizer goes through positions in random order.
                'src/System.h', 'src/Chain.h', 'src/Position.h' -defineRotamerSamplingLevels takes a reference
1.1.2.6    May 13, 2013    sabs
                'myProgs/sabs/CAHTM.cpp' -Removed position dependent Rotamer loading. This is CAHTM v.0.0.10.
1.1.2.5    April 30, 2013    sabs
                'src/SelfConsistentMeanField.cpp', 'src/SelfConsistentMeanField.h', 'src/MslTools.cpp', 'src/MslTools.h', 'src/MonteCarloManager.cpp',
                 'src/MonteCarloManager.h' -Created a const variable in MslTools for the molar gas constant. Updated all objects
                 that use this constant
1.1.2.4    April 29, 2013    sabs
                'src/MslTools.cpp', 'src/MslTools.h' -Methods to compute Boltzmann energies
                'src/FormatConverter.cpp' -Updated PDB3 conversions
                'programs/createEnergyTable.cpp' -The system needs to specified before reading parameter files in the HydrogenBondBuilder.
                
1.1.2.3    April 24, 2013    asenes
                'programs/alignMolecules.cpp' -Added support for specifying a model for NMR files with
1.1.2.2    April 22, 2013    sabs
                'Makefile', 'src/PDBWriter.cpp', 'src/FormatConverter.cpp' -Bug fixes in FormatConverter and PDBWriter. Added pdb2crd
                 to the Makefile
1.1.2.1    April 16, 2013    asenes
                'programs/pdb2crd.h' -A program for converting a PDB file to CHARMM format (header file, forgotten in previous
                 commit
1.1.2.0    April 16, 2013    asenes
                'src/FormatConverter.h', 'src/FormatConverter.cpp' -Improved conversion, added support for current PDB format (v.3),
                 added a function for converting a PDBFormat::atomData line
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added buildSystemFromPDB from AtomPointerVector, not
                 just from file
                'src/PDBWriter.h', 'src/PDBWriter.cpp' -Changed the way it integrates with the FormatConverted for changing the
                 atom/residue names, now it has a setter setConvertFormat to set to change for example CHARMM22 to PDB3
                'myProgs/sabs/CAHTM.cpp' -Tweaked for the changes in PDBWriter API for format conversion
                'programs/pdb2crd.cpp' -A program for converting a PDB file to CHARMM format
                'RELEASE_NOTES.txt' -Added the improved format conversion and the pdb2crd program to the release notes
1.1.1.10    April 14, 2013    asenes
                'src/PDBFragments.cpp', 'src/SurfaceAreaAndVolume.cpp', 'src/EnergeticAnalysis.cpp', 'src/Line.cpp', 'src/BackRub.cpp',
                 'src/CartesianGeometry.cpp', 'tests/gold/testCharmmBuild.cpp', 'tests/gold/testCharmmEnergies.cpp', 'tests/gold/testEZpotential.cpp',
                 'tests/gold/testRMSDalignment.cpp' -Fixed bug with abs(int) function used with double variables
                'testLevels/submit.level' -Restored testEZpotential testRMSDalignment in the list of GOLD tests since they work
                
                'scripts/submit.py', 'scripts/mslBuildTools.py' -Turned on GSL by default for submission tests and added code to
                 print the exception in case submit.py fails
1.1.1.9    April 14, 2013    asenes
                'programs/renumberResidues.cpp' -Fixed bug, abs(unsigned int) does not work on some compilers
1.1.1.8    April 13, 2013    asenes
                'src/Clustering.cpp' -Fixed bug, was using integer function abs instead of double function fab in double comparision
                
1.1.1.7    April 03, 2013    asenes
                'README.txt', 'RELEASE_NOTES.txt' -Updated in preparation for release v.1.2
1.1.1.6    March 23, 2013    asenes
                'src/SysEnv.cpp' -Corrected the location of MSL_PDB_TOP and added an variable to the file with hydrogens, MSL_PDBH_TOP.
                
                'src/PDBTopologyBuilder.cpp' -Fixed small bug in buildSystemFromPDB, if an unknow residues occurred it would fail
                 to build but it would not return a false flag
                'toppar/top_pdb2.3_H.inp', 'toppar/top_pdb2.3_noH.inp' -Added a RESI entry for water molecules
1.1.1.5    March 22, 2013    jedonald
                'src/RandomNumberGenerator.h', 'src/RandomNumberGenerator.cpp' -Allow for random integers to not include upper
                 limit (corrected version)
1.1.1.4    March 22, 2013    jedonald
                'src/Quench.cpp' -Fix bug related to shuffling repacking order (second change)
1.1.1.3    March 22, 2013    jedonald
                'src/Position.cpp' -Check for hiddenIdentityIndeces vector values, otherwise will segfault
                'src/RandomNumberGenerator.h', 'src/RandomNumberGenerator.cpp' -Allow for random integers to not include upper
                 limit
                'src/Quench.cpp' -Fix bug related to shuffling repacking order
1.1.1.2    March 20, 2013    sabs
                'myProgs/sabs/connectWithFragments.cpp' -A program to stitch the TM and CC domains of FtsB using fragments from
                 real proteins.
1.1.1.1    March 19, 2013    sabs
                'myProgs/sabs/predictHelixOligomer.cpp', 'myProgs/sabs/filterOligomerByConstraints.cpp' -Program filterOligomerBYConstraints
                 filters oligomer geometries based on hydrogen bond constraints. Program predictHelixOligomer takes a list of geometries
                 and repacks a sequence onto each geometry in order to predict the oligomeric state.
1.1.1.0    February 28, 2013    asenes
                'README.txt', 'var/header.txt' -Converted trunk to development version past 1.1 and toward 1.2. The naming conventions
                 for the trunk is now 1.1.1.x. The branch 1.1 was created in beta (bug fixing) with v.1.0.3.1 toward a 1.1.0.0
                 release after the code is stabilized.
1.0.3.1    February 28, 2013    sabs
                'src/CharmmSystemBuilder.cpp' -Solvation terms should not be built for 1-2 or 1-3 bonded atoms. Broken in version
                 1.0.2.4, fixed by this commit.
1.0.3.0    February 28, 2013    sabs
                'src/Residue.h', 'src/AtomGroup.cpp', 'src/Residue.cpp', 'src/Atom.h', 'src/AtomGroup.h', 'src/Position.cpp', 'src/Atom.cpp',
                 'src/Position.h', 'src/SelfPairManager.cpp' -Added support to hide identities. Fixed a bug in hideAtomConformation.
                
                'programs/repackSideChains.cpp', 'programs/repackSideChains.h' -Added solvation options.
                'myProgs/sabs/designViaSequenceMC.cpp' -A program that performs design by a monte carlo search in the sequence
                 space.
1.0.2.4    February 14, 2013    sabs
                'src/BaselineEnergyBuilder.cpp' -Facility to add identity specific baselines
                'src/CharmmEnergy.h', 'src/CharmmSystemBuilder.cpp', 'src/CharmmEnergy.cpp', 'src/CharmmSystemBuilder.h' -Added
                 IMM1 solvation model
                'src/CharmmIMM1Interaction.h', 'src/CharmmIMM1Interaction.cpp', 'src/CharmmIMM1RefInteraction.h', 'src/CharmmIMM1RefInteraction.cpp'
                 -Added IMM1 solvation model
                'programs/energyOptimizations.h' -Added solvation parameters
                'myProgs/sabs/sabs.mk', 'Makefile' -Makefile changes
                'myProgs/sabs/genHomoUniverse.cpp', 'myProgs/sabs/CAHTM.cpp' -Changes to my programs.
                'src/SysEnv.cpp' -Changed MSL_PDB_2.3_DOF to MSL_PDB_2_3_DOF since . can not be a part of csh env variable name
                
1.0.2.3    February 05, 2013    dwkulp
                'src/LinearProgrammingOptimization.h', 'src/LinearProgrammingOptimization.cpp', 'programs/optimizeLP.cpp', 'programs/optimizeMC.cpp'
                 -Updated API to use System setRotamerState
                'programs/renumberResidues.cpp', 'programs/renumberResidues.h' -Renumber based on fasta file
                'src/CharmmEnergyCalculator.cpp' -Zero out terms
                'src/MslExceptions.h' -AtomNotExist and RmsdAlignment Exceptions added
                'src/MslTools.h', 'src/MslTools.cpp' -parseMutationId function and new regex function
                'src/PDBFragments.h', 'src/PDBFragments.cpp' -Added a search for 2 fragments from same PDB file
                'src/PolymerSequence.h', 'src/PolymerSequence.cpp' -Added a PDBNamesFlag for HIS residues
                'src/SasaCalculator.h' -MslNotFoundException when no radii found for atom
                'src/SysEnv.cpp' -Added SaltBridge propensities and Prosite dataset
                'src/System.h', 'src/System.cpp' -assignCoordinates function can take a name-conversion map; atomExists funciton
                 modification to skip chain level if desired
                'myProgs/dwkulp/compareRosettaModels.cpp' -Resfile no longer required
                'myProgs/dwkulp/designCheck.cpp', 'myProgs/dwkulp/designCheck.h' -Added reportEnergyMetrics function to report
                 a per-position energy analysis
                'myProgs/dwkulp/resurfaceSaltBridges.cpp', 'myProgs/dwkulp/resurfaceSaltBridges.h', 'myProgs/dwkulp/findMotifMSA.cpp',
                 'myProgs/dwkulp/findMotifMSA.h', 'myProgs/dwkulp/dwkulp.mk' -New programs: resurfaceSaltBridges, findMotifMSA
                
                'myProgs/dwkulp/refinePotentialFusions.cpp', 'myProgs/dwkulp/refinePotentialFusions.h' -added PDBFragment searching
                
                'myProgs/dwkulp/setupRosettaMSA.cpp' -Added a trim for PDB file name inside MSA
                'myProgs/dwkulp/superRotamerExtraction.cpp', 'myProgs/dwkulp/superRotamerExtraction.h' -Added neighbor distance
                 arguement
                'tables/prosite.dat', 'tables/sb_prop_table_12_14_12.txt' -Adding Prosite database and Salt-Bridge Propensity database
                
1.0.2.2    December 12, 2012    dwkulp
                'src/PDBFormat.cpp' -fix to element extraction logic from atomname
1.0.2.1    December 12, 2012    dwkulp
                'programs/designSideChains.cpp' -updated use of energyOptimization API
                'programs/designSideChains.cpp' -updated use of energyOptimization API
1.0.2.0    December 12, 2012    dwkulp
                'Makefile' -CoiledCoilFitter RosettaScoredPDBReader Clustering, associated tests
                'examples/example_pdbfrag.cpp' -PDBFragment API change
                'myProgs/dwkulp/designLinearSegment.cpp' -reformat output pdb name
                'myProgs/dwkulp/dwkulp.mk' -added metalRotamers patchAnalysis compareRosettaModels setupRosettaFixbb getConditionalPhiPsi
                 refinePotentialFusions createTertFragDB buildRotamers getLoopLengths
                'myProgs/dwkulp/getVectorPairs.cpp', 'myProgs/dwkulp/getVectorPairs.h' -interChainOnly option
                'myProgs/dwkulp/insertLoopIntoTemplate.cpp', 'myProgs/dwkulp/insertLoopIntoTemplate.h' -fragmentChain option
                'myProgs/dwkulp/searchForFusions.cpp' -priority_queue, error checking code
                'myProgs/dwkulp/setupRosettaMSA.cpp', 'myProgs/dwkulp/setupRosettaMSA.h' -fixed some bugs with the way positions
                 are indexed
                'programs/createFragmentDatabase.cpp' -keep track of the total number of matches, write out options in output
                'programs/getSelection.cpp', 'programs/getSelection.h' -'list' option implemented
                'programs/optimizeMC.cpp', 'programs/optimizeLP.cpp', 'programs/energyOptimizations.h' -Fixed createSystem API
                 change from last check in, now still supports structureConfiguration + system as input
                'src/Chain.h' -Exception handling in getPosition
                'src/Predicate.h', 'src/LogicalParser.cpp', 'src/LogicalParser.h', 'src/Tree.h' -modifications to parsing boolean
                 expressions, almost all functional
                'src/MslExceptions.h' -MslNotFoundException implemented
                'src/MslTools.cpp', 'src/MslTools.h' -hamming_distance between 2 strings implemented
                'src/PDBFragments.cpp', 'src/PDBFragments.h' -MSLOUT functionality
                'src/RegEx.cpp' -Initialized 'stype' variable in constructor
                'src/PDBFormat.cpp' -Fixed the missing element issue
                'tests/sandbox/testAtomSelection.cpp' -Added checks for number of atoms selected, test will bail out if the number
                 of atoms is not what is expected.
                'tests/sandbox/testCoiledCoils.cpp' -Re-purposed this test to generate a coiled-coil with known parameters , then
                 try to fit it
                'myProgs/dwkulp/buildRotamers.cpp', 'myProgs/dwkulp/buildRotamers.h' -simple program to build a specified number
                 of rotamers at a given position on an input structure...no energy calcs
                'myProgs/dwkulp/compareRosettaModels.cpp', 'myProgs/dwkulp/compareRosettaModels.h' -compare position energies from
                 two scored rosetta PDB files, print out PyMOL scripts for analysis
                'myProgs/dwkulp/metalRotamers.cpp', 'myProgs/dwkulp/metalRotamers.h' -build metal sites with proper geometry
                'myProgs/dwkulp/patchAnalysis.cpp', 'myProgs/dwkulp/patchAnalysis.h' -get SASA between two patches on a list of
                 PDBs
                'myProgs/dwkulp/refinePotentialFusions.cpp', 'myProgs/dwkulp/refinePotentialFusions.h' -sampling two trimeric proteins
                 for potentially fusing them together
                'myProgs/dwkulp/setupRosettaFixbb.cpp', 'myProgs/dwkulp/setupRosettaFixbb.h' -easy way to setup a rosetta resfile
                 using known mutations 'A425P' + neighbor flexibility
                'src/Clustering.h', 'src/Clustering.cpp' -Initial implementation of clustering algorithms: SingleLinkage, Average-Linkage,
                 Complete-Linkage and PAM/Kmedoids
                'src/CoiledCoilFitter.h', 'src/CoiledCoilFitter.cpp' -Fitting a coiled-coil to get the parameters
                'src/HelixFit.h', 'src/HelixFit.cpp' -Fitting helical bundles
                'src/RosettaScoredPDBReader.h', 'src/RosettaScoredPDBReader.cpp' -Parsing a Rosetta Scored PDB file
                'tests/sandbox/testClustering.cpp' -A test for clustering code in MSL
1.0.1.14    December 11, 2012    sabs
                'programs/energyOptimizations.h', 'programs/designSideChains.cpp' -added baseline and hydrogen bond builders. Fixed
                 bugs with flexibleNeighbor addition.
1.0.1.13    November 24, 2012    asenes
                'src/AtomContainer.h', 'src/AtomContainer.cpp' -Added support for multiple atom conformations. When reading a PDB
                 or adding atoms, atoms that already exist will be added as alternative conformations, just like for the System.
                 For backward compatibility THIS IS NOT THE DEFAULT. One needs to set the setAddAtomsAsAltCoors(true) function.
                 By default all atoms are still added as separate atoms.
                'tests/sandbox/testAtomContainer.cpp' -Changed the sandbox test to test the support for alt conformations of AtomContainer
                
1.0.1.12    November 23, 2012    asenes
                'src/Atom.h' -Now it returns a bool if setActiveConformation is set to an alternative coordinate that does not
                 exist
                'src/AtomContainer.h' -Added support for switching alt conf, but we need to add support for reading a multi model
                 PDB as alt-conf, right now all atoms are added multiple times
                'var/header.txt' -Updated reference of paper
                'programs/calculateDistanceOrAngle.cpp' -Added support for calculating distance and angles of specific models.
                 To do so I switched intenally from AtomContainer to System.
1.0.1.11    October 16, 2012    sabs
                'src/Scwrl4HBondInteraction.cpp' -Fixed function to select stronger interaction among e1 and e2 - previously was
                 transitioning to e2 only if e1 energy was less than zero.
1.0.1.10    September 28, 2012    bhannigan
                'scripts/submit.py' -Fixing script to copy over branches into tags too.
1.0.1.9    September 27, 2012    sabs
                'src/SelfPairManager.cpp' -onTheFly computation was not working with LinkedPositions because the correct rotamers
                 were not being set, fixed it. Also Removed erroneous check in subdivideInteractions.
1.0.1.8    September 20, 2012    dwkulp
                'Makefile' -GSLDEFAULT to T in trunk
                'testLevels/submit.level' -testCharmmBuild to GOLD and testEZ to LEAD
1.0.1.7    September 18, 2012    dwkulp
                'examples/example_coiled_coils_and_symmetric_bundles.cpp', 'examples/examples.mk' -Update coiled coil example to
                 new CoiledCoil object API
1.0.1.6    August 29, 2012    sabs
                'src/Scwrl4HBondInteraction.cpp', 'src/Scwrl4HBondInteraction.h' -Check if hbond energy will be zero based on distance
                 alone and return 0 immediately. Some code rearrangements to save time.
                'myProgs/sabs/CAHTM.cpp' -Removed EZPotential - v0.0.8
1.0.1.5    August 29, 2012    bhannigan
                'scripts/submit.py' -Fixing bug which prevented new tags from being copied.
1.0.1.4    August 29, 2012    bhannigan
                'scripts/submit.py' -Testing submit script.
1.0.1.3    August 29, 2012    jedonald
                'programs/fillInSideChains.cpp' -Additional program switching from AtomicPairwiseEnergy to CharmmEnergyCalculator
                
1.0.1.2    August 29, 2012    jedonald
                'src/EnergeticAnalysis.cpp', 'src/EnergeticAnalysis.h' -Change from AtomicPairwiseEnergy to CharmmEnergyCalculator.
                 Add a pymolOutput option to specify whether pymol files will be created.
                'src/GSLMinimizer.h', 'src/CharmmEnergyCalculator.cpp', 'programs/energyOptimizations.h', 'programs/tableEnergies.cpp',
                 'programs/energyTable.cpp', 'programs/mutate.cpp', 'programs/findClashes.cpp', 'Makefile' -Change from PairwiseEnergyCalculator
                 to OnTheFlyManager and AtomicPairwiseEnergy to CharmmEnergyCalculator.
                'src/Quench.h', 'src/Quench.cpp' -Change from AtomicPairwiseEnergy to CharmmEnergyCalculator. Add ability to set
                 rotamer levels.
                'programs/runQuench.h', 'programs/runQuench.cpp' -Update program code to have it give more flexible options and
                 output, adding option for rotamer levels.
                'programs/analEnergy.h', 'programs/analEnergy.cpp' -Update program code to make it more flexible options and output.
                 Change from AtomicPairwiseEnergy to CharmmEnergyCalculator. Make pymol output optional.
                'programs/repackSideChains.h' -Use static SysEnv SYSENV found in several other programs.
                'src/AtomicPairwiseEnergy.h', 'src/AtomicPairwiseEnergy.cpp', 'src/PairwiseEnergyCalculator.h', 'src/PairwiseEnergyCalculator.cpp'
                 -Replace AtomicPairwiseEnergy and PairwiseEnergyCalculator with the newer, improved versions (CharmmEnergyCalculator
                 and OnTheFlyManager) as discussed a few months ago.
                'src/BBQTable.cpp', 'src/EnergySet.cpp', 'src/MslTools.cpp', 'src/SelfPairManager.cpp', 'src/TwoBodyDistanceDependentPotentialTable.cpp',
                 'src/PSSMCreator.cpp' -Change ints to doubles and vice versa to remove compile warnings
                'src/FastaReader.cpp' -Move variable inside BOOST ifdef to remove compilation warning
1.0.1.1    July 11, 2012    bhannigan
                'scripts/submit.py' -The submit script only looks at the very top directory when deciding what tree its submitting
                 to. That works fine for trunk, as this is in mslib/trunk. However, branches live in say branches/v.1.0, while
                 currently submit would just see v.1.0. So quick fix that says if not trunk, then assume you're in branches/<whatever>.
                
1.0.1.0    July 11, 2012    bhannigan
                'scripts/submit.py' -Changing submit script so that it will look at what msl directory you are currently using,
                 not assuming trunk anymore.
1.0.0.0    July 08, 2012    asenes
                'README.txt' -Revised the README.txt file to mark this release as v.1.0
                'toppar/charmm22.top', 'toppar/charmm22.par' -Added a link to the official charmm parameter repository
0.23.1.20    July 07, 2012    gevorg
                'src/CrystalLattice.cpp', 'src/CrystalLattice.h' -A little performance enhnacement for CrystalLattice class
                'Makefile' -Added -fopenmp to CCDEBUG flags
0.23.1.19    July 05, 2012    gevorg
                'src/CrystalLattice.cpp', 'src/CrystalLattice.h' -Changed src/CrystalLattice.cpp and scr/CrystalLattice.h - generateCrystal
                 now takes optional parameters to generate only units in contact with the original one. Default behavior is unchanged.
                
                'src/PDBWriter.cpp' -A slight change in PDBWriter.cpp to make sure that when EITHER chain names OR segment names
                 differ, a TER is produced. Also, make sure to compare only the first PDBFormat::L_CHAIN_ID letters of the chain
                 name (often it is useful to use longer chain names to make sure atoms get sorted out correctly, but obviously
                 only the first character of the name is written).
0.23.1.18    July 05, 2012    jedonald
                'scripts/submit.py' -Fix MSL_DIR variable to have trunk/
0.23.1.17    July 05, 2012    jedonald
                'myProgs/myProgs.mk.RENAME_ME' -Add more instructions on paths in myProgs.mk, USERNAME.mk
                'programs/repackSideChains.h', 'programs/repackSideChains.cpp' -Allow SysEnv defaults for rotlib, topfile, and
                 parfiles. For rotamer library levels that do not include ALA, GLY, or PRO, set numRots to zero for these aa's
                
                'src/File.cpp' -Implement doesFileExist function
                'myProgs/jedonald/jedonald.mk', 'myProgs/jedonald/moveProteinCenter.cpp', 'myProgs/jedonald/moveProteinCenter.h'
                 -Add a quick program to translate a pdb
0.23.1.16    June 25, 2012    sabs
                'myProgs/sabs/CAHTM.cpp' -accomodating changes to MonteCarloManager
0.23.1.15    June 22, 2012    sabs
                'src/RotamerLibrary.h', 'src/RotamerLibrary.cpp' -Added methods to removeRotamers, trim the library to a specified
                 level.
                'Makefile', 'programs/trimConformerLibrary.cpp' -The distributed conformer libraries are large and take considerable
                 time to load into memory. This program may be used to trim an MSL conformer library file to a specified level.
                
0.23.1.14    June 19, 2012    sabs
                'programs/getDihedrals.h', 'programs/getDihedrals.cpp' -Added functionality to read a bbdep rotamer library file
                 and assign prob to each sidechain in the pdb.
0.23.1.13    June 19, 2012    dwkulp
                'myProgs/dwkulp/getPairEnergy.cpp', 'myProgs/dwkulp/getPairEnergy.h', 'myProgs/dwkulp/dwkulp.mk' -getPairEnergy
                 program
                'myProgs/dwkulp/compareStructures.cpp' -API change update for findProteinInterfacePositions
0.23.1.12    June 15, 2012    sabs
                'src/ALNReader.cpp', 'src/ALNReader.h', 'src/AtomAngleRelationship.cpp', 'src/AtomAngleRelationship.h', 'src/AtomBondBuilder.h',
                 'src/AtomContainer.cpp', 'src/AtomContainer.h', 'src/Atom.cpp', 'src/AtomDihedralRelationship.cpp', 'src/AtomDihedralRelationship.h',
                 'src/AtomDistanceRelationship.cpp', 'src/AtomDistanceRelationship.h', 'src/AtomGeometricRelationship.cpp', 'src/AtomGeometricRelationship.h',
                 'src/AtomGroup.cpp', 'src/AtomGroup.h', 'src/Atom.h', 'src/AtomicPairwiseEnergy.cpp', 'src/AtomicPairwiseEnergy.h',
                 'src/AtomPointerVector.cpp', 'src/AtomPointerVector.h', 'src/AtomSelection.cpp', 'src/AtomSelection.h', 'src/BackRub.cpp',
                 'src/BackRub.h', 'src/BaselineEnergyBuilder.cpp', 'src/BaselineEnergyBuilder.h', 'src/BaselineInteraction.cpp',
                 'src/BaselineInteraction.h', 'src/BBQTable.cpp', 'src/BBQTable.h', 'src/BBQTableReader.cpp', 'src/BBQTableReader.h',
                 'src/BBQTableWriter.cpp', 'src/BBQTableWriter.h', 'src/CartesianGeometry.cpp', 'src/CartesianGeometry.h', 'src/CartesianPoint.cpp',
                 'src/CartesianPoint.h', 'src/CCD.cpp', 'src/CCD.h', 'src/Chain.cpp', 'src/Chain.h', 'src/CharmmAngleInteraction.cpp',
                 'src/CharmmAngleInteraction.h', 'src/CharmmBondInteraction.cpp', 'src/CharmmBondInteraction.h', 'src/CharmmDihedralInteraction.cpp',
                 'src/CharmmDihedralInteraction.h', 'src/CharmmEEF1Interaction.cpp', 'src/CharmmEEF1Interaction.h', 'src/CharmmEEF1ParameterReader.cpp',
                 'src/CharmmEEF1ParameterReader.h', 'src/CharmmEEF1RefInteraction.cpp', 'src/CharmmEEF1RefInteraction.h', 'src/CharmmElectrostaticInteraction.cpp',
                 'src/CharmmElectrostaticInteraction.h', 'src/CharmmEnergyCalculator.cpp', 'src/CharmmEnergyCalculator.h', 'src/CharmmEnergy.cpp',
                 'src/CharmmEnergy.h', 'src/CharmmImproperInteraction.cpp', 'src/CharmmImproperInteraction.h', 'src/CharmmParameterReader.cpp',
                 'src/CharmmParameterReader.h', 'src/CharmmSystemBuilder.cpp', 'src/CharmmSystemBuilder.h', 'src/CharmmTopologyReader.cpp',
                 'src/CharmmTopologyReader.h', 'src/CharmmTopologyResidue.cpp', 'src/CharmmTopologyResidue.h', 'src/CharmmUreyBradleyInteraction.cpp',
                 'src/CharmmUreyBradleyInteraction.h', 'src/CharmmVdwInteraction.cpp', 'src/CharmmVdwInteraction.h', 'src/ChiStatistics.cpp',
                 'src/ChiStatistics.h', 'src/CoiledCoils.cpp', 'src/CoiledCoils.h', 'src/ConformationEditor.cpp', 'src/ConformationEditor.h',
                 'src/CoordAxes.h', 'src/CRDFormat.cpp', 'src/CRDFormat.h', 'src/CRDReader.cpp', 'src/CRDReader.h', 'src/CRDWriter.cpp',
                 'src/CRDWriter.h', 'src/CrystalLattice.cpp', 'src/CrystalLattice.h', 'src/DeadEndElimination.cpp', 'src/DeadEndElimination.h',
                 'src/DegreeOfFreedomReader.cpp', 'src/DegreeOfFreedomReader.h', 'src/EnergeticAnalysis.cpp', 'src/EnergeticAnalysis.h',
                 'src/EnergySet.cpp', 'src/EnergySet.h', 'src/Enumerator.cpp', 'src/Enumerator.h', 'src/EnvironmentDatabase.cpp',
                 'src/EnvironmentDatabase.h', 'src/EnvironmentDescriptor.cpp', 'src/EnvironmentDescriptor.h', 'src/EZpotentialBuilder.cpp',
                 'src/EZpotentialBuilder.h', 'src/EZpotentialInteraction.cpp', 'src/EZpotentialInteraction.h', 'src/FastaReader.cpp',
                 'src/FastaReader.h', 'src/File.cpp', 'src/File.h', 'src/FormatConverter.cpp', 'src/FormatConverter.h', 'src/FourBodyInteraction.cpp',
                 'src/FourBodyInteraction.h', 'src/Frame.cpp', 'src/Frame.h', 'src/FuseChains.cpp', 'src/FuseChains.h', 'src/GSLMinimizer.cpp',
                 'src/GSLMinimizer.h', 'src/Hash.h', 'src/Helanal.cpp', 'src/Helanal.h', 'src/HelixFusion.cpp', 'src/HelixFusion.h',
                 'src/HelixGenerator.cpp', 'src/HelixGenerator.h', 'src/HydrogenBondBuilder.cpp', 'src/HydrogenBondBuilder.h',
                 'src/IcEntry.cpp', 'src/IcEntry.h', 'src/IcTable.cpp', 'src/IcTable.h', 'src/Interaction.cpp', 'src/Interaction.h',
                 'src/InterfaceResidueDescriptor.cpp', 'src/InterfaceResidueDescriptor.h', 'src/LinearProgrammingOptimization.cpp',
                 'src/LinearProgrammingOptimization.h', 'src/Line.cpp', 'src/Line.h', 'src/LogicalCondition.cpp', 'src/LogicalCondition.h',
                 'src/LogicalParser.cpp', 'src/LogicalParser.h', 'src/Matrix.cpp', 'src/Matrix.h', 'src/MIDReader.cpp', 'src/MIDReader.h',
                 'src/Minimizer.cpp', 'src/Minimizer.h', 'src/MoleculeInterfaceDatabase.cpp', 'src/MoleculeInterfaceDatabase.h',
                 'src/MonteCarloManager.cpp', 'src/MonteCarloManager.h', 'src/MonteCarloOptimization.cpp', 'src/MonteCarloOptimization.h',
                 'src/MslExceptions.h', 'src/MslOut.cpp', 'src/MslOut.h', 'src/MslTools.cpp', 'src/MslTools.h', 'src/OneBodyInteraction.cpp',
                 'src/OneBodyInteraction.h', 'src/OnTheFlyManager.cpp', 'src/OnTheFlyManager.h', 'src/OptimalRMSDCalculator.cpp',
                 'src/OptimalRMSDCalculator.h', 'src/OptionParser.cpp', 'src/OptionParser.h', 'src/PairwiseEnergyCalculator.cpp',
                 'src/PairwiseEnergyCalculator.h', 'src/PDBFormat.cpp', 'src/PDBFormat.h', 'src/PDBFragments.cpp', 'src/PDBFragments.h',
                 'src/PDBReader.cpp', 'src/PDBReader.h', 'src/PDBTopologyBuilder.cpp', 'src/PDBTopologyBuilder.h', 'src/PDBTopology.cpp',
                 'src/PDBTopology.h', 'src/PDBWriter.cpp', 'src/PDBWriter.h', 'src/PhiPsiReader.cpp', 'src/PhiPsiReader.h', 'src/PhiPsiStatistics.cpp',
                 'src/PhiPsiStatistics.h', 'src/PhiPsiWriter.cpp', 'src/PhiPsiWriter.h', 'src/PolymerSequence.cpp', 'src/PolymerSequence.h',
                 'src/Position.cpp', 'src/Position.h', 'src/PotentialTable.cpp', 'src/PotentialTable.h', 'src/Predicate.cpp', 'src/Predicate.h',
                 'src/PrincipleComponentAnalysis.cpp', 'src/PrincipleComponentAnalysis.h', 'src/PrositeReader.cpp', 'src/PrositeReader.h',
                 'src/PSFReader.cpp', 'src/PSFReader.h', 'src/PSSMCreator.cpp', 'src/PSSMCreator.h', 'src/PyMolVisualization.cpp',
                 'src/PyMolVisualization.h', 'src/PythonMSL.cpp', 'src/Quaternion.cpp', 'src/Quaternion.h', 'src/Quench.cpp', 'src/Quench.h',
                 'src/RandomNumberGenerator.cpp', 'src/RandomNumberGenerator.h', 'src/RandomSeqGenerator.cpp', 'src/RandomSeqGenerator.h',
                 'src/Reader.cpp', 'src/Reader.h', 'src/Real.h', 'src/RegEx.cpp', 'src/RegEx.h', 'src/release.h', 'src/Residue.cpp',
                 'src/Residue.h', 'src/ResiduePairTable.cpp', 'src/ResiduePairTable.h', 'src/ResiduePairTableReader.cpp', 'src/ResiduePairTableReader.h',
                 'src/ResidueSelection.cpp', 'src/ResidueSelection.h', 'src/ResidueSubstitutionTable.cpp', 'src/ResidueSubstitutionTable.h',
                 'src/ResidueSubstitutionTableReader.cpp', 'src/ResidueSubstitutionTableReader.h', 'src/RotamerLibrary.cpp', 'src/RotamerLibrary.h',
                 'src/RotamerLibraryReader.cpp', 'src/RotamerLibraryReader.h', 'src/SasaCalculator.cpp', 'src/SasaCalculator.h',
                 'src/Scwrl4HBondInteraction.cpp', 'src/Scwrl4HBondInteraction.h', 'src/Selectable.h', 'src/SelfConsistentMeanField.cpp',
                 'src/SelfConsistentMeanField.h', 'src/SelfPairManager.cpp', 'src/SelfPairManager.h', 'src/SidechainOptimizationManager.cpp',
                 'src/SidechainOptimizationManager.h', 'src/SphericalPoint.cpp', 'src/SphericalPoint.h', 'src/SpringConstraintInteraction.cpp',
                 'src/SpringConstraintInteraction.h', 'src/Symmetry.cpp', 'src/Symmetry.h', 'src/SysEnv.cpp', 'src/SysEnv.h', 'src/System.cpp',
                 'src/System.h', 'src/SystemRotamerLoader.cpp', 'src/SystemRotamerLoader.h', 'src/TBDReader.cpp', 'src/TBDReader.h',
                 'src/ThreeBodyInteraction.cpp', 'src/ThreeBodyInteraction.h', 'src/Timer.cpp', 'src/Timer.h', 'src/Transforms.cpp',
                 'src/Transforms.h', 'src/Tree.cpp', 'src/Tree.h', 'src/triple.h', 'src/TwoBodyDistanceDependentPotentialTable.cpp',
                 'src/TwoBodyDistanceDependentPotentialTable.h', 'src/TwoBodyInteraction.cpp', 'src/TwoBodyInteraction.h', 'src/UserDefinedEnergy.cpp',
                 'src/UserDefinedEnergy.h', 'src/UserDefinedEnergySetBuilder.cpp', 'src/UserDefinedEnergySetBuilder.h', 'src/UserDefinedInteraction.cpp',
                 'src/UserDefinedInteraction.h', 'src/VectorHashing.cpp', 'src/VectorHashing.h', 'src/VectorPair.cpp', 'src/VectorPair.h',
                 'src/Writer.cpp', 'src/Writer.h' -Updating header with final paper citation
                'examples/example_add_atoms_to_System_and_AtomContainer.cpp', 'examples/example_add_identity_to_position.cpp',
                 'examples/example_AtomPointerVector.cpp', 'examples/example_backrub.cpp', 'examples/example_bbq.cpp', 'examples/example_ccd.cpp',
                 'examples/example_coiled_coils_and_symmetric_bundles.cpp', 'examples/example_looping_over_Chain_Residues_Atoms.cpp',
                 'examples/example_measurements.cpp', 'examples/example_molecular_alignment.cpp', 'examples/example_multipleAtomsCoordinates.cpp',
                 'examples/example_multiple_coordinates_from_NMR_multiModel_PDB.cpp', 'examples/example_multipleResidueIdentities.cpp',
                 'examples/example_mutation_rotamers.cpp', 'examples/example_pdbfrag.cpp', 'examples/example_read_write_PDBs_with_the_AtomContainer.cpp',
                 'examples/example_read_write_PDBs_with_the_System.cpp', 'examples/example_regular_expressions.cpp', 'examples/example_SasaCalculator_usage.cpp',
                 'examples/example_selecting_atoms.cpp' -Updating header with final paper citation
                'README.txt', 'programs/alignMolecules.cpp', 'programs/analEnergy.cpp', 'programs/analEnergy.h', 'programs/backrubPdb.cpp',
                 'programs/backrubPdb.h', 'programs/calculateDistanceOrAngle.cpp', 'programs/calculateSasa.cpp', 'programs/coiledCoilBuilder.cpp',
                 'programs/createEBL.cpp', 'programs/createEnergyTable.cpp', 'programs/createFragmentDatabase.cpp', 'programs/createFragmentDatabase.h',
                 'programs/designSideChains.cpp', 'programs/energyOptimizations.h', 'programs/energyTable.cpp', 'programs/fillInSideChains.cpp',
                 'programs/fillInSideChains.h', 'programs/findClashes.cpp', 'programs/findClashes.h', 'programs/generateCoiledCoils.cpp',
                 'programs/generateCoiledCoils.h', 'programs/generateCrystalLattice.cpp', 'programs/generateCrystalLattice.h',
                 'programs/getChiRecovery.cpp', 'programs/getChiRecovery.h', 'programs/getDihedrals.cpp', 'programs/getDihedrals.h',
                 'programs/getSelection.cpp', 'programs/getSelection.h', 'programs/getSphericalCoordinates.cpp', 'programs/getSphericalCoordinates.h',
                 'programs/getSurroundingResidues.cpp', 'programs/getSurroundingResidues.h', 'programs/grepSequence.cpp', 'programs/grepSequence.h',
                 'programs/insertLoopIntoTemplate.cpp', 'programs/insertLoopIntoTemplate.h', 'programs/minimize.cpp', 'programs/minimize.h',
                 'programs/mutate.cpp', 'programs/mutate.h', 'programs/optimizeLP.cpp', 'programs/optimizeMC.cpp', 'programs/printSequence.cpp',
                 'programs/printSequence.h', 'programs/renumberResidues.cpp', 'programs/renumberResidues.h', 'programs/repackSideChains.cpp',
                 'programs/repackSideChains.h', 'programs/runKBQuench.cpp', 'programs/runKBQuench.h', 'programs/runQuench.cpp',
                 'programs/runQuench.h', 'programs/searchFragmentDatabase.cpp', 'programs/searchFragmentDatabase.h', 'programs/setConformation.cpp',
                 'programs/tableEnergies.cpp', 'programs/tableEnergies.h' -Updating header with final paper citation
                'tests/failSafeTests/utility_createAminoAcids.cpp', 'tests/gold/testCharmmBuild.cpp', 'tests/gold/testCharmmEnergies.cpp',
                 'tests/gold/testEZpotential.cpp', 'tests/gold/testRMSDalignment.cpp' -Updating header with final paper citation
                
                'tests/sandbox/testAddCharmmIdentity.cpp', 'tests/sandbox/testALNReader.cpp', 'tests/sandbox/testAtomAndResidueId.cpp',
                 'tests/sandbox/testAtomBondBuilder.cpp', 'tests/sandbox/testAtomContainer.cpp', 'tests/sandbox/testAtomGroup.cpp',
                 'tests/sandbox/testAtomPointerVector.cpp', 'tests/sandbox/testAtomSelection.cpp', 'tests/sandbox/testBackRub.cpp',
                 'tests/sandbox/testBBQ2.cpp', 'tests/sandbox/testBBQ.cpp', 'tests/sandbox/testBoost.cpp', 'tests/sandbox/testCCD.cpp',
                 'tests/sandbox/testCharmmEEF1ParameterReader.cpp', 'tests/sandbox/testCharmmEnergies.cpp', 'tests/sandbox/testCharmmTopologyReader.cpp',
                 'tests/sandbox/testCoiledCoils.cpp', 'tests/sandbox/testConformationEditor.cpp', 'tests/sandbox/testCRDIO.cpp',
                 'tests/sandbox/testData.h', 'tests/sandbox/testDeleteBondedAtom.cpp', 'tests/sandbox/testDerivatives.cpp', 'tests/sandbox/testEEF1_2.cpp',
                 'tests/sandbox/testEEF1.cpp', 'tests/sandbox/testEnergeticAnalysis.cpp', 'tests/sandbox/testEnergySet.cpp', 'tests/sandbox/testEnvironmentDatabase.cpp',
                 'tests/sandbox/testEnvironmentDescriptor.cpp', 'tests/sandbox/testFormatConverter.cpp', 'tests/sandbox/testFrame.cpp',
                 'tests/sandbox/testGenerateCrystalLattice.cpp', 'tests/sandbox/testHelixFusion.cpp', 'tests/sandbox/testHelixGenerator.cpp',
                 'tests/sandbox/testIcBuilding.cpp', 'tests/sandbox/testLinkedPositions.cpp', 'tests/sandbox/testLoopOverResidues.cpp',
                 'tests/sandbox/testMinimization.cpp', 'tests/sandbox/testMolecularInterfaceDatabase.cpp', 'tests/sandbox/testMslOut2.cpp',
                 'tests/sandbox/testMslOut.cpp', 'tests/sandbox/testMslToolsFunctions.cpp', 'tests/sandbox/testNonBondedCutoff.cpp',
                 'tests/sandbox/testOptimalRMSDCalculator.cpp', 'tests/sandbox/testPDBFragments.cpp', 'tests/sandbox/testPDBIO.cpp',
                 'tests/sandbox/testPDBTopologyBuild.cpp', 'tests/sandbox/testPDBTopology.cpp', 'tests/sandbox/testPhiPsi.cpp',
                 'tests/sandbox/testPolymerSequence.cpp', 'tests/sandbox/testPSFReader.cpp', 'tests/sandbox/testPSFReader.h', 'tests/sandbox/testQuench.cpp',
                 'tests/sandbox/testQuench.h', 'tests/sandbox/testRandomNumberGenerator.cpp', 'tests/sandbox/testRandomSeqGenerator.cpp',
                 'tests/sandbox/testRegEx.cpp', 'tests/sandbox/testResiduePairTable.cpp', 'tests/sandbox/testResidueSelection.cpp',
                 'tests/sandbox/testResidueSubstitutionTable.cpp', 'tests/sandbox/testRInterface.cpp', 'tests/sandbox/testRotamerLibraryWriter.cpp',
                 'tests/sandbox/testRotamerOptimization.cpp', 'tests/sandbox/testSasaCalculator.cpp', 'tests/sandbox/testSaveAtomAltCoor.cpp',
                 'tests/sandbox/testSharedPointers2.cpp', 'tests/sandbox/testSurfaceAreaAndVolume.cpp', 'tests/sandbox/testSymmetry.cpp',
                 'tests/sandbox/testSysEnv.cpp', 'tests/sandbox/testSystemCopy.cpp', 'tests/sandbox/testSystemIcBuilding.cpp',
                 'tests/sandbox/testTokenize.cpp', 'tests/sandbox/testTransformBondAngleDiheEdits.cpp', 'tests/sandbox/testTransforms.cpp',
                 'tests/sandbox/testTree.cpp', 'tests/sandbox/testVectorPair.cpp' -Updating header with final paper citation
0.23.1.11    June 14, 2012    sabs
                'src/IcEntry.h', 'src/Atom.cpp', 'src/IcEntry.cpp', 'src/Atom.h' -Changed logic for building from internal coordiantes.
                 IcEntries that were already tried during build and failed should not be tried again. The old code did this based
                 on the first atom in the icentry:changed to base it on the icentry pointer.
                'src/ConformationEditor.cpp' -Commented out some debug couts.
                'programs/repackSideChains.h', 'programs/repackSideChains.cpp' -Added an option to specify which positions should
                 not be repacked.
0.23.1.10    June 11, 2012    bhannigan
                'src/BBQTable.h', 'src/CartesianPoint.h' -Removing lt and gt operators for CartesianPoint as they don't make a
                 lot of sense. Instead, creating a custom CartesianPointCompare class so it can be used with std::map
0.23.1.9    June 07, 2012    sabs
                'library/EBRL_11-2011_CHARMM22.txt', 'library/EBRL_11-2011_PDB2.3.txt' -Changed level names in EBRL to SLxx.xx
                 format to conform with EBL.
0.23.1.8    June 05, 2012    gevorg
                'tests/sandbox/testOptimalRMSDCalculator.cpp' -A somewhat more fair comparison between Transform::rmsdAlignment
                 and OptimalRMSDCalculator.
0.23.1.7    June 05, 2012    gevorg
                'src/OptimalRMSDCalculator.h', 'src/OptimalRMSDCalculator.cpp', 'tests/sandbox/testOptimalRMSDCalculator.cpp' -Updated
                 API for OptimalRMSDCalculator a bit.
0.23.1.6    June 05, 2012    gevorg
                'Makefile' -Forgot to include an updated Makefile when previously committed OptimalRMSDAligner.
0.23.1.5    June 04, 2012    gevorg
                'src/OptimalRMSDCalculator.h', 'src/OptimalRMSDCalculator.cpp' -Made the third parameter in src/OptimalRMSDCalculator::align()
                 optional.
0.23.1.4    June 04, 2012    gevorg
                'src/OptimalRMSDCalculator.h', 'src/OptimalRMSDCalculator.cpp' -Added a class for very rapid calculation of RMSD
                 upon optimal superposition.
                'tests/sandbox/testOptimalRMSDCalculator.cpp' -Test for OptimalRMSDCalculator class that demonstrates its efficiency
                 relative to the current method of calculating RMSD.
0.23.1.3    June 04, 2012    bhannigan
                'scripts/submit.py', 'scripts/mslBuildTools.py', 'testLevels/submit.level' -Changing submit script to look for
                 GOLD or LEAD.
0.23.1.2    June 01, 2012    asenes
                'programs/energyOptimizations.h' -Changed the obsolete reference to the environmental variable MSL_EBL to MSL_ROTLIB
                
0.23.1.1    June 01, 2012    sabs
                'programs/repackSideChains.cpp' -Removed unnecessary call to setVariablePositions
                'myProgs/sabs/discardSimilarConformers.cpp', 'myProgs/sabs/selectCavities.cpp', 'myProgs/sabs/buildDunbrackLibrary.cpp',
                 'myProgs/sabs/makeShettyLibrary.cpp', 'myProgs/sabs/jobDistributor.cpp', 'myProgs/sabs/makePdbFromRotamerLibrary.cpp',
                 'myProgs/sabs/getInteractionGraph.cpp', 'myProgs/sabs/hbondRecovery.cpp', 'myProgs/sabs/makeHonigLibrary.cpp',
                 'myProgs/sabs/genHeteroUniverse.cpp', 'myProgs/sabs/genHomoUniverse.cpp', 'myProgs/sabs/CAHTM.cpp', 'myProgs/sabs/clusterHeteroCandidates.cpp',
                 'myProgs/sabs/clusterCandidates.cpp', 'myProgs/sabs/tmHeteroRulesCreator.cpp', 'myProgs/sabs/tmRulesCreator.cpp',
                 'myProgs/sabs/chainizeMissingLoops.cpp', 'myProgs/sabs/analyzeTMStructures.cpp', 'myProgs/sabs/convertToPdbNames.cpp',
                 'myProgs/sabs/getUniqueChains.cpp', 'myProgs/sabs/sabs.mk' -Moved some programs to myProgs/sabs directory
0.23.1.0    May 31, 2012    asenes
                'toppar/charmm22.par', 'toppar/charmm22.top', 'toppar/charmm22.origIC.top' -Revised CHARMM topology and parameter
                 files (derived from top_all22_prot.inp). The topology file has the IC set in the conformation of the top conformer
                 of the Energy-Based library
                'tests/gold/testCharmmEnergies.cpp' -Revised to work with the new toppar files charmm22.top charmm22.par (since
                 the IC changes the conformation of the molecule built in the test also changes) and the EBL library. Added also
                 an option for writing PDB files
                'tests/sandbox/testCharmmBuild.cpp' -Revised and moved to the gold directory
                'tests/gold/testCharmmBuild.cpp' -New gold test for building from polymer sequence and from PDB file. This also
                 assumes to be used with toppar/charmm22.par toppar/charmm22.top
                'testLevels/submit.level' -Removed some old tests and added all the gold files except the one that requires GSL
                 (testRMSDalignment)
                'Makefile' -Added testCharmmBuild.cpp to the gold files and removed from the sandbox
                'src/SysEnv.h', 'src/SysEnv.cpp' -Changed the default toppar files to /toppar/charmm22.top toppar/charmm22.par
                 and the default rotamer library to library/EBL_11-2011_CHARMM22.txt
0.23.0.4    May 25, 2012    asenes
                'tests/sandbox/testRMSDalignment.cpp', 'tests/lead/testCharmmEnergies.cpp' -Modified the tests and moved them to
                 the test/gold directory
                'tests/gold/testRMSDalignment.cpp' -Modified test and added to the test/gold directory, checks for coordinates
                 and RMSD after alingments
                'tests/gold/testCharmmEnergies.cpp' -Modified test and added to the test/gold directory, run 4 test using System
                 energies, SelfPairManager and OnTheFlyManager and compares with expected energies
                'src/OnTheFlyManager.cpp', 'src/OnTheFlyManager.h' -Added string getSummary and reformatted the output of printSummary
                 to align the number and use 6 decimal digits
                'Makefile' -Added new gold test testRMSDalignment and removed the same from the tests/sandbox
0.23.0.3    May 24, 2012    asenes
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added functions setBuildTerm(term, bool), setBuildAllTerms
                 and setBuildNoTerms to selectively decide what energy terms are built in case some won't be used.
                'Makefile' -Fixed bug, it was not removing tests/lead binary on make clean
                'src/System.h', 'src/EnergySet.h', 'src/EnergySet.cpp', 'src/SelfPairManager.h', 'src/SelfPairManager.cpp', 'tests/lead/testCharmmEnergies.cpp'
                 -Added optional argument to specify precision in printSummary and getSummary and getEnergySummary functions (the
                 defauil is back to 6, it was changed to 15 in the previous commit)
0.23.0.2    May 24, 2012    dwkulp
                'programs/optimizeMC.cpp' -Enum usage for INIT and MCSHAPE
0.23.0.1    May 24, 2012    dwkulp
                'programs/repackSideChains.h' -MC shape arguement now uses enum and not string compare
0.23.0.0    May 24, 2012    dwkulp
                'Makefile' -PhiPsiWriter OnTheFlyManager CharmmEnergyCalculator designSideChains generateCoiledCoils
                'myProgs/dwkulp/compareStructures.cpp' -Added findProteinInterfacePositions that was removed from System
                'myProgs/dwkulp/designLinearSegment.cpp', 'myProgs/dwkulp/designLinearSegment.h' -Added Linear,Stem and Spot type
                 searches for PDBFragments
                'src/PSSMCreator.h', 'src/PSSMCreator.cpp', 'myProgs/dwkulp/querySeqCons.cpp', 'myProgs/dwkulp/querySeqCons.h'
                 -entropy vs freq vs mostfreq sequence conservation added
                'myProgs/dwkulp/dwkulp.mk', 'myProgs/dwkulp/addSequenceConservation.cpp', 'myProgs/dwkulp/addSequenceConservation.h',
                 'myProgs/dwkulp/calcSasaInterface.h', 'myProgs/dwkulp/calcSasaInterface.cpp', 'myProgs/dwkulp/checkCloseTerminii.h',
                 'myProgs/dwkulp/checkCloseTerminii.cpp', 'myProgs/dwkulp/circularPermutant.cpp', 'myProgs/dwkulp/circularPermutant.h',
                 'myProgs/dwkulp/clashCheck.cpp', 'myProgs/dwkulp/clashCheck.h', 'myProgs/dwkulp/designCheck.cpp', 'myProgs/dwkulp/designCheck.h',
                 'myProgs/dwkulp/domainSasa.cpp', 'myProgs/dwkulp/domainSasa.h', 'myProgs/dwkulp/getVectorPairs.cpp', 'myProgs/dwkulp/getVectorPairs.h',
                 'myProgs/dwkulp/helixWindowAnalysis.h', 'myProgs/dwkulp/helixWindowAnalysis.cpp', 'myProgs/dwkulp/setupRemodelFusions.cpp',
                 'myProgs/dwkulp/setupRemodelFusions.h', 'myProgs/dwkulp/setupRosettaMSA.cpp', 'myProgs/dwkulp/setupRosettaMSA.h'
                 -New programs
                'programs/calculateSasa.cpp' -Computes Buried hydrophobic percent
                'programs/designSideChains.cpp' -Equivalent program to RepackSideChains, except allow identity changes
                'programs/generateCoiledCoils.cpp' -updated to use new CoiledCoil API
                'programs/energyOptimizations.h' -updated to use new loadRotamers call, added hbond support, MonteCarlo options,
                 AnalysisOptions added
                'programs/getDihedrals.h', 'programs/getDihedrals.cpp' -Option to write out phi-psi table
                'programs/printSequence.cpp', 'programs/printSequence.h' -Option to print in FASTA format
                'programs/insertLoopIntoTemplate.cpp', 'programs/insertLoopIntoTemplate.h' -Option to use a single chain from a
                 template structure and chain for fragment structure, also it now can output Rosetta input files for a program
                 called fusionDesign (not checked into RosettaSVN yet)..still useful for printing fixbb resfiles
                'programs/renumberResidues.cpp', 'programs/renumberResidues.h' -Option to useIcodes and outPdb
                'src/CharmmParameterReader.cpp' -pre-computing pairs of parameters had a bug, using push_back instead of direct
                 indexing
                'src/EnergySet.h', 'src/EnergySet.cpp', 'src/CharmmSystemBuilder.cpp', 'src/CharmmSystemBuilder.h' -remove createPairwiseTable
                 hack
                'src/FuseChains.cpp', 'src/FuseChains.h' -Keep track of inserted positions
                'src/MonteCarloManager.h' -Added enum ANNEALTYPES
                'src/MonteCarloOptimization.h', 'src/MonteCarloOptimization.cpp' -Added different types of initialization algorithms;
                 get rid of verbose, use MSLOUT.stream(); each move is only 1 position, 1 rotamer each step... this has to be the
                 default.
                'src/MslOut.h', 'src/MslOut.cpp' -Added a static lookup map, for objects to populate with output data. This way
                 you can check for internal consistency
                'src/MslTools.cpp', 'src/MslTools.h' -Added string replace function
                'src/OnTheFlyManager.h', 'src/OnTheFlyManager.cpp', 'src/CharmmEnergyCalculator.h', 'src/CharmmEnergyCalculator.cpp'
                 -New on-the
                'src/PDBFragments.h', 'src/PDBFragments.cpp' -Search in linear,stem or spot mode
                'src/PhiPsiStatistics.h', 'src/PhiPsiStatistics.cpp' -set grid size, getPhiPsiBin is public
                'src/PhiPsiWriter.h', 'src/PhiPsiWriter.cpp' -new object
                'src/PyMolVisualization.cpp', 'src/PyMolVisualization.h' -createSelection function added
                'src/SasaCalculator.h', 'src/SasaCalculator.cpp' -getResidueSasa function was implemented
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Use MonteCarloManager:ANNEALTYPES, added computeSelfE and computeBestPairE
                
                'src/SysEnv.cpp' -Added MSL_PHIPSI_TABLE
                'src/System.h', 'src/System.cpp' -removed the findProteinInterfacePositions
                'src/VectorPair.h', 'src/VectorPair.cpp', 'src/VectorHashing.cpp' -Added more distances and angles to the vector
                 pair
                'src/PrositeReader.h', 'src/PrositeReader.cpp' -Reader for prosite database
                'tests/lead/testCharmmEnergies.cpp' -Test to insure that energies are computed the same using different mechanisms
                
                'src/AtomicPairwiseEnergy.h', 'src/AtomicPairwiseEnergy.cpp' -removed calls to get bonded terms, we have a different
                 mechanism now, this object is old and will retire soon.
0.22.8.1    May 22, 2012    gevorg
                'src/CharmmSystemBuilder.cpp' -CharmmSystemBuilder::updateNonBonded() does not pre-build the atom box information
                 if non-bond cutoffs are specified as 0.
0.22.8.0    May 21, 2012    asenes
                'src/EZpotentialInteraction.h', 'src/EZpotentialInteraction.cpp' -EZ potential single body interaction.
                'src/EZpotentialBuilder.h', 'src/EZpotentialBuilder.cpp' -A builder to add EZ potential interactions to the System
                
                'tests/gold/testEZpotential.cpp' -A test for EZ potential. This is a GOLD/LEAD test (the first one contributed
                 to MSL)
                'src/Atom.h', 'src/Atom.cpp', 'src/AtomGroup.h', 'src/AtomGroup.cpp', 'src/Residue.h', 'src/Residue.cpp', 'src/Position.h',
                 'src/Position.cpp', 'src/Chain.h', 'src/Chain.cpp' -Added functions isPositionNterminal and isPositionCterminal.
                 The function is performed by the Chain: the Atom, Residue, Position access it calling their parent. This check
                 requires a System (because it requires a Chain), it does not work if the atom is in a AtomContainer at this point
                
                'Makefile' -Added EZpotentialBuilder EZpotentialInteraction and testEZpotential. Also added a section for the GOLD
                 tests
0.22.7.9    May 19, 2012    gevorg
                'src/Atom.cpp', 'src/Atom.h' -Atom::getGroupGeometricPoint() now returns CartesianPoint&, rather than CartesianPoint.
                
                'src/AtomPointerVector.cpp', 'src/AtomPointerVector.h' -AtomPointerVector::getGeometricCenter() now returns a reference
                 to CartesianPoint, since the geometric center is stored as a CartesianPoint in the AtomPointerVector object, so
                 why not return a reference and save on a copy constructor call if the user wants that.
                'src/CharmmSystemBuilder.cpp', 'src/CharmmSystemBuilder.h' -Added group-based cutoffs as an option. If this is
                 set as true, CharmmSystemBuilder::updateNonBonded() will use group-based cutoffs (which was the only type of cutoff
                 used before).
0.22.7.8    May 19, 2012    asenes
                'src/EnergySet.h', 'src/EnergySet.cpp', 'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp', 'tests/sandbox/testCharmmEnergies.cpp'
                 -Added function setCreatePairwiseTable(bool _flag) that allows to DISABLE the creation of the pairwise tables
                 in the EnergySet used for On-The-Fly calculations. These tables should be soon moved to a separate object for
                 improved performance.
                'Makefile' -Removed testTree from the Makefile (commented out) since it segfaults
                'testLevels/submit.level' -Removed testTree from the tests that run, added testCharmmEnergies
0.22.7.7    May 19, 2012    asenes
                'tests/sandbox/testLinkedPositions.cpp' -Improved the comments in the test. Added also a threshold for the comparison
                 of the energies, now differences <1E-10 are considered OK
0.22.7.6    May 18, 2012    asenes
                'src/EnergySet.cpp' -Fixed bug in function eraseTerm, it was not deleting the pairInteractions maps
                'Makefile' -Changed the Makefile according to the fact that the old messy tests have been moved to tests/sandbox.
                 Now the tests are referred as SANDBOX in the Makefile
0.22.7.5    May 17, 2012    bhannigan
                'scripts/submit.py', 'scripts/miscUtils.py', 'scripts/mslBuildTools.py', 'testLevels/submit.level', 'tests/testBBQ2.cpp'
                 -Changing submit script to now run tests and report on results. Also, submit now requires double dashes and not
                 single dashes.
0.22.7.4    May 16, 2012    gevorg
                'src/CharmmSystemBuilder.cpp', 'src/CharmmSystemBuilder.h' -Updated two things in CharmmSystemBuilder::updateNonBonded().
                 1) Added an optional boolean argument, which when set as true (defaults to false), will not create interaction
                 objects for atom pairs that are both fixed (i.e. have only one alternative conformation). This saves a lot of
                 time and memory in cases where a large portion of the molecule is not changing during the simulation (e.g. the
                 'template' in design). Because the argument defaults to false, the default behavior does not change. 2) Implemented
                 'smart cutoff' for deciding which atom pairs to create interactions for. Before, the cutoffs were based on atomic
                 positions in the 'current' state, which can give artifacts in design which depend on the starting conformation.
                 Now, the set of resulting interactions is guaranteed to contain ALL atom pairs that are ever within the interaction
                 cutoff distance of each other (in any of their conformations). It will, however, be a super
0.22.7.3    May 15, 2012    dwkulp
                'Makefile', 'examples/examples.mk', 'tests/testBackRub.cpp', 'tests/testPDBFragments.cpp' -Re
0.22.7.2    May 13, 2012    asenes
                'src/Atom.h', 'src/PDBTopologyBuilder.cpp', 'src/SelfPairManager.cpp' -Added (Atom*) declaration for NULL pointer
                 instantiation such as vector<Atom*>(<num>, (Atom*)NULL) and similar for compiler compatibility (some compilers
                 take NULL as an int)
                'programs/optimizeMC.cpp' -Was broken and did not compile, fixed typo (: instead of ; in for loop)
                'Makefile' -Temporarily disabled GLPK, it does not compile in Ubuntu, needs to be looked into. Changed header to
                 better explain use of external libraries
0.22.7.1    May 09, 2012    gevorg
                'src/AtomContainer.cpp' -This complets an earlier commit, where I forgot to include AtomContainer.cpp. Fixed cthe
                 opy constructor for AtomContainner and enabled construction from from stringstream.
0.22.7.0    May 07, 2012    asenes
                'Makefile' -Put objects PDBFragments and HelixFusion objects under GSL requirement. Added test program testRMSDalignment.
                
                'README.txt', 'var/header.txt' -Updated MSL reference (still in press but with DOI).
                'src/System.cpp', 'src/System.h' -Deprecated function findProteinInterfacePositions, needs to be moved to its own
                 object. Now the function cerr a warning.
                'src/HelixFusion.cpp', 'src/HelixFusion.h', 'src/PDBFragments.cpp', 'src/PDBFragments.h' -Now these object requires
                 GSL to compile. Update headers.
                'src/MslOut.h', 'src/MslOut.cpp' -MslOut is now turned off by default.
                'src/Transforms.h', 'src/Transforms.cpp' -Now the RMSD functions (rmsdAlignment and smartRmsdAlignment) require
                 GSL. Now the default option of smartRmsdAlignment is to use atomids (before was atom names) and removed option
                 to use the reference address (not used and possibly bogus). Added a getLastRMSD function, returns the RMSD obtained
                 with the last alignment.
                'src/HelixGenerator.h', 'src/HelixGenerator.cpp' -Moved under GSL requirement the function getHelixAxis. NOTE:
                 this is a global function, not an object function. Should be changed (its own object or MSL tools?). Does not
                 seem to be used.
                'src/Quaternion.h', 'src/Quaternion.cpp' -Put functions makeQuaternion(AtomPointerVector &_align, AtomPointerVector
                 &_ref) and getPrincipalAxes dependent on GSL.
                'tests/testRMSDalignment.cpp' -Added a test for rmsdAlignment and smartRmsdAlignment
0.22.6.2    May 07, 2012    gevorg
                'src/Transforms.cpp' -Commented out two lines in Transforms::rmsdAlignment that used to print to stdout the transformation
                 and rotation matrices upon alignment
0.22.6.1    May 04, 2012    asenes
                'src/PDBReader.cpp' -Fixed bug in read, section parsing the missing residues
0.22.6.0    May 03, 2012    sabs
                'src/SysEnv.cpp' -Changed MSL_EBL to the committed library file.
                'src/FormatConverter.h', 'src/FormatConverter.cpp' -Removed Terminal Flags for charmm to pdb conversion. Created
                 a function to take just atom and res names and provide conversion.Changed the overall object to be independant
                 of System
                'src/PDBWriter.h', 'src/PDBWriter.cpp' -Added functionality to convert and print CRD names.
                'tests/testFormatConverter.cpp' -Updated the test for the new interface.
                'src/CRDWriter.cpp' -Print the number of atoms in the first line.
0.22.5.9    April 26, 2012    sabs
                'src/CRDWriter.cpp' -The absolute residue number was not updated correctly. Fixed it.
0.22.5.8    April 25, 2012    sabs
                'src/CRDFormat.h', 'src/CRDFormat.cpp' -Fixed the format string: there was slight misalignment.
                'src/CRDReader.h', 'src/CRDReader.cpp' -Removed call to deletePointers from setup, initialized member pointer to
                 NULL.
                'src/CRDWriter.h', 'src/CRDWriter.cpp', 'tests/testCRDIO.cpp' -Added writeREMARKS function.
0.22.5.7    April 20, 2012    sabs
                'src/Scwrl4HBondInteraction.cpp' -printParameters function was printing wrong values.
                'programs/createEBL.cpp' -Added the configfile option.
                'toppar/scwrl4hb/par_hbond_2.txt', 'toppar/scwrl4hb/par_hbond_CA_1.txt', 'toppar/scwrl4hb/par_hbond_CA_2.txt' -New
                 hydrogen bond parameter files.
0.22.5.6    April 13, 2012    sabs
                'toppar/scwrl4hb/par_hbond_CA_1.txt' -Hydrogen bond parameters file for both canonical and noncanonical (CAH) hydrogen
                 bonds.
                'toppar/scwrl4hb/par_hbond_1.txt' -Hydrogen bond parameters file for canonical hydrogen bonds.
                'src/SysEnv.cpp', 'tests/testSysEnv.cpp' -Added MSL_HBOND_CA_PAR
                'library/par_hbond_1.txt' -Moved to toppar directory.
0.22.5.5    April 11, 2012    sabs
                'programs/createEnergyTable.cpp' -Remove loaded rotamers after each environment is evaluated and set environment
                 to crystal rotamer
0.22.5.4    April 06, 2012    bkmueller
                'src/Chain.h', 'src/Chain.cpp' -Fixed bug in renumberChain, now properly iterate through the positionMap
0.22.5.3    April 03, 2012    sabs
                'programs/createEBL.cpp', 'programs/createEnergyTable.cpp' -createEnergyTable creates the energy table from a set
                 of environments. createEBL sorts conformers in a rotamerlibrary using an energytable.
                'src/MslTools.cpp' -Added the charmm names of HIS : HSD,HSE,HSP to threeToOneLetterMap.
                'Makefile' -Added the programs createEBL and createEnergyTable
0.22.5.2    March 21, 2012    sabs
                'programs/getChiRecovery.cpp' -Bug fix: adjusting for angles near the limits 0,360.
0.22.5.1    March 20, 2012    sabs
                'src/CharmmTopologyReader.h' -Added an interface to getMass.
                'programs/getDihedrals.h', 'programs/getDihedrals.cpp' -Added a required option to specify the degrees_of_freedom_file.
                
                'programs/getChiRecovery.cpp' -Trim each line while reading the list of pdbs.
0.22.5.0    March 12, 2012    asenes
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Added removeIdentity functions to remove a certain
                 identity from a Position. This fix some outstanding issue that were present when an atom was deleted, because
                 pointers to it would remain in the bonded information of other atoms, in the ic table and in the table of interactions
                
                'src/Atom.h', 'src/Atom.cpp' -Revised setUnboundFromAll and added removeFromIc, functions needed to cleanup the
                 bond information from other atoms and cleanup the IcTable before this atom is deleted
                'src/IcEntry.h', 'src/IcEntry.cpp' -Added removeAtom function to be called to remove references to Atom pointers
                 before they go out of scope; added isValid function, to allow purging of IcEntries that no longer can be used
                 to build because key atoms are NULL; added getParentTable to call functions to the parent IcTable (again, for
                 removing obsolete pointers)
                'src/Position.h', 'src/Position.cpp' -Updated removeIdentity function
                'src/EnergySet.h', 'src/EnergySet.cpp' -Implemented deleteInteractionsWithAtom functions to remove interactions
                 with atom pointers that are about to go out of scope
                'src/IcTable.h', 'src/IcTable.cpp' -Added removeAtom function to remove the pointer to an atop before it is deleted.
                 Also added a lot of safety check for NULL atoms
                'src/System.h', 'src/System.cpp' -Added purgeIcTable function to remove all non valid ic entries
                'src/Interaction.h', 'src/Interaction.cpp' -Added hasAtom function to check if an interaction contains a certain
                 atom pointer (needed for the removal of interactions that contain atoms that will go out of scope)
                'src/AtomContainer.h', 'src/AtomContainer.cpp' -Fixed small bug in updateAtomMap and removed old commented out
                 code
                'tests/testDeleteBondedAtom.cpp' -A quick test to check if atoms correctly propagate the fact that a bond has been
                 deleted
                'tests/testAddCharmmIdentity.cpp' -Added a third test in which an amino acid type is actually removed
                'Makefile' -Added testDeleteBondedAtom
0.22.4.0    March 05, 2012    asenes
                'src/Atom.h', 'src/Atom.cpp' -Added function getHiddenCoor() to get an AtomPointerVector of the hidden coordinates
                 (needed by Transforms)
                'src/Transforms.h', 'src/Transforms.cpp' -If setTransformAllCoors(true) is given, now transformations are applied
                 to all alternative coordinates including the hidden coordinates
                'src/Residue.h', 'src/Residue.cpp' -Added hideRotamerRelIndex hideRotamerAbsIndex hideAllRotamersButOneRelIndex
                 hideAllRotamersButOneAbsIndex hideAllRotamersButFirstN unhideRotamerAbsIndex unhideAllRotamers, based on the Atom
                 functions that hide
                'tests/testSaveAtomAltCoor.cpp' -Small changes
0.22.3.0    March 01, 2012    asenes
                'src/Atom.h', 'src/Atom.cpp' -Added the possibility of selectingly hiding any of the alternative coordinates. Once
                 hidden, the atom will behave like it does not have them. The intended use is with rotamers, to enable a use to
                 load a large number of rotamers and chose later how many, and even which to use. The new functions are hideAltCoorRelIndex,
                 hideAltCoorAbsIndex, hideAllAltCoorsButOneRelIndex, hideAllAltCoorsButOneAbsIndex, hideAllAltCoorsButFirstN, unhideAltCoorAbsIndex,
                 unhideAllAltCoors. Also, added a flag to change the toString output (setToStringFormat function) to print all
                 alt coors, with our without the hidden (the default priting is the same).
                'tests/testSaveAtomAltCoor.cpp' -Modified the test to test hiding and unhiding of coordinates
0.22.2.13    February 27, 2012    sabs
                'src/System.h' -Fixed residueExists function.Changed header.
0.22.2.12    February 18, 2012    gevorg
                'src/AtomContainer.h', 'src/AtomContainer.h' -AtomContainer can now construct from stringstream.
                'src/PDBReader.cpp' -PDBReader's constructor that constructed from stringstream was calling read()
0.22.2.11    February 17, 2012    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Fixed runGreedyOptimizer : erase history before starting a
                 new optimization run. Updated header in SelfPairManager.h
0.22.2.10    February 16, 2012    sabs
                'programs/repackSideChains.h', 'programs/repackSideChains.cpp' -Added rotlevel functionality. Removed chi stats
                 and sasa printing : use getChiRecovery / getDihedrals instead
0.22.2.9    February 14, 2012    bhannigan
                'scripts/submit.py' -Fixing bug so that submit will now create new directories.
0.22.2.8    February 14, 2012    bhannigan
                'test/test.txt' -Another add test.
0.22.2.7    February 14, 2012    bhannigan
                'test/test.txt' -This is a test
                'test/test.txt' -This is a test
0.22.2.6    February 14, 2012    sabs
                'src/FormatConverter.h', 'src/FormatConverter.cpp' -Added new release header. Removed some unnecessary if conditions.
                
                'programs/getChiRecovery.cpp' -Removed hard coded tolerance value (used to be 20).
0.22.2.5    February 13, 2012    sabs
                'Makefile' -Added getChiRecovery.
0.22.2.4    February 13, 2012    sabs
                'src/DegreeOfFreedomReader.cpp', 'src/DegreeOfFreedomReader.h' -Added getChiAtoms(resName) function
                'src/SysEnv.cpp' -Added MSL_PDB_2.3_DOF and MSL_CHARMM_22_DOF env variables.
                'src/ChiStatistics.h', 'src/ChiStatistics.cpp' -Now uses DegreeOfFreedomReader instead of PDBTopology. Added function
                 getChis.
                'toppar/CHARMM_22_DegOfFreedoms.txt', 'toppar/pdb_2.3_DegOfFreedoms.txt' -Renamed HIS to HSD,HSE,HSP corrected
                 SER HG to HG1. Corrected typo in CYS
                'programs/getChiRecovery.h', 'programs/getChiRecovery.cpp' -Program to compute chi recovery statistics over a set
                 of pdbs.
0.22.2.3    February 11, 2012    asenes
                'toppar/CHARMM_22_DegOfFreedoms.txt' -File with the definitions of degrees of freedom (such as chi1 or phi and
                 psi) for CHARMM 22)
                'programs/repackSideChains.h' -Minor change. Added
                'README.txt' -Reformatted
0.22.2.2    February 11, 2012    asenes
                'src/DegreeOfFreedomReader.h', 'src/DegreeOfFreedomReader.cpp', 'toppar/pdb_2.3_DegOfFreedoms.txt' -New object,
                 returns the atoms in a degree of freedom, for example giving LEU and chi1 it returns N,CA,CB,CG. It reads the
                 definitions from a file: toppar/pdb_2.3_DegOfFreedoms.txt (2nd resubmission after fail)
                'var/header.txt' -Source code file sample header (2nd resubmission after fail)
                'toppar/top_pdb2.3_H.inp', 'toppar/top_pdb2.3_noH.inp' -Charmm topology style files for PDB, NOT containing energetic
                 information, only for building (resubmission after fail)
0.22.2.1    February 11, 2012    asenes
                'src/DegreeOfFreedomReader.h', 'src/DegreeOfFreedomReader.cpp', 'toppar/pdb_2.3_DegOfFreedoms.txt' -New object,
                 returns the atoms in a degree of freedom, for example giving LEU and chi1 it returns N,CA,CB,CG. It reads the
                 definitions from a file: toppar/pdb_2.3_DegOfFreedoms.txt (resubmission after fail)
                'var/header.txt' -Source code file sample header (resubmission after fail)
                'toppar/top_pdb2.3_H.inp', 'toppar/top_pdb2.3_noH.inp' -Charmm topology style files for PDB, NOT containing energetic
                 information, only for building
0.22.2.0    February 10, 2012    asenes
                'src/ConformationEditor.h', 'src/ConformationEditor.cpp' -New object. Allows to edit the conformation of a molecule
                 by entering a series of IC edits, and applying the changes all at once. It can use definitions such as chi1 or
                 phi
                'src/DegreeOfFreedomReader.h', 'src/DegreeOfFreedomReader.cpp', 'toppar/pdb_2.3_DegOfFreedoms.txt' -New object,
                 returns the atoms in a degree of freedom, for example giving LEU and chi1 it returns N,CA,CB,CG. It reads the
                 definitions from a file: toppar/pdb_2.3_DegOfFreedoms.txt
                'README.txt' -New file with MSL credits
                'scripts/submit.py' -Added the executable (/usr/bin/python) at the head of the file
                'tests/testConformationEditor.cpp' -A test for the new ConformationEditor
                'src/PDBTopologyBuilder.h', 'src/PDBTopologyBuilder.cpp' -Now it adds bonded information. It also fills the IC
                 table from the current coordinates
                'var/header.txt' -Source code file sample header
                'Makefile' -Added new objects ConformationEditor and DegreeOfFreedomReader and the testConformationEditor test
                
                'src/CharmmSystemBuilder.cpp' -Now it fills the IC table from the current coordinates
                'src/IcEntry.h', 'src/IcEntry.cpp' -Fixed small bug in fillFromCoor().
                'src/IcTable.h', 'src/IcTable.cpp' -Implemented missing printTable() function
0.22.1.5    February 08, 2012    sabs
                'src/FormatConverter.h', 'src/FormatConverter.cpp', 'tests/testFormatConverter.cpp', 'Makefile' -Added an object
                 for conversion between pdb and charmm names
                'src/Atom.h', 'src/Atom.cpp' -Added copyAllCoor to copy alternate coordinates.
0.22.0.5    February 04, 2012    asenes
                'src/AtomContainer.h' -Bug fix: the readPdb function was not closing the file
0.22.0.4    January 05, 2012    dwkulp
                'Makefile' -move PDBFragments under BOOST flag
                'myProgs/dwkulp/dwkulp.mk', 'myProgs/dwkulp/designLinearSegment.h', 'myProgs/dwkulp/designLinearSegment.cpp', 'myProgs/dwkulp/conformationalSampling.h',
                 'myProgs/dwkulp/conformationalSampling.cpp', 'myProgs/dwkulp/glycineSearch.h', 'myProgs/dwkulp/glycineSearch.cpp'
                 -new programs designLinearSegment checkRMS_EOD conformationalSampling glycineSearch
                'myProgs/dwkulp/generateRotamerLibrary.cpp' -Atom pointer/Atom reference change
                'myProgs/dwkulp/getTripletCaMeasurements.cpp', 'myProgs/dwkulp/getTripletCaMeasurements.h' -Chain broken logic
                 error fixed
                'myProgs/dwkulp/insertSelectionIntoTemplate.cpp' -Using new API: getIndex vs getIndexInSystem
                'myProgs/dwkulp/superRotamerExtraction.cpp', 'myProgs/dwkulp/superRotamerExtraction.h' -includeAllNeighbors flag
                 added
                'programs/createFragmentDatabase.cpp', 'programs/createFragmentDatabase.h' -Use regex to select parts of PDB for
                 fragments
                'src/PDBFormat.cpp' -change ERROR 34918 to print the line of the pdb that got the error!
                'src/PDBFragments.cpp', 'src/PDBFragments.h' -add regex for picking PDBFragments; also implemented a PDB Dir as
                 an alternative to BBQ, so you can go back to original file and pull backbone coordinates out
                'src/PhiPsiStatistics.h', 'src/PhiPsiStatistics.cpp' -added getFreqInBin function
                'programs/renumberResidues.cpp', 'programs/renumberResidues.h' -new program
                'src/PSSMCreator.h', 'src/PSSMCreator.cpp', 'src/FastaReader.h', 'src/FastaReader.cpp' -Create mulitple sequence
                 aligment from set of sequences
0.22.0.3    December 22, 2011    sabs
                'src/HydrogenBondBuilder.cpp' -Avoided building some unnecessary hydrogen bond interactions.
                'library/par_hbond_1.txt' -The OG acceptor line for SER had the HG1 atom wrong (it was entered HG).
0.22.0.2    December 14, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp', 'programs/repackSideChains.h', 'programs/repackSideChains.cpp'
                 -Added code for greedy sidechain optimization.
0.22.0.1    December 13, 2011    sabs
                'src/RotamerLibraryReader.cpp' -Added code to read files with the ICDEF format
                'src/AtomContainer.h' -Added Set and get NameSpace.
                'library/EBL_11-2011_CHARMM22.txt', 'library/EBL_11-2011_PDB2.3.txt', 'library/EBRL_11-2011_CHARMM22.txt', 'library/EBRL_11-2011_PDB2.3.txt'
                 -Energy-based conformer and rotamer libraries
0.22.0.0    December 03, 2011    dwkulp
                'Makefile' -Added FastaReader,PSSMCreator objects, renumberResidues program
                'myProgs/dwkulp/dwkulp.mk', 'myProgs/dwkulp/querySeqCons.cpp', 'myProgs/dwkulp/querySeqCons.h', 'myProgs/dwkulp/calcSasaAll.cpp',
                 'myProgs/dwkulp/calcSasaAll.h', 'myProgs/dwkulp/compareStructures.cpp', 'myProgs/dwkulp/compareStructures.h' -Added
                 compareStrutures, querySeqCons and calcSasaAll programs
                'myProgs/dwkulp/findPositionsWithRotamers.cpp', 'myProgs/dwkulp/findPositionsWithRotamers.h' -added glycan clash
                 check, config file and output file
                'myProgs/dwkulp/multiSearchDM.cpp' -modified to remove binary database code, it is broken for now
                'programs/calculateSasa.cpp' -added a selection option, so you can calculate SASA of only a selection
                'programs/energyOptimizations.h' -re
                'programs/findClashes.cpp', 'programs/findClashes.h' -added more options to this program, for setting a clash tolerance
                 and selecting clashes between specific atom types or element types
                'programs/getSelection.cpp' -fixed usage statement
                'programs/insertLoopIntoTemplate.h', 'programs/insertLoopIntoTemplate.cpp', 'src/FuseChains.h', 'src/FuseChains.cpp'
                 -added includeTemplateStems and checkCaCaDistanes options, includeTemplateStems added to FuseChains object
                'programs/optimizeMC.cpp' -fixed output, System is printed as 1 AA code
                'programs/searchFragmentDatabase.cpp' -setSegID to blank
                'src/PairwiseEnergyCalculator.cpp' -fixed getIndex to getIndexInSystem
                'src/PolymerSequence.cpp' -added Icode support and symmetry support in constructor that takes in variablePositionMap
                 (for design)
                'src/Position.h', 'src/Position.cpp' -setActiveIdentity/setActiveRotamer are moved to the .cpp file and symmetry
                 conditions are more carefully checked. linkedPositions.size() == 0, specificically
                'src/SysEnv.cpp' -Added MSL_EBL variable by default, which sets the EnergyBasedLibrary path
                'src/System.h', 'src/System.cpp' -added findProteinInterfacePositions(chain1,chain2) function
                'src/SystemRotamerLoader.cpp' -setActiveIdentity call shouldn't apply to linked positions, so now it does not;
                 also MSLOUT support added
                'src/FastaReader.cpp', 'src/FastaReader.h', 'src/PSSMCreator.cpp', 'src/PSSMCreator.h' -new objects
0.21.0.2    November 21, 2011    dwkulp
                'myProgs/dwkulp/multiSearchDM.cpp' -update to work with current MSL
                'myProgs/dwkulp/multiSearchDM.cpp' -update to work with current MSL
0.21.0.1    November 15, 2011    sabs
                'programs/repackSideChains.h' -Corrected a log message.
                'src/SelfPairManager.cpp', 'src/SelfPairManager.h', 'src/EnergySet.cpp', 'src/EnergySet.h' -Support for term weights.
                
0.21.0.0    November 08, 2011    dwkulp
                'myProgs/myProgs.mk.RENAME_ME' -new myProgs/dwkulp dir
0.20.3.0    November 07, 2011    sabs
                'Makefile', 'src/SidechainOptimizationManager.h', 'src/SidechainOptimizationManager.cpp' -Added the SidechainOptimizationManager
                 object which will replace SelfPairManager eventually. For now, it is only a copy.
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Removed code for self energy cutoff. Added onTheFly energy
                 computation.Replaced calls to getNumberOfAltConformations with getNumberOfRotamers so that rotamer limits may
                 be used.
                'src/Position.h' -Enabled setting number of rotamers by level.Replaced calls to getNumberOfAltConformations with
                 getNumberOfRotamers so that rotamer limits may be used.
                'src/RotamerLibraryReader.cpp', 'src/RotamerLibraryWriter.h', 'src/RotamerLibraryWriter.cpp' -Added code to read
                 rotamer sampling levels in LEVRES and LEVEL lines.
                'src/MonteCarloOptimization.cpp', 'src/MonteCarloOptimization.h' -Now takes a pointer to SelfPairManager to compute
                 energies in OnTheFly mode.
                'src/SystemRotamerLoader.h', 'src/SystemRotamerLoader.cpp' -loadRotamers can now take the sampling level as a string.
                 The defineRotamerSamplingLevels method assigns to each residue, the number of rotamers at each sampling level
                 based on the RotamerLibrary.
                'src/RotamerLibrary.h', 'src/RotamerLibrary.cpp' -Added code to support rotamer levels.
                'src/Residue.h', 'src/Residue.cpp', 'src/Chain.h', 'src/System.h' -Added code to support rotamer levels and maxRotamers.
                
                'src/MslTools.h', 'src/MslTools.cpp', 'src/OptionParser.cpp', 'src/OptionParser.h' -Added methods to convert unsigned
                 ints.
                'programs/repackSideChains.cpp', 'programs/repackSideChains.h' -Added onTheFly support and removed selfEnergyCutOff.
                
0.20.2.0    November 01, 2011    dwkulp
                'programs/createFragmentDatabase.cpp', 'programs/searchFragmentDatabase.cpp', 'programs/searchFragmentDatabase.h'
                 -update to include regex,sasa and sse criteria
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp' -get functions for maximum and minimum number of alternate
                 conformations
                'src/AtomSelection.h', 'src/AtomSelection.cpp' -helper function to allow inverseSelections, this is easier/cleaner
                 than users having to create 'not' strings
                'programs/backrubPdb.cpp', 'programs/backrubPdb.h', 'src/BackRub.cpp' -skip PRO residues, added MSLOUT, distance
                 check, new backrubPdb program
                'src/Frame.h', 'src/Frame.cpp' -compute frames from points function added
                'src/MslTools.h', 'src/MslTools.cpp' -stringf function added. stringf is the equivalent of sprintf but for std::string
                
                'src/PDBTopology.h', 'src/PDBTopology.cpp' -getGenericResidue function added. Builds any residue type in a default
                 orientation, and will create numRotamers number of rotamers in that orientation
                'src/PhiPsiStatistics.cpp', 'src/PhiPsiStatistics.h' -getOmega utility function added
                'src/RegEx.cpp', 'src/RegEx.h' -added string type, now you can search by segment_id or primary sequence. I commonly
                 put secondary structure info in segment_id field
                'src/Residue.cpp' -neighbors function, added an optional skippable atom set when detecting neighbors
                'src/System.h', 'src/System.cpp' -added writeMulitplePdbs, getSizes, getPDBReader, getPDBWriter functions
                'src/Transforms.h', 'src/Transforms.cpp' -added smartRmsdAlignment function that will fish through the atomvectors
                 to find a set of atoms to align; added an enum called MatchType to support smartRmsdAlignment
                'src/VectorHashing.h', 'src/VectorHashing.cpp', 'src/VectorPair.h', 'src/VectorPair.cpp' -updates to still
                'Makefile' -removed hydrogen bond builder test, added MSL_STATIC (defaulted to true) for building static vs shared
                 library, added R compile code to help with compliation (removed hardcoded paths).
0.20.1.0    October 31, 2011    sabs
                'src/RotamerLibraryBuilder.cpp', 'src/RotamerLibraryReader.cpp', 'src/RotamerLibrary.h', 'src/RotamerLibrary.cpp',
                 'src/RotamerLibraryBuilder.h' -Added an optional rotamer bin column at the end of each CONF line. This could represent
                 the dunbrack02 bin for each rotamer. Default bin is 0 and if the bin value is absent it is assumed to be zero.
                
                'src/MslTools.h', 'src/MslTools.cpp' -Added a toUnsignedInt function to convert string to unsigned int.
0.20.0.10    October 27, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Changed the default Monte Carlo Parameters. Added an interface
                 to perform LP/MIP
                'src/CharmmSystemBuilder.cpp' -BugFix: Ensured that both C and N termini are patched for a single residue chain
                
                'src/RotamerLibrary.cpp' -BugFix: RemoveRotamer now works with _libName=
                'src/Atom.h' -Changed %u to %d for printing residue number
                'programs/repackSideChains.h', 'programs/repackSideChains.cpp' -Changed default number of rotamers to match dunbrack
                 numbers. Made runUnbiasedMC true by default. Bug: cuton,cuthb,minDeltaE are now set properly. Not using SelfECutoff
                 anymore.
0.20.0.9    September 26, 2011    asenes
                'src/CartesianGeometry.cpp' -Removed debugging cout statement that I forgot in the build() function
0.20.0.8    September 23, 2011    asenes
                'src/CartesianGeometry.cpp' -Function buildRadians(): renames some variables, added comments, no actual change
                 of code
0.20.0.7    September 02, 2011    sabs
                'programs/repackSideChains.cpp' -Bug System::buildAllAtoms should be called before CharmmSystemBuilder::updateNonBonded.
                
                'src/Atom3DGrid.cpp' -Bug add another bin to catch overflow due to loss of precision. for eg. int bin = int(6.9999999999999999999999500)
                 turns out to be 7 instead of 6.
0.20.0.6    August 19, 2011    bhannigan
                'src/SurfaceAreaAndVolume.cpp', 'src/SurfaceAreaAndVolume.h', 'tests/testSurfaceAreaAndVolume.cpp' -Fixed some
                 bugs in the analytical calculation of SASA and excluded volume. Also added a few more tests to the test suite
                 for this functionality.
0.20.0.5    August 17, 2011    sabs
                'src/EnergySet.cpp' -Changed behaviour of saveEnergySubset. Now, all terms will be saved while the subset is created,
                 only those that are turned on at the time of calcEnergySubset will be considered for energy calculation. Also,
                 getTermEnergy() returns 0.0 instead of -1.0 when a term does not exist.
                'tests/testSasaCalculator.cpp', 'src/SasaCalculator.cpp' -Fixed the way sasaTable was created when byAtom flag
                 is false.
0.20.0.4    July 13, 2011    sabs
                'programs/repackSideChains.cpp' -Bug: Calculate Chi2 and altChi2 correctly for all residues.
0.20.0.3    July 11, 2011    sabs
                'Makefile', 'programs/repackSideChains.h', 'programs/repackSideChains.cpp' -Added options to excludeEnergyTerms.
                 Bug: Chi values need to be measured before rotamers are loaded.
0.20.0.2    July 11, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added method in SelfPairManager to set MCO params
                'programs/repackSideChains.cpp', 'programs/repackSideChains.h' -Program to repack sidechains
0.20.0.1    July 05, 2011    sabs
                'tests/testRotamerOptimization.cpp' -Use updated interfaces in MonteCarloOptimization Object.
                'src/LinearProgrammingOptimization.h', 'src/LinearProgrammingOptimization.cpp' -Works with lower triangular energy
                 table.
                'src/MonteCarloOptimization.cpp', 'src/MonteCarloOptimization.h' -Uses MonteCarloManger to implement the MonteCarlo
                 Search. Some interface have been updated.
                'src/SelfPairManager.cpp' -Print out only when verbose flag is set.
0.20.0.0    June 25, 2011    asenes
                'src/Position.h' -Changes in linked positions. Linked positions are now created only by acting on the MASTER using
                 the addLinkedPosition(Position &_pos) function. No need to set the type with setLinkedPositionType (removed).
                 Now setActiveIdentity and setActiveRotamer when applied to a MASTER position also apply by default also to the
                 SLAVE linked positions.
                'src/System.h', 'src/System.cpp' -Adjusted for linked positions. Added getMasterPositions to get a list of all
                 MASTER linked positions, and getLinkedPositionType(index) to check the type of position at a certain index. The
                 setLinkedPositions now requires a proper positionId.
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added support for linked positions. Linked positions are treated
                 as a single position in the calculations. These position need to have the same identities and number of rotamers
                 or the SelfPairManager will exit with an error
                'tests/testLinkedPositions.cpp' -Rewrote from scratch the test to make sure it works with linked positions and
                 the SelfPairManager
0.19.1.2    June 25, 2011    sabs
                'src/MonteCarloOptimization.cpp', 'src/MonteCarloOptimization.h' -Uses MonteCarloManager to runMC.
                'src/RandomNumberGenerator.cpp' -Fixed a memory leak
                'src/MonteCarloManager.cpp' -changed the SIGMOIDAL and EXPONENTIAL cooling functions
                'src/SelfPairManager.cpp', 'src/SelfPairManager.h' -Added an unbiased MC search using the MonteCarloOptimization
                 object.Made storing interaction counts optional to save some memory.
                'src/SelfConsistentMeanField.cpp' -Made SCMF biased MC search part of the SCMF object.
                'programs/energyOptimizations.h', 'programs/optimizeMC.cpp' -Updated the program to use the new methods in MonteCarloOptimization
                
                'examples/example_pdbfrag.cpp' -Fixed a typo.
0.19.1.1    June 24, 2011    sabs
                'src/SystemRotamerLoader.cpp' -Fixing a bug in SystemRotamerLoader::loadRotamers. If _keepOldRotamers is false,
                 all existing rotamers must be removed. This bug was introduced in version 0.18.7.2
0.19.1.0    June 24, 2011    asenes
                'src/Transforms.h', 'src/Transforms.cpp', 'tests/testTransformBondAngleDiheEdits.cpp' -Added a flag (function setNaturalMovements)
                 for the setBond, setAngle, setImproper and setDihedral functions, it rotates/translates both sides but propotionally
                 to the number of atoms and movements so that the change is more natural
0.19.0.0    June 14, 2011    dwkulp
                'examples/example_backrub.cpp', 'examples/example_ccd.cpp', 'examples/example_pdbfrag.cpp', 'examples/examples.mk'
                 -backbone motion algorithms
                'programs/generateCoiledCoils.cpp' -doesn't compile anymore because of CoiledCoil API changes, commented out for
                 now
                'src/BackRub.cpp', 'src/BackRub.h', 'src/CCD.h', 'src/CCD.cpp' -API changes to store the Atoms internally to the
                 object
                'src/PDBFragments.h', 'src/PDBFragments.cpp' -Change API to use a System, return system with getSystem, some slight
                 changes to make sure this could be used without BBQ
                'src/PDBReader.h', 'src/PDBReader.cpp' -analogous to PDBWriter, the read function can be told to skip the hydrogens
                
                'src/SysEnv.cpp' -Added BBQ_TABLE, PDB_FRAG_TABLE varirables
                'src/System.h', 'src/System.cpp' -Added writeAllModels the writePdb function which goes through all the models
                 and writes them out as an NMR
                'src/Transforms.cpp' -clean up old comments about MSLOUT
                'src/CharmmEEF1Interaction.h' -did not compile, needed a local distance defined
                'tests/testBBQ.cpp' -Uses SYSENV now
                'Makefile' -testBoost was added to tests
                'exampleFiles/example0008.pdb', 'exampleFiles/example0008_caOnly.pdb' -example pdbs for example_ccd example_backrub
                 and example_pdbfrag
0.18.8.1    June 13, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfConsistentMeanField.cpp', 'src/SelfConsistentMeanField.h', 'src/SpringConstraintInteraction.cpp',
                 'src/SelfPairManager.cpp', 'programs/coiledCoilBuilder.cpp' -Moved SCMF based monte carlo into the SCMF object.
                 Added flags to run unbiased MC and SCMF biased MC. Updated interfaces in SelfPairManager and coiledCoilBuilder.cpp.
                
0.18.7.3    June 07, 2011    sabs
                'src/CharmmImproperInteraction.h', 'src/UserDefinedInteraction.cpp', 'src/SpringConstraintInteraction.h', 'src/CharmmElectrostaticInteraction.h',
                 'src/CharmmVdwInteraction.cpp', 'src/CharmmEEF1RefInteraction.h', 'src/CharmmAngleInteraction.h', 'src/BaselineInteraction.cpp',
                 'src/CharmmImproperInteraction.cpp', 'src/BaselineInteraction.h', 'src/CharmmElectrostaticInteraction.cpp', 'src/Scwrl4HBondInteraction.cpp',
                 'src/CharmmVdwInteraction.h', 'src/SpringConstraintInteraction.cpp', 'src/CharmmDihedralInteraction.h', 'src/CharmmAngleInteraction.cpp',
                 'src/CharmmUreyBradleyInteraction.h', 'src/Scwrl4HBondInteraction.h', 'src/Interaction.h', 'src/UserDefinedInteraction.h',
                 'src/CharmmDihedralInteraction.cpp', 'src/CharmmBondInteraction.cpp', 'src/CharmmEEF1Interaction.h', 'src/CharmmBondInteraction.h'
                 -Added a getEnergy(vector<double>* _gradient) function which will be used to compute gradients during minimization.
                 If the pointer is NULL this function computes the energy without the switching function even if cutoffs are in
                 place.
                'src/EnergySet.cpp', 'src/CharmmEnergy.cpp', 'src/CharmmEnergy.h', 'src/GSLMinimizer.cpp', 'src/EnergySet.h' -The
                 gradient computation for minimization doesnot consider the switching function. A new method calcEnergyWithoutSwitchingFunction
                 is added to EnergySet to be able to compute energies without applying the switching function. Removed code from
                 CharmmEnergy that computed gradients with the switching function.
                'programs/minimize.cpp' -Reads env variables for default topology and parameter files.
0.18.7.2    June 01, 2011    sabs
                'src/EnergySet.cpp', 'src/SystemRotamerLoader.cpp', 'src/Interaction.h' -Fixed a bug in SystemRotamerLoader. Added
                 notes to combine the partialDerivative and getEnergy functions in the interaction hierarchy.
0.18.7.1    May 26, 2011    sabs
                'tests/testNonBondedCutoff.cpp', 'tests/testCharmmBuild.cpp', 'tests/testMinimization.cpp', 'tests/testEEF1_2.cpp',
                 'tests/testRegEx.cpp', 'tests/testVectorPair.cpp', 'tests/testDerivatives.cpp' -Fixed compile issues by using
                 the updated methods in System,CharmmSystemBuilder, etc.,
                'src/CharmmImproperInteraction.h', 'src/CharmmElectrostaticInteraction.h', 'src/CharmmVdwInteraction.cpp', 'src/CharmmEEF1RefInteraction.h',
                 'src/CharmmUreyBradleyInteraction.cpp', 'src/CharmmAngleInteraction.h', 'src/BaselineInteraction.cpp', 'src/CharmmImproperInteraction.cpp',
                 'src/BaselineInteraction.h', 'src/CharmmElectrostaticInteraction.cpp', 'src/Scwrl4HBondInteraction.cpp', 'src/CharmmVdwInteraction.h',
                 'src/CharmmEEF1RefInteraction.cpp', 'src/CharmmEEF1Interaction.cpp', 'src/CharmmDihedralInteraction.h', 'src/CharmmAngleInteraction.cpp',
                 'src/CharmmUreyBradleyInteraction.h', 'src/Scwrl4HBondInteraction.h', 'src/Interaction.h', 'src/UserDefinedInteraction.h',
                 'src/CharmmDihedralInteraction.cpp', 'src/CharmmBondInteraction.cpp', 'src/CharmmEEF1Interaction.h', 'src/CharmmBondInteraction.h'
                 -typeName is back to being a const static member. partialDerivative function has moved from EnergySet to the Interaction
                 hierarchy.
                'src/GSLMinimizer.h', 'src/GSLMinimizer.cpp', 'programs/minimize.h', 'programs/minimize.cpp' -INTERFACES TO GSLMinimizer
                 HAVE CHANGED. GSLMinimizer can now be initialized using a reference to a System. Changed the Minimize function
                 to minimize.Constrained and fixed atoms can be specified using selections. Cycles of constrained minimization
                 may be performed where the spring is reset after each cycle. The minimize program has been modified to reflect
                 this
                'src/CartesianGeometry.cpp', 'src/CharmmSystemBuilder.cpp', 'src/EnergySet.cpp', 'src/HydrogenBondBuilder.cpp',
                 'src/Atom.h', 'src/EnergySet.h', 'programs/getDihedrals.cpp', 'Makefile' -EnergySet::resetTerm is renamed to eraseTerm.
                 Added comment in Atom.h to note that minimization index should start from 1. Removed tests that dont compile from
                 the Makefile
                'src/SpringConstraintInteraction.h', 'src/SpringConstraintInteraction.cpp' -Interaction to model the spring in
                 constraint minimization.
0.18.6.2    May 26, 2011    asenes
                'src/SysEnv.cpp' -Added an MSL_EXAMPLE_FILE_DIR variable
                'tests/testSysEnv.cpp' -Updated the test with the MSL_EXAMPLE_FILE_DIR variable
0.18.6.1    May 26, 2011    asenes
                'src/SysEnv.cpp' -Added an MSL_PDB_TOP variable
                'tests/testSysEnv.cpp' -Updated the test with the MSL_PDB_TOP variable
0.18.6.0    May 25, 2011    asenes
                'src/SysEnv.h', 'src/SysEnv.cpp' -Changed the name of the enviromental variables for paratmers and topology files
                 to have an MSL prefix, also reduced the number: MSL_CHARMM_TOP, MSL_CHARMM_PAR, MSL_ROTLIB, MSL_HBOND_PAR
                'src/Quench.cpp', 'tests/testNonBondedCutoff.cpp', 'tests/testQuench.cpp', 'tests/testCharmmTopologyReader.cpp',
                 'tests/testCharmmBuild.cpp', 'tests/testEnergeticAnalysis.cpp', 'tests/testEEF1.cpp', 'tests/testSysEnv.cpp',
                 'tests/testCharmmEnergies.cpp', 'tests/testLinkedPositions.cpp', 'tests/testEEF1_2.cpp', 'tests/testPDBTopology.cpp',
                 'tests/testSurfaceAreaAndVolume.cpp' -Changed the name of the enviromental variables for paratmers and topology
                 files to MSL_CHARMM_TOP, MSL_CHARMM_PAR, MSL_ROTLIB, MSL_HBOND_PAR
                'src/System.h' -The buildAllAtoms function now takes backbone atoms for a copy (defaulted to N CA C O HN, charmm
                 names), no need to copy bb cordinates over anymore
                'tests/testAddCharmmIdentity.cpp' -Put it in line with the new building functionalities of the system (seed, buildAllAtoms)
                
                'src/Position.cpp', 'src/Position.h' -The copyCoordinatesOfAtoms now by default does not copy unless the atoms
                 has no coordinates (added a bool flag for it)
                'src/CharmmSystemBuilder.h', 'src/CharmmSystemBuilder.cpp' -Removed addIdentity functions that take the list of
                 backbone atoms as a vector
                'src/SelfPairManager.cpp' -Fixed a bug in the calculation of the self energies. Corrected the printing for when
                 the individual terms are not saved
                'src/GSLMinimizer.cpp' -Reintroduced the offending GSL variables that require v.1.14 but put their function compilation
                 under a new enviromental variable MSL_GSL_OLD
                'Makefile' -Added an environmental variable to flag that GSL is old, i.e. prior v.1.14, to fix a problem with the
                 minimizer not compiling, MSL_GSL_OLD
0.18.5.1    May 25, 2011    jedonald
                'src/Quench.cpp' -Do not double count template energies
                'src/AtomicPairwiseEnergy.cpp' -Do not double count template energies
                'src/Transforms.cpp' -Add a function to revert transformations
                'src/TwoBodyDistanceDependentPotentialTable.cpp' -Do not double count template energies
                'src/Transforms.h' -Add a function to revert transformations
                'src/GSLMinimizer.cpp' -Comment out unrecognized random number generators from GSL
0.18.5.0    May 23, 2011    dwkulp
                'src/Quench.cpp', 'tests/testAddCharmmIdentity.cpp', 'tests/testCharmmBuild.cpp', 'tests/testCharmmEnergies.cpp',
                 'tests/testCharmmTopologyReader.cpp', 'tests/testEEF1.cpp', 'tests/testEEF1_2.cpp', 'tests/testEnergeticAnalysis.cpp',
                 'tests/testLinkedPositions.cpp', 'tests/testNonBondedCutoff.cpp', 'tests/testPDBTopology.cpp', 'tests/testPolymerSequence.cpp',
                 'tests/testQuench.cpp', 'tests/testSurfaceAreaAndVolume.cpp', 'Makefile', 'src/SysEnv.h', 'src/SysEnv.cpp', 'tests/testSysEnv.cpp'
                 -Implementation of user environement; SysEnv stores environment variables with defaults for use throughout MSL
                
0.18.4.1    May 23, 2011    asenes
                'src/AtomSelection.h', 'src/AtomSelection.cpp' -Added selectionSize() function to get the number of atoms in a
                 selection
                'src/IcTable.h', 'src/IcTable.cpp' -Fixed bug, now seed() will attempt to seed multiple chains, not just the first
                 one
                'src/System.h', 'src/System.cpp' -Added the missing seed() function, also inlined all seed functions
                'src/CharmmSystemBuilder.h' -Now the addIdentity has default backbone atoms (when they are given in string format
                 only). The default is CHARMM22 default N CA C O HN
                'src/PDBTopologyBuilder.h' -Now the addIdentity has default backbone atoms (when they are given in string format
                 only). The default is PDB 2.3 default N CA C O H
                'tests/testPDBTopologyBuild.cpp' -Addded test for PDBTopologyBuilder
                'Makefile' -Added test testPDBTopologyBuild
0.18.4.0    May 21, 2011    asenes
                'toppar/top_pdb2.3_H.inp' -Topology file with PDB v.2.3 atom names, including hydrogen atoms
                'toppar/top_pdb2.3_noH.inp' -Topology file with PDB v.2.3 atom names, without the hydrogen atoms
                'exampleFiles/example0000.pdb', 'exampleFiles/example0001.pdb' -Fixed terminal oxygen names (O and OXT instead
                 of OT1 OT2)
                'src/PDBTopologyBuilder.h', 'src/PDBTopologyBuilder.cpp' -Derived from CharmmSystemBuilder, this objects can create
                 molecules with PDB naming conventions (using top_pdb2.3_H.inp and top_pdb2.3_noH.inp), and can also add identies
                 to a pre
                'examples/example_add_identity_to_position.cpp' -Example, adding a new identity to a PDB using the PDBTopologyBuilder.cpp
                
                'examples/examples.mk' -Updated with example_add_identity_to_position.cpp
                'Makefile' -Added the PDBTopologyBuilder
0.18.3.0    May 19, 2011    asenes
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp', 'src/AtomContainer.h', 'src/Residue.h', 'src/Position.h',
                 'src/Chain.h', 'src/System.h' -Added support to save and restore alt coor to all molecular containers
                'src/RotamerLibrary.h', 'src/RotamerLibrary.cpp', 'src/RotamerLibraryReader.cpp', 'src/RotamerLibraryWriter.cpp',
                 'src/RotamerLibraryBuilder.cpp', 'src/SystemRotamerLoader.cpp', 'src/PDBTopology.cpp', 'tests/testRotamerLibraryWriter.cpp'
                 -Changed the tag for the atoms that need to be rebuilt from INIT to the more understandable MOBI (for mobile).
                 The old format of the rotamer library file with INIT is still read, the writer uses MOBI. The functions such as
                 setInitAtoms and getInitAtoms have been renamed setMobileAtoms (the old ones are still active but deprecated)
                
                'Makefile' -Remove duplicated test program entry, removed all remaining references to the objs/flags file
0.18.2.0    May 16, 2011    asenes
                'programs/Minimize.h', 'programs/Minimize.cpp' -Removed, renamed as minimize.h minimize.cpp as programs are normally
                 lowercase
                'programs/minimize.cpp', 'programs/minimize.h' -Renamed from formerly uppercase Minimize
                'src/Position.h' -Added getTotalNumberOfRotamers of an identity (by index or identity name)
                'src/Transforms.h', 'src/Transforms.cpp' -Now align and orient return a bool value. Removed unused functions RotatePdbAboutZYX
                 and TranslateRigidBodyPdbResidue. Added a setter for transforming also the alt coors when a rotation
                'src/CoiledCoils.h', 'src/CoiledCoils.cpp' -Deprecated function primarySequenceToCoiledCoil and renamed to setSystemToCoiledCoil.
                 Allowed to change the default backbone atom names with setters. Removed a unused System pointer. Added check for
                 existance of atoms during coiled coil building. Removed commented out functions.
                'programs/coiledCoilBuilder.cpp' -Changed primarySequenceToCoiledCoil function to setSystemToCoiledCoil, and reactivated
                 write pdb function
                'src/Atom.h', 'src/Atom.cpp' -Removing all alt conf keeps the current conf, not the first. Now we can save to buffer
                 also a whole set of alt coor (with the saveAltCoor function). Fixed also a bug, saving coor could have potentially
                 generated memory leaks.
                'src/AtomPointerVector.h', 'src/AtomPointerVector.cpp' -Added the saveAltCoor function.
                'src/System.h' -Added the saveAltCoor function. Also added a function to the the total number of rotamers by identity.
                
                'programs/energyOptimizations.h' -Removed hard coded name of the rotamer library
                'tests/testSaveAtomAltCoor.cpp' -Added test for saving to buffer all alt coors
                'tests/testAtomAndResidueId.cpp' -Added a check, when the identity is asked if it exists, the position can also
                 be obtained with a lastFoundPosition
                'Makefile' -Added and removed objects and programs
0.18.1.0    May 14, 2011    asenes
                'src/System.h', 'src/System.cpp' -Added seed() function without atom pointers, it looks for atoms in the IC table
                 by itself. Removed unused seed function with chainid, resnum, icode and atom name. Moved the functionality of
                 the seed functions into the IcTable.
                'src/IcTable.h', 'src/IcTable.cpp' -Added seeding functionality from the System, added seed() function that does
                 not take atoms, searches seedable atoms on its own
                'src/IcEntry.h', 'src/IcEntry.cpp' -Added some functions to check if atoms are those that specify distance1 and
                 2 and angle 1 and 2. Removed MSLOUT debugging output that was on by default
                'tests/testIcBuilding.cpp' -Revised to test the new seed() function
0.18.0.2    May 13, 2011    sabs
                'tests/testMinimization.cpp' -Modified to work with new interface for GSLMinimizer
                'src/GSLMinimizer.h', 'src/GSLMinimizer.cpp' -GSLMinimizer::Minimize() returns a bool.
                'src/HydrogenBondBuilder.cpp' -Commented out debug output.
                'programs/Minimize.cpp', 'programs/Minimize.h' -Added parameter tolerance.
                'library/par_hbond_1.txt' -Initial version of hydrogen bond parameter file.
0.18.0.1    May 11, 2011    sabs
                'src/CharmmImproperInteraction.h', 'src/CharmmAngleInteraction.h', 'src/CharmmDihedralInteraction.h', 'src/EnergySet.cpp'
                 -Changed getEnergy function to get input angle in radians instead of degrees.
                'programs/Minimize.h', 'programs/Minimize.cpp' -Added stepsize parameter and removed cycles parameter.
0.18.0.0    May 11, 2011    sabs
                'Makefile' -Changed Minimization to the correct file name Minimize
                'programs/Minimize.h', 'programs/Minimize.cpp' -Added options to set number of cycles and the algorithm to be used.
                
                'scripts/submit.py' -Added code to enable file deletion.
                'src/EnergySet.h', 'src/EnergySet.cpp' -Added the method getEnergyTermsSubsets.Updated calls to distanceDerivative,angleDerivative
                 and dihedralDerivative to use pointers.
                'src/GSLMinimizer.h', 'src/GSLMinimizer.cpp' -Added constraint minimization.
                'src/CharmmImproperInteraction.h', 'src/UserDefinedInteraction.cpp', 'src/CharmmElectrostaticInteraction.h', 'src/CharmmVdwInteraction.cpp',
                 'src/CharmmEEF1RefInteraction.h', 'src/CharmmUreyBradleyInteraction.cpp', 'src/CharmmAngleInteraction.h', 'src/FourBodyInteraction.h',
                 'src/Interaction.cpp', 'src/BaselineInteraction.cpp', 'src/CharmmImproperInteraction.cpp', 'src/BaselineInteraction.h',
                 'src/CharmmElectrostaticInteraction.cpp', 'src/CharmmVdwInteraction.h', 'src/CharmmEEF1RefInteraction.cpp', 'src/CharmmEEF1Interaction.cpp',
                 'src/TwoBodyInteraction.h', 'src/CharmmDihedralInteraction.h', 'src/CharmmAngleInteraction.cpp', 'src/CharmmUreyBradleyInteraction.h',
                 'src/Interaction.h', 'src/UserDefinedInteraction.h', 'src/OneBodyInteraction.h', 'src/CharmmDihedralInteraction.cpp',
                 'src/CharmmBondInteraction.cpp', 'src/CharmmEEF1Interaction.h', 'src/CharmmBondInteraction.h', 'src/ThreeBodyInteraction.h'
                 -Removed members distance/angle and energy from the Interaction hierarchy. Changed the typeName member from a
                 static const to an ordinary member. Fixed the getEnergy functions of vdw,elec,eef1 interactions
                'src/Scwrl4HBondInteraction.h', 'src/Scwrl4HBondInteraction.cpp', 'src/HydrogenBondBuilder.h', 'src/HydrogenBondBuilder.cpp'
                 -Preprocess to create lists of acceptors and donors. This fixes a bug in the previous implementation. Added the
                 update option to hydrogen bond builder. Removed hard coded values from Scwrl4HBondInteraction.
                'src/CartesianGeometry.h', 'src/CartesianGeometry.cpp' -All derivative functions take a vector<double>* arguement.
                
                'src/MslTools.cpp' -Removed unnecessary output.
                'src/CharmmSystemBuilder.cpp' -Removed repeated assignment.
                'src/CharmmEnergy.cpp', 'src/CharmmEnergy.h' -Gradient computation with nonbonded cutoffs.
0.17.0.1    May 09, 2011    sabs
                'src/CharmmImproperInteraction.h', 'src/UserDefinedInteraction.cpp', 'src/CharmmElectrostaticInteraction.h', 'src/CharmmEEF1RefInteraction.h',
                 'src/CharmmAngleInteraction.h', 'src/FourBodyInteraction.h', 'src/BaselineInteraction.cpp', 'src/BaselineInteraction.h',
                 'src/Scwrl4HBondInteraction.cpp', 'src/CharmmVdwInteraction.h', 'src/TwoBodyInteraction.h', 'src/CharmmDihedralInteraction.h',
                 'src/CharmmUreyBradleyInteraction.h', 'src/Scwrl4HBondInteraction.h', 'src/Interaction.h', 'src/UserDefinedInteraction.h',
                 'src/OneBodyInteraction.h', 'src/CharmmEEF1Interaction.h', 'src/CharmmBondInteraction.h', 'src/ThreeBodyInteraction.h'
                 -Fixed bug in getEnergy function of Vdw,Electrostatics and EEF1 interactions. Also changed the signature of getEnergy
                 function, instead of a reference to a double, it takes a double.
0.17.0.0    April 29, 2011    dwkulp
                'Makefile', 'tests/testMinimization.cpp', 'programs/Minimize.h', 'programs/Minimize.cpp' -Added tests/testMinimization
                 and programs/Minimize
                'src/BaselineInteraction.cpp', 'src/BaselineInteraction.h', 'src/CartesianGeometry.h', 'src/CartesianGeometry.cpp',
                 'src/CharmmAngleInteraction.cpp', 'src/CharmmAngleInteraction.h', 'src/CharmmBondInteraction.cpp', 'src/CharmmBondInteraction.h',
                 'src/CharmmDihedralInteraction.h', 'src/CharmmDihedralInteraction.cpp', 'src/CharmmEEF1Interaction.h', 'src/CharmmEEF1RefInteraction.h',
                 'src/CharmmElectrostaticInteraction.cpp', 'src/CharmmElectrostaticInteraction.h', 'src/CharmmEnergy.h', 'src/CharmmEnergy.cpp',
                 'src/CharmmImproperInteraction.cpp', 'src/CharmmImproperInteraction.h', 'src/CharmmUreyBradleyInteraction.h',
                 'src/CharmmVdwInteraction.cpp', 'src/CharmmVdwInteraction.h', 'src/EnergySet.h', 'src/EnergySet.cpp', 'src/FourBodyInteraction.h',
                 'src/GSLMinimizer.h', 'src/GSLMinimizer.cpp', 'src/Interaction.h', 'src/OneBodyInteraction.h', 'src/Scwrl4HBondInteraction.cpp',
                 'src/Scwrl4HBondInteraction.h', 'src/ThreeBodyInteraction.h', 'src/TwoBodyInteraction.h', 'src/UserDefinedInteraction.cpp',
                 'src/UserDefinedInteraction.h' -Energy Minimization edits, mostly adding energy gradient related functions
0.16.5.1    April 22, 2011    dwkulp
                'Makefile' -new MSL_MSLOUT_DEBUG_OFF flag, default behavior is for debug information to be turned off
                'scripts/mslBuildTools.py' -number of cores to build is defaulted to 1
                'src/MslOut.h' -added debug() stream
                'src/OptionParser.cpp' -default to MslOut output to be turned off
                'tests/testMslOut2.cpp' -added some new ways to put MslOut output in objects (WARNINGS,ERRORS) and a performance
                 test
0.16.5.0    April 20, 2011    dwkulp
                'Makefile' -testSharedPointers mutate findClashes generateCoiledCoils
                'programs/getDihedrals.cpp', 'programs/getDihedrals.h' -added a SASA calculator, so in addition to dihedral angles
                 it prints out SASA/deltaSASA
                'programs/getSelection.cpp', 'programs/getSelection.h' -added a length option to print out the number of residues
                 that matched a given selection
                'src/Atom.h' -formated toString function so it is readable
                'src/CartesianGeometry.h' -found extra character in GPL statement
                'src/IcEntry.cpp' -MSLOUT functionality added
                'src/MslTools.cpp', 'src/MslTools.h' -Added parseRotamerId functionality
                'src/PDBTopology.cpp', 'src/PDBTopology.h' -Added protonated HIS (HSP) , fixed L/D improper dihedrals coming from
                 IcEntry constructor
                'src/Position.h' -Added getRotamerId function
                'src/PyMolVisualization.h', 'src/PyMolVisualization.cpp' -createCylinder function uses double not int for rgb values
                 and added a random name generator for naming objects
                'src/RotamerLibraryBuilder.h', 'src/RotamerLibraryBuilder.cpp' -added an additional addRotamer function that can
                 add a rotamer to a new library and get DEFI values from an old one
                'src/SasaCalculator.h', 'src/SasaCalculator.cpp' -getTotalSasa function added
                'src/SurfaceAreaAndVolume.cpp' -added the ARVO algorithm reference
                'src/VectorHashing.cpp', 'src/VectorHashing.h' -new files
                'src/VectorPair.cpp', 'src/VectorPair.h' -new files
                'tests/testSharedPointers2.cpp' -new files
                'tests/testVectorPair.cpp' -new files
                'programs/findClashes.h', 'programs/findClashes.cpp' -new files
                'programs/mutate.h', 'programs/mutate.cpp' -new files
0.16.4.5    April 20, 2011    brettth
                'scripts/submit.py' -This is the script for submitting code.
0.16.4.4    April 08, 2011    brettth
                'src/IcTable.h', 'src/Symmetry.cpp', 'src/PDBTopology.cpp', 'src/IcTable.cpp' -Fixing a few memory leaks. IcTable
                 now exposes a deletePointers() method.
0.16.4.3    April 07, 2011    sabs
                'src/DeadEndElimination.h', 'src/DeadEndElimination.cpp' -Fixed GoldSteinSingles after pairs have been flagged.
                 Added runSimpleGoldsteinPairsOnce.
                'src/SelfPairManager.cpp' -Calling runSimpleGoldsteinPairsOnce instead of runSimpleGoldsteinPairs from runOptimizer.
                 Also added code to ensure that DEE is stopped when number of combinations is within the enumerationLimit
0.16.4.2    April 06, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added options to run SimpleGoldsteinPairs and code to runSimpleGoldsteinPairs.
                 Added method to set enumerationLimit and changed default enumerationLimit to 50000.
                'src/Enumerator.h', 'src/Enumerator.cpp' -Initialised valueSet_flag to false. Modified operator[] to return enumeratedValues
                 when values were supplied for states and enumerations are returned otherwise. Added operator () to return states.
                 Added getState and getValue.
                'src/DeadEndElimination.cpp' -Fixed runSimpleGoldSteinPairs.
0.16.4.1    April 01, 2011    sabs
                'src/Scrwl4HBondInteraction.h', 'src/Scrwl4HBondInteraction.cpp' -removing
0.16.4.0    April 01, 2011    sabs
                'Makefile' -Added BaselineEnergyBuilder, BaselineInteraction and changed Scrwl4 to Scwrl4.
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added code to apply a selfEnergyCutoff. Any conformer whose
                 energy is selfECutoff greater than the best conformer for that position is flagged. Any flagged conformer is left
                 out of the selfE and pairE tables, thus saving memory. This functionality is off by default.The IdRotAbsIndex
                 member of SelfPairManager was not being used, so it has been commented out.
                'src/BaselineInteraction.cpp', 'src/BaselineInteraction.h' -This interaction is used to specify a baseline energy
                 for each residue type in a protein. It is a oneBodyInteraction.
                'src/Scwrl4HBondInteraction.cpp', 'src/Scwrl4HBondInteraction.h' -Changed only the name from Scrwl4 to the correct
                 one Scwrl4.
                'src/BaselineEnergyBuilder.h', 'src/BaselineEnergyBuilder.cpp' -Creates BaselineInteractions for the given system.
                
                'src/HydrogenBondBuilder.h', 'src/HydrogenBondBuilder.cpp' -Changed only the name from Scrwl4 to the correct one
                 Scwrl4.
0.16.3.2    March 16, 2011    asenes
                'programs/calculateDistanceOrAngle.cpp' -Now it can take multiple degrees of freedoms at once (print on separate
                 lines). An option printAtoms spits the atom names/numbers after the value
                'Makefile' -
0.16.3.1    March 15, 2011    asenes
                'programs/calculateDistanceOrAngle.cpp' -A quick program that takes a PDB and 2, 3 or 4 atoms (atom ids or atom
                 number) and spits out distance, angle or dihedral
0.16.3.0    March 15, 2011    sabs
                'src/SelfPairManager.h', 'src/SelfPairManager.cpp' -Added getters to get the saved states and energies in SelfPairManager.
                 Deprecated getMCState, replaced with getMCfinalState. SelfPairManager was saving energies per term , disabled
                 this option by default. It can be turned on using the saveEnergiesByTerm method.
0.16.2.0    March 10, 2011    bkmueller
                'src/CoiledCoils.h', 'src/CoiledCoils.cpp' -Renamed primarySequenceToCoiledCoil to setSystemToCoiledCoil (old one
                 present but deprecated), added setBackboneAtomNames to set specific names for the CA C N O atoms, removed comments
                 and unused code, changed internals in setSystemToCoiledCoil
                'programs/coiledCoilBuilder.cpp' -
0.16.1.3    March 07, 2011    jedonald
                'src/PolymerSequence.cpp', 'src/PolymerSequence.h' -Replaced a deprecated function to allow energyTable program
                 to work (from Dan)
                'programs/energyOptimizations.h', 'programs/optimizeMC.cpp' -Updated code to allow optimizeMC to run (from Dan)
                
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
