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

#define MSLVERSION "0.0.2.2"
#define MSLDATE "October 08, 2009"

/*
HISTORY:
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
