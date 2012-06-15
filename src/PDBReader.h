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

#ifndef PDBREADER_H
#define PDBREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"
#include "PDBFormat.h"

// Storage Formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"


// STL Includes
#include <vector>

/**
 * This class will provide an object which is able
 * to read in and interpret PDB files.
 */
namespace MSL { 
class PDBReader : public Reader {

	public:
		// Constructors/Destructors
		PDBReader();
		PDBReader(const std::string &_filename);
		PDBReader(const PDBReader & _reader);
		PDBReader(std::stringstream &_stream);
		~PDBReader();

		// this function assigns coordinates from the atoms of the
		// PDB to an external AtomPointerVector as long as chainId, resnum, resname
		// and atom name are identical.  No errors are assigned for mismatches
		void assignCoordinates(AtomPointerVector & _av);
		struct MissingResidue {
			int model;
			std::string resName;		
			std::string chainId;
			int resNum;
			std::string resIcode;
		};
		struct MissingAtoms {
			int model;
			std::string resName;		
			std::string chainId;
			int resNum;
			std::string resIcode;
			std::vector<std::string> atoms;
		};

		bool read(bool _noHydrogens=false);
		bool read(std::string &_inputString);

		AtomPointerVector & getAtomPointers(); 
		unsigned int size() const;
		Atom * operator[](unsigned int _n);

		std::vector<Matrix *> & getSymmetryRotations();
		std::vector<CartesianPoint *> & getSymmetryTranslations();

		std::vector<Matrix *> & getBioUnitRotations();
		std::vector<CartesianPoint *> & getBioUnitTranslations();

		Matrix & getScaleRotation();
		CartesianPoint & getScaleTranslation();
			
		std::vector<double> & getUnitCellParameters();

		std::map<std::string, double> & getBoundingCoordinates();

		/**
		  Only choose a single alt location for each atom, uses occupancy to decide.
		   -- within a residue can you mix atoms from occupation A and B?
		   ----> I don't think so.  For now once an occupation is choosen for a residue it is always that.
		   ------> Normally they are matched, though (all A are 60%, all B are 40%...)
		 */
		bool getSingleAltLocationFlag();
		void setSingleAltLocationFlag(bool _flag);

		void reset();

		const std::vector<PDBReader::MissingResidue> & getMissingResidues() const;
		const std::vector<PDBReader::MissingAtoms> & getMissingAtoms() const;

		unsigned int getNumberOfModels() const;


		// Operators
		friend PDBReader& operator>>(PDBReader& pdbReader, std::vector<CartesianPoint> &_cv) { return pdbReader;};

	protected:		
	private:
		void deletePointers();
		void parsePDBLine(std::string _pdbline);
		std::string discoverProperResidueAltLoc(std::map<std::string, std::vector<PDBFormat::AtomData> > &_currentResidueAtoms);
		void addAtoms(std::map<std::string, std::vector<PDBFormat::AtomData> > &_currentResidueAtoms, std::string _altLoc);


//		bool allowSloppyPDB;

		AtomPointerVector atoms;

		bool singleAltLocFlag;

		std::vector<Matrix *> symmetryRotations;
		std::vector<CartesianPoint *> symmetryTranslations;

		std::vector<Matrix *> biounitRotations;
		std::vector<CartesianPoint *> biounitTranslations;

		Matrix *scaleRotation;
		CartesianPoint *scaleTranslation;
		
		std::vector<double> unitCellParams;

		std::map<std::string,double> boundingCoords;

		std::vector<MissingResidue> misRes;
		std::vector<MissingAtoms> misAtoms;

		unsigned int numberOfModels;

};

//Inlines go HERE
/**
* This method will return a std::vector of atoms found in this PDB file.
*
* @return A std::vector or atoms from the PDB file.
*/
inline AtomPointerVector& PDBReader::getAtomPointers() { return atoms; }

/**
 * This method will delete all data held in the PDBReader.  All
 * Atom pointers that were previously saved off will no longer be valid
 * after this reset.
 */
inline void PDBReader::reset() {deletePointers();}

inline void PDBReader::assignCoordinates(AtomPointerVector & _av) {
	for (AtomPointerVector::iterator AVatoms=_av.begin(); AVatoms!=_av.end(); AVatoms++) {
		//bool found = false;
		for (AtomPointerVector::iterator PDBatoms=atoms.begin(); PDBatoms!=atoms.end(); PDBatoms++) {
			// in the order of most likely to be different: resnum, resname, atom name, chain id, insertion code
		//		std::cout << "UUU     Check atom " << (*PDBatoms)->getResidueNumber()  << " " << (*PDBatoms)->getResidueName()  << " " << (*PDBatoms)->getName()  << " " << (*PDBatoms)->getChainId()  << " " << (*PDBatoms)->getResidueIcode() << std::endl;
			if ((*PDBatoms)->getResidueNumber() == (*AVatoms)->getResidueNumber() && (*PDBatoms)->getResidueName() == (*AVatoms)->getResidueName() && (*PDBatoms)->getName() == (*AVatoms)->getName() && (*PDBatoms)->getChainId() == (*AVatoms)->getChainId() && (*PDBatoms)->getResidueIcode() == (*AVatoms)->getResidueIcode()) {
				(*AVatoms)->setCoor((*PDBatoms)->getCoor());
				break;
		//		found = true;
			}
		}
		//if (found) {
		//	std::cout << "UUU assigned coor for atom " << (*AVatoms)->getResidueNumber()  << " " << (*AVatoms)->getResidueName()  << " " << (*AVatoms)->getName()  << " " << (*AVatoms)->getChainId()  << " " << (*AVatoms)->getResidueIcode() << std::endl;
		//} else {
		//	std::cout << "UUU NOT assigned coor for atom " << (*AVatoms)->getResidueNumber()  << " " << (*AVatoms)->getResidueName()  << " " << (*AVatoms)->getName()  << " " << (*AVatoms)->getChainId()  << " " << (*AVatoms)->getResidueIcode() << std::endl;
		//}
	}
}


/*
  Why do I need this even though I have one in Reader.h
 */
inline bool PDBReader::read(std::string &_inputString){
	fileName = "string";
	stringStreamPtr = new std::stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}


inline bool PDBReader::getSingleAltLocationFlag() { return singleAltLocFlag; }
inline void PDBReader::setSingleAltLocationFlag(bool _flag){ singleAltLocFlag = _flag;}

inline std::vector<Matrix *> & PDBReader::getSymmetryRotations() { return symmetryRotations;}
inline std::vector<CartesianPoint *> & PDBReader::getSymmetryTranslations() { return symmetryTranslations;}
inline std::vector<Matrix *> & PDBReader::getBioUnitRotations() { return biounitRotations; }
inline std::vector<CartesianPoint *> & PDBReader::getBioUnitTranslations() { return biounitTranslations; }
inline Matrix & PDBReader::getScaleRotation() { return *scaleRotation; }
inline CartesianPoint & PDBReader::getScaleTranslation() { return *scaleTranslation; }
inline std::vector<double> & PDBReader::getUnitCellParameters() { return unitCellParams; }

inline std::map<std::string,double> & PDBReader::getBoundingCoordinates() { return boundingCoords; }
inline const std::vector<PDBReader::MissingResidue> & PDBReader::getMissingResidues() const {return misRes;}
inline const std::vector<PDBReader::MissingAtoms> & PDBReader::getMissingAtoms() const{return misAtoms;}
inline unsigned int PDBReader::getNumberOfModels() const {return numberOfModels;}


}

#endif
