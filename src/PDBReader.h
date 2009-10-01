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
#include "AtomVector.h"


// STL Includes
#include <vector>
using namespace std;

/**
 * This class will provide an object which is able
 * to read in and interpret PDB files.
 */
class PDBReader : public Reader {

	public:
		// Constructors/Destructors
		PDBReader();
		PDBReader(const string &_filename);
		PDBReader(const PDBReader & _reader);
		PDBReader(stringstream &_stream);
		virtual ~PDBReader();

		// this function assigns coordinates from the atoms of the
		// PDB to an external AtomVector as long as chainId, resnum, resname
		// and atom name are identical.  No errors are assigned for mismatches
		void assignCoordinates(AtomVector & _av);
		struct MissingResidue {
			int model;
			string resName;		
			string chainId;
			int resNum;
			string resIcode;
		};
		struct MissingAtoms {
			int model;
			string resName;		
			string chainId;
			int resNum;
			string resIcode;
			vector<string> atoms;
		};

		bool read();
		bool read(string &_inputString);

		AtomVector & getAtoms(); 
		unsigned int size() const;
		Atom * operator[](unsigned int _n);

		vector<Matrix *> & getSymmetryRotations();
		vector<CartesianPoint *> & getSymmetryTranslations();

		vector<Matrix *> & getBioUnitRotations();
		vector<CartesianPoint *> & getBioUnitTranslations();

		Matrix & getScaleRotation();
		CartesianPoint & getScaleTranslation();
			
		vector<double> & getUnitCellParameters();

		map<string, double> & getBoundingCoordinates();

		/**
		  Only choose a single alt location for each atom, uses occupancy to decide.
		   -- within a residue can you mix atoms from occupation A and B?
		   ----> I don't think so.  For now once an occupation is choosen for a residue it is always that.
		   ------> Normally they are matched, though (all A are 60%, all B are 40%...)
		 */
		bool getSingleAltLocationFlag();
		void setSingleAltLocationFlag(bool _flag);

		void reset();

		const vector<PDBReader::MissingResidue> & getMissingResidues() const;
		const vector<PDBReader::MissingAtoms> & getMissingAtoms() const;

		// Operators
		friend PDBReader& operator>>(PDBReader& pdbReader, vector<CartesianPoint> &_cv) { return pdbReader;};

	protected:		
	private:
		void deletePointers();
		void parsePDBLine(string _pdbline);
		string discoverProperResidueAltLoc(map<string, vector<PDBFormat::AtomData> > &_currentResidueAtoms);
		void addAtoms(map<string, vector<PDBFormat::AtomData> > &_currentResidueAtoms, string _altLoc);



		AtomVector atoms;

		bool singleAltLocFlag;

		vector<Matrix *> symmetryRotations;
		vector<CartesianPoint *> symmetryTranslations;

		vector<Matrix *> biounitRotations;
		vector<CartesianPoint *> biounitTranslations;

		Matrix *scaleRotation;
		CartesianPoint *scaleTranslation;
		
		vector<double> unitCellParams;

		map<string,double> boundingCoords;

		vector<MissingResidue> misRes;
		vector<MissingAtoms> misAtoms;

};

//Inlines go HERE
/**
 * Simple constructor.
 */
inline PDBReader::PDBReader() : Reader() { singleAltLocFlag = false; scaleTranslation = new CartesianPoint(0.0,0.0,0.0); scaleRotation = new Matrix(3,3,0.0);(*scaleRotation)[0][0] = 1.0;(*scaleRotation)[1][1] = 1.0;(*scaleRotation)[2][2] = 1.0;}
/**
 * With this constructor the user specifies the filename
 * of the PDB to be read.
 *
 * @param _filename  The name of the PDB file to be read.
 */
inline PDBReader::PDBReader(const string &_filename) : Reader(_filename) { singleAltLocFlag = false; scaleTranslation = new CartesianPoint(0.0,0.0,0.0); scaleRotation = new Matrix(3,3,0.0);(*scaleRotation)[0][0] = 1.0;(*scaleRotation)[1][1] = 1.0;(*scaleRotation)[2][2] = 1.0;}
/**
 * A copy constructor.  All of the atoms from the given PDBReader are
 * copied into the new PDBReader.
 *
 * @param _reader The PDBReader to be copied.
 */
inline PDBReader::PDBReader(const PDBReader & _reader) {
	for (AtomVector::const_iterator k=_reader.atoms.begin(); k!= _reader.atoms.end(); k++) {
		atoms.push_back(new Atom(**k));
	}
	singleAltLocFlag = _reader.singleAltLocFlag;
	scaleTranslation = new CartesianPoint(0.0,0.0,0.0); 
	scaleRotation = new Matrix(3,3,0.0);
}
/**
 * A constructor which will read input data from a stringstream.
 *
 * @param _ss The stringstream to get data from.
 */
inline PDBReader::PDBReader(stringstream &_ss) : Reader(_ss)     {read();}

/**
 * The deconstructor.  All data will be deleted, so any Atom pointers
 * that were previously saved off will no longer be valid after the PDBReader
 * object has been destroyed.
 */
inline PDBReader::~PDBReader() { deletePointers(); close();}

/**
* This method will return a vector of atoms found in this PDB file.
*
* @return A vector or atoms from the PDB file.
*/
inline AtomVector& PDBReader::getAtoms() { return atoms; }

/**
 * This method will delete all data held in the PDBReader.  All
 * Atom pointers that were previously saved off will no longer be valid
 * after this reset.
 */
inline void PDBReader::reset() {deletePointers();}

inline void PDBReader::assignCoordinates(AtomVector & _av) {
	for (AtomVector::iterator AVatoms=_av.begin(); AVatoms!=_av.end(); AVatoms++) {
		//bool found = false;
		for (AtomVector::iterator PDBatoms=atoms.begin(); PDBatoms!=atoms.end(); PDBatoms++) {
			// in the order of most likely to be different: resnum, resname, atom name, chain id, insertion code
		//		cout << "UUU     Check atom " << (*PDBatoms)->getResidueNumber()  << " " << (*PDBatoms)->getResidueName()  << " " << (*PDBatoms)->getName()  << " " << (*PDBatoms)->getChainId()  << " " << (*PDBatoms)->getResidueIcode() << endl;
			if ((*PDBatoms)->getResidueNumber() == (*AVatoms)->getResidueNumber() && (*PDBatoms)->getResidueName() == (*AVatoms)->getResidueName() && (*PDBatoms)->getName() == (*AVatoms)->getName() && (*PDBatoms)->getChainId() == (*AVatoms)->getChainId() && (*PDBatoms)->getResidueIcode() == (*AVatoms)->getResidueIcode()) {
				(*AVatoms)->setCoor((*PDBatoms)->getCoor());
				break;
		//		found = true;
			}
		}
		//if (found) {
		//	cout << "UUU assigned coor for atom " << (*AVatoms)->getResidueNumber()  << " " << (*AVatoms)->getResidueName()  << " " << (*AVatoms)->getName()  << " " << (*AVatoms)->getChainId()  << " " << (*AVatoms)->getResidueIcode() << endl;
		//} else {
		//	cout << "UUU NOT assigned coor for atom " << (*AVatoms)->getResidueNumber()  << " " << (*AVatoms)->getResidueName()  << " " << (*AVatoms)->getName()  << " " << (*AVatoms)->getChainId()  << " " << (*AVatoms)->getResidueIcode() << endl;
		//}
	}
}


/*
  Why do I need this even though I have one in Reader.h
 */
inline bool PDBReader::read(string &_inputString){
	fileName = "string";
	stringStreamPtr = new stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}


inline bool PDBReader::getSingleAltLocationFlag() { return singleAltLocFlag; }
inline void PDBReader::setSingleAltLocationFlag(bool _flag){ singleAltLocFlag = _flag;}

inline vector<Matrix *> & PDBReader::getSymmetryRotations() { return symmetryRotations;}
inline vector<CartesianPoint *> & PDBReader::getSymmetryTranslations() { return symmetryTranslations;}
inline vector<Matrix *> & PDBReader::getBioUnitRotations() { return biounitRotations; }
inline vector<CartesianPoint *> & PDBReader::getBioUnitTranslations() { return biounitTranslations; }
inline Matrix & PDBReader::getScaleRotation() { return *scaleRotation; }
inline CartesianPoint & PDBReader::getScaleTranslation() { return *scaleTranslation; }
inline vector<double> & PDBReader::getUnitCellParameters() { return unitCellParams; }

inline map<string,double> & PDBReader::getBoundingCoordinates() { return boundingCoords; }
inline const vector<PDBReader::MissingResidue> & PDBReader::getMissingResidues() const {return misRes;}
inline const vector<PDBReader::MissingAtoms> & PDBReader::getMissingAtoms() const{return misAtoms;}

#endif
