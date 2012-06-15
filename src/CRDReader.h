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

#ifndef CRDREADER_H
#define CRDREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"
#include "CRDFormat.h"
#include "CharmmTopologyReader.h"

// Storage Formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"


// STL Includes
#include <vector>

/**
 * This class will provide an object which is able
 * to read in and interpret CRD files.
 */
namespace MSL { 
class CRDReader : public Reader {

	public:
		// Constructors/Destructors
		CRDReader();
		CRDReader(const std::string &_filename);
		CRDReader(const std::string &_filename, std::string _topologyFile);
		CRDReader(const CRDReader & _reader);
		CRDReader(std::stringstream &_stream);
		CRDReader(std::stringstream &_stream,std::string _topologyFile);
		virtual ~CRDReader();

		// this function assigns coordinates from the atoms of the
		// CRD to an external AtomPointerVector as long as chainId, resnum, resname
		// and atom name are identical.  No errors are assigned for mismatches
		void assignCoordinates(AtomPointerVector & _av);
		
		bool read();
		bool read(std::string &_inputString);

		AtomPointerVector & getAtomPointers(); 
		unsigned int size() const;
		Atom * operator[](unsigned int _n);

		bool readTopology(std::string _topologyFile);

		void reset();

		// Operators
		friend CRDReader& operator>>(CRDReader& pdbReader, std::vector<CartesianPoint> &_cv) { return pdbReader;};

	protected:		
	private:
		void setup(std::string _topologyFile);
		void deletePointers();
		void addAtoms(std::map<std::string, std::vector<CRDFormat::AtomData> > &_currentResidueAtoms, std::string _altLoc);

		AtomPointerVector atoms;
		CharmmTopologyReader *pTopReader;

};

//Inlines go HERE
/**
 * Simple constructor.
 */
inline CRDReader::CRDReader() : Reader() { setup(""); }
/**
 * With this constructor the user specifies the filename
 * of the CRD to be read.
 *
 * @param _filename  The name of the CRD file to be read.
 */
inline CRDReader::CRDReader(const std::string &_filename) : Reader(_filename) { setup("");}
inline CRDReader::CRDReader(const std::string &_filename,std::string _topologyFile) : Reader(_filename) { setup(_topologyFile); }
/**
 * A copy constructor.  All of the atoms from the given CRDReader are
 * copied into the new CRDReader.
 *
 * @param _reader The CRDReader to be copied.
 */
inline CRDReader::CRDReader(const CRDReader & _reader) {
	deletePointers();
	for (AtomPointerVector::const_iterator k=_reader.atoms.begin(); k!= _reader.atoms.end(); k++) {
		atoms.push_back(new Atom(**k));
	}
	if (_reader.pTopReader != NULL) {
		pTopReader = new CharmmTopologyReader(*(_reader.pTopReader));
	}
}
/**
 * A constructor which will read input data from a std::stringstream.
 *
 * @param _ss The std::stringstream to get data from.
 */
inline CRDReader::CRDReader(std::stringstream &_ss) : Reader(_ss)     {
	setup("");
	read();
}
inline CRDReader::CRDReader(std::stringstream &_ss, std::string _topologyFile) : Reader(_ss)     {
	setup(_topologyFile);
	read();
}

inline void CRDReader::setup(std::string _topologyFile) {
	pTopReader = NULL;
	if(_topologyFile != "") {
		readTopology(_topologyFile);
	}
}

inline bool CRDReader::readTopology(std::string _topologyFile) {
	if(pTopReader) {
		pTopReader->reset();
	} else {
		pTopReader = new CharmmTopologyReader;
	}

	if (!pTopReader->open(_topologyFile)) {
		return false;
	}
	bool out = false;
	if (pTopReader->read()) {
		out = true;
	}
	pTopReader->close();
	return out;
}

/**
 * The deconstructor.  All data will be deleted, so any Atom pointers
 * that were previously saved off will no longer be valid after the CRDReader
 * object has been destroyed.
 */
inline CRDReader::~CRDReader() { deletePointers(); close();}

/**
* This method will return a std::vector of atoms found in this CRD file.
*
* @return A std::vector or atoms from the CRD file.
*/
inline AtomPointerVector& CRDReader::getAtomPointers() { return atoms; }

/**
 * This method will delete all data held in the CRDReader.  All
 * Atom pointers that were previously saved off will no longer be valid
 * after this reset.
 */
inline void CRDReader::reset() {deletePointers();}

inline void CRDReader::assignCoordinates(AtomPointerVector & _av) {
	for (AtomPointerVector::iterator AVatoms=_av.begin(); AVatoms!=_av.end(); AVatoms++) {
		//bool found = false;
		for (AtomPointerVector::iterator CRDatoms=atoms.begin(); CRDatoms!=atoms.end(); CRDatoms++) {
			// in the order of most likely to be different: resnum, resname, atom name, chain id, insertion code
		//		std::cout << "UUU     Check atom " << (*CRDatoms)->getResidueNumber()  << " " << (*CRDatoms)->getResidueName()  << " " << (*CRDatoms)->getName()  << " " << (*CRDatoms)->getChainId()  << " " << (*CRDatoms)->getResidueIcode() << std::endl;
			if ((*CRDatoms)->getResidueNumber() == (*AVatoms)->getResidueNumber() && (*CRDatoms)->getResidueName() == (*AVatoms)->getResidueName() && (*CRDatoms)->getName() == (*AVatoms)->getName() && (*CRDatoms)->getChainId() == (*AVatoms)->getChainId() && (*CRDatoms)->getResidueIcode() == (*AVatoms)->getResidueIcode()) {
				(*AVatoms)->setCoor((*CRDatoms)->getCoor());
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



inline bool CRDReader::read(std::string &_inputString){
	fileName = "string";
	stringStreamPtr = new std::stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}

 
}


#endif

