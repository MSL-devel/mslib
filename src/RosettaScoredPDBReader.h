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

#ifndef ROSETTASCOREDPDBREADER_H
#define ROSETTASCOREDPDBREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

// STL Includes
#include <vector>
#include <map>
using namespace std;

/**
 * This class will provide an object which is able
 * to read in and interpret PDB files.
 */
namespace MSL { 
class RosettaScoredPDBReader : public Reader {

	public:
		// Constructors/Destructors
		RosettaScoredPDBReader();
		RosettaScoredPDBReader(const std::string &_filename);
		RosettaScoredPDBReader(const RosettaScoredPDBReader & _reader);
		~RosettaScoredPDBReader() {}

		bool read();
		bool read(std::string &_inputString);

		// CHAIN_RESNUM_ICODE, TYPE_OF_SCORE = SCORE_VALUE
		map<string,map<string,double> > & getResidueScores();
		
	protected:		
	private:
		vector<string> score_labels;
		map<string, map<string,double> > residueScores;
		vector<pair<string,string> > aaSeq; // Store type of AA + positionId

};

//Inlines go HERE

  inline map<string, map<string,double> > & RosettaScoredPDBReader::getResidueScores() { return residueScores;}



/*
  Why do I need this even though I have one in Reader.h
 */
inline bool RosettaScoredPDBReader::read(std::string &_inputString){
	fileName = "string";
	stringStreamPtr = new std::stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}

}

#endif
