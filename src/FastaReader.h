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

#ifndef FASTAREADER_H
#define FASTAREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"
#include "Chain.h"

// Storage Formats


// STL Includes
#include <vector>

/**
 * This class will provide an object which is able
 * to read in and interpret FASTA files.
 *
 *
 * >FOSB_MOUSE
 * ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS 60
 * >FOSB_HUMAN
 * ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPVVDPYDMPGTSYSTPGMSGYSSGGASGS 60
 *
 */
namespace MSL { 
class FastaReader : public Reader {

	public:
		// Constructors/Destructors
		FastaReader();
		FastaReader(const std::string &_filename);
		FastaReader(const FastaReader & _reader);
		virtual ~FastaReader();

		bool read(std::string _regex="");

		std::map<std::string,std::string>& getSequences();
		std::string getSequence(std::string _key);
		int getIndex(Chain &_ch, std::string _posId, std::string _key);
		std::string getPositionId(Chain &_ch, int _index, std::string _key);
	private:
		
		std::map<std::string, std::string> sequences;
		
};

inline std::map<std::string,std::string>& FastaReader::getSequences(){
  return sequences;
}
inline std::string FastaReader::getSequence(std::string _key){
  std::map<std::string,std::string>::iterator it;
  it = sequences.find(_key);
  if (it == sequences.end()){
    return "";
  } else {
    return it->second;
  }
}

}

#endif

