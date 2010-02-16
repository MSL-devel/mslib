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

#ifndef ALNREADER_H
#define ALNREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

// Storage Formats


// STL Includes
#include <vector>
using namespace std;

/**
 * This class will provide an object which is able
 * to read in and interpret ALN(CLUSTAL-W) files.
 *
 *
 * http://www.ebi.ac.uk/help/formats.html
 *
 * ALN/ClustalW2 format:
 * ALN format was originated in the alignment program ClustalW2. The file starts with word "CLUSTAL" and then some information about which clustal program was run and the version of clustal used
 * e.g. "CLUSTAL W (2.1) multiple sequence alignment"
 * The type of clustal program is "W" and the version is 2.1.
 * The alignment is written in blocks of 60 residues. 
 * Every block starts with the sequence names, obtained from the input sequence, and a count of the total number of residues is shown at the end of the line. 
 * The information about which residues match is shown below each block of residues:
 * "*" means that the residues or nucleotides in that column are identical in all sequences in the alignment.
 * ":" means that conserved substitutions have been observed.
 * "." means that semi-conserved substitutions are observed.
 * 
 * An example is shown below.
 * 
 * 
 * CLUSTAL W 2.1 multiple sequence alignment
 * 
 * 
 * FOSB_MOUSE      ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS 60
 * FOSB_HUMAN      ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPVVDPYDMPGTSYSTPGMSGYSSGGASGS 60
 *                 ********************************.***************:*.**:******
 *
 *
 * CLUSTAL W 2.1 multiple sequence alignment
 * 
 * 
 * seq1            ------------MRVLLAALGLLFLGAL--------------RAFPQDRPFEDTCHGNPS 34
 * seq2            PVAEERGLMSQPLMETCHSVGAAYLESLPLQDASPAGGPSSPRDLPEPRVSTEHTNNKIE 60
 *                             :     ::*  :* :*              * :*: *   :  :.: .
 * 
 * seq1            HYYDKAVRRCCYRCPMGLFPTQQ---CPQRP---TDCRKQCEPDYYLDEADR----CTAC 84
 * seq2            KIYIMKADTVIVGTVKAELPEGRGLAGPAEPELEEELEADHTPHYPEQETEPPLGSCSDV 120
 *                 : *   .         . :*  :    * .*    : . :  *.*  :*::     *:  
 * 
 * seq1            VTCSRDD------------- 91
 * seq2            MLSVEEEGKEDPLPTAASGK 140
 *                 : . .::             
 */
class ALNReader : public Reader {

	public:
		// Constructors/Destructors
		ALNReader();
		ALNReader(const string &_filename);
		ALNReader(const ALNReader & _reader);
		virtual ~ALNReader();

		bool read();

		map<string,string>& getSequences();
		string getSequence(string _key);
	private:
		
		map<string, string> sequences;
		string version;
		string remark;
		
};

inline map<string,string>& ALNReader::getSequences(){
  return sequences;
}
inline string ALNReader::getSequence(string _key){
  map<string,string>::iterator it;
  it = sequences.find(_key);
  if (it == sequences.end()){
    return "";
  } else {
    return it->second;
  }
}

#endif

