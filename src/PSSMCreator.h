/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2011 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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

#ifndef PSSMCREATOR_H
#define PSSMCREATOR_H



// STL Includes


// MSL Includes
#include "System.h"
#include "MslTools.h"
#include "FastaReader.h"

namespace MSL { 
class PSSMCreator {

	public:
		PSSMCreator();
		~PSSMCreator();

		enum PSSMType { logodds=0, freq=1};

		void addMultipleSequenceAlignment(string _fastaFile,string _regex="");
		void setSequences(map<string,string> &_seqs);
		void addReferenceCounts(map<string,double> &_referenceCounts);
		void readReferenceCounts(string _filename);
		void create(PSSMType _type=freq);

		
		vector<double>  getScoreFunction(string _sequence, string _nameRefSeq="", int _beginOffset=0, string _specificScoreType="");
		vector<map<string,double> > getFrequencies(string _nameRefSeq, int _resiBegin, int _resiEnd);
		vector<map<string,double> > getFrequencies();
		string getSequence(string _nameRefSeq, int _resiStart, int _resiEnd);



	private:
		int getInternalNumbering(string _nameRefSeq, int _resi);
		string getSequence(string _name);

		FastaReader fin;
		string fastaFile;
		string referenceFile;
		map<string,string> sequences;
		map<string,double> expectedValues;
		vector<map<string,double> > observedValues;
		vector<map<string,double> > scoreFunction;

};


}

#endif
