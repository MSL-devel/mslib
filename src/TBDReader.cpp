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


#include "TBDReader.h"
#include "TwoBodyDistanceDependentPotentialTable.h"
#include "MslExceptions.h"

using namespace MSL;
using namespace std;


bool TBDReader::read() { return true;}
bool TBDReader::read(TwoBodyDistanceDependentPotentialTable *_tbd) {
	/*
	  TBD Format:
                    1
	  0123456789012345678901234567890123456789

	  NAME      POTNAME
	  BINS      0-2,2-2.5,2.5-3
	  BINS      8-9,9-10,10-11
	  RESSKIP   7
	  MINDIST   2,100
	  MAXDIST   15,0
	  POTS      ILE_CB,LEU_CB,2,-0.75342

	 */

	cout << "READING FILE: "<<fileName<<endl;

	if (!is_open()) {
		return false;
	}
	
	try{
		while (!endOfFileTest()){

			string line = Reader::getLine();
			if (line.length() < 10) continue;

			string header = line.substr(0,10);
			string value  = line.substr(10);
			vector<string> tokens = MslTools::tokenize(value, ",");


			if (header == "NAME      "){
				cout << "Setting Name: "<<MslTools::trim(value)<<endl;
				_tbd->setPotentialName(MslTools::trim(value));
			}
			if (header == "BINS      "){

				for (uint i = 0 ; i < tokens.size();i++){
					vector<string> range = MslTools::tokenizeAndTrim(tokens[i],"-");
					cout << "Adding bin: "<<range[0]<<" -- "<<range[1]<<endl;
					_tbd->addBin(MslTools::toDouble(range[0], "TBDReader::read(), range[0] not double"),MslTools::toDouble(range[1],"TBDReader::read(), range[1] not double"));
				}
			}



			if (header == "RESSKIP   "){
				cout << "Residue Skip: "<<tokens[0]<<endl;
				_tbd->setResidueSkippingNumber(MslTools::toInt(tokens[0], "TBDReader::read(),resskip -> tokens[0] not int"));
			}

			if (header == "MINDIST   "){
				if (tokens.size() != 2){
					cerr << "ERROR 3253 TBDReader::read(), MINDIST line expecting 2 tokens got: "<<tokens.size()<<" From "<<value<<endl;
				}
				cout << "MINDIST: "<<tokens[0]<<" "<<tokens[1]<<endl;
				_tbd->setMinDistCutoffAndValue(MslTools::toDouble(tokens[0], "TBDReader::read(), mindist -> tokens[0] not double"),MslTools::toDouble(tokens[1], "TBDReader::read(), mindist -> tokens[1] not double"));
			}

			if (header == "MAXDIST   "){
				if (tokens.size() != 2){
					cerr << "ERROR 3253 TBDReader::read(), MAXDIST line expecting 2 tokens got: "<<tokens.size()<<" From "<<value<<endl;
				}
				cout << "MAXDIST: "<<tokens[0]<<" "<<tokens[1]<<endl;
				_tbd->setMaxDistCutoffAndValue(MslTools::toDouble(tokens[0], "TBDReader::read(), maxdist -> tokens[0] not double"),MslTools::toDouble(tokens[1], "TBDReader::read(), maxdist -> tokens[1] not double"));
			}



			if (header == "POTS      "){
				try {
					_tbd->addPotential(tokens[0], tokens[1], MslTools::toInt(tokens[2],"TBDReader::read(), pots ->tokens[2] not int"),MslTools::toDouble(tokens[3],"TBDReader::read(), pots -> tokens[3] not double"));
				} catch (MslSizeException mse){
					cerr << "ERROR TBDReader::read(), asking for a distance bin that does not exist ;MSG; = " << tokens[2]<<"."<<endl; // add msg somehow from mse.
					
				}

			}


		}
	} catch (...){
		cerr << "ERROR 5555 TBDReader::read() got an unexpected exeception, please don't cry. Hey atleast I told you what function it happened in.."<<endl;
		exit(5555);
	}


	return true;
}
