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


#include "PSFReader.h"

using namespace MSL;
using namespace std;




bool PSFReader::read(){
	if (!is_open()) {
		return false;
	}



	while (!endOfFileTest()){



		string line = MslTools::trim(Reader::getLine());

		if (line.size() < 3) continue;

		vector<string> toks = MslTools::tokenizeAndTrim(line);

		int exclIndex = -1;
		for (uint i = 0; i < toks.size();i++){
			if (toks[0] == "!"){
				exclIndex = i;
				break;
			}
		}
		if (exclIndex == 0){
			string type = toks[1].substr(1);
			int numberOfRecords = MslTools::toInt(toks[0], "PSFReader::read() , numberOfRecords not an INT on line: "+line);
		
			fprintf(stdout, "Reading Section: %-10s %10d\n",type.c_str(),numberOfRecords);

			if (type == "NATOM"){
				for (uint i = 0 ; i < numberOfRecords;i++){
					vector<string> fields = MslTools::tokenizeAndTrim(MslTools::trim(Reader::getLine()));

					Atom *a = new Atom();
					(*a).setSegID(fields[1]);
					(*a).setResidueNumber(MslTools::toInt(fields[2]));
					(*a).setResidueName(fields[3]);
					(*a).setName(fields[4]);
					(*a).setType(fields[5]);
					(*a).setCharge(MslTools::toDouble(fields[6]));
					(*a).setTempFactor(MslTools::toDouble(fields[7])); // really this is the MASS
					(*a).setElement(fields[8]);
					atoms.push_back(a);

				}

				// Next line/section
				continue;
			}

			if (type == "NBOND"){

				for (uint i = 0; i < numberOfRecords;i++){
					vector<string> fields = MslTools::tokenizeAndTrim(MslTools::trim(Reader::getLine()));

					for (uint j = 0; j < fields.size();j+=2){
						int index1 = MslTools::toInt(fields[j])+1;
						int index2 = MslTools::toInt(fields[j+1])+1;


						atoms(index1).setBoundTo(atoms[index2]);
						atoms(index2).setBoundTo(atoms[index1]);
						
						vector<int> tmp;
						tmp.push_back(index1);
						tmp.push_back(index2);
						bonds.push_back(tmp);
						

						i++;
					}
				}

				// Next line/section
				continue;
			}


			if (type == "NTHETA"){

				for (uint i = 0; i < numberOfRecords;i++){
					vector<string> fields = MslTools::tokenizeAndTrim(MslTools::trim(Reader::getLine()));

					for (uint j = 0; j < fields.size();j+=3){
						int index1 = MslTools::toInt(fields[j])+1;
						int index2 = MslTools::toInt(fields[j+1])+1;						
						int index3 = MslTools::toInt(fields[j+2])+1;						

						vector<int> tmp;
						tmp.push_back(index1);
						tmp.push_back(index2);
						tmp.push_back(index3);
						angles.push_back(tmp);
						

						i++;
					}
				}

				// Next line/section
				continue;
			}

			if (type == "NPHI"){

				for (uint i = 0; i < numberOfRecords;i++){
					vector<string> fields = MslTools::tokenizeAndTrim(MslTools::trim(Reader::getLine()));

					for (uint j = 0; j < fields.size();j+=4){
						
						int index1 = MslTools::toInt(fields[j])+1;
						int index2 = MslTools::toInt(fields[j+1])+1;						
						int index3 = MslTools::toInt(fields[j+2])+1;						
						int index4 = MslTools::toInt(fields[j+3])+1;						
						
						vector<int> tmp;
						tmp.push_back(index1);
						tmp.push_back(index2);
						tmp.push_back(index3);
						tmp.push_back(index4);
						dihedrals.push_back(tmp);
						

						i++;
					}
				}

				// Next line/section
				continue;
			}

			if (type == "NIMPHI"){

				for (uint i = 0; i < numberOfRecords;i++){
					vector<string> fields = MslTools::tokenizeAndTrim(MslTools::trim(Reader::getLine()));


					for (uint j = 0; j < fields.size();j+=4){


						int index1 = MslTools::toInt(fields[j])+1;
						int index2 = MslTools::toInt(fields[j+1])+1;						
						int index3 = MslTools::toInt(fields[j+2])+1;						
						int index4 = MslTools::toInt(fields[j+3])+1;						
						
						vector<int> tmp;
						tmp.push_back(index1);
						tmp.push_back(index2);
						tmp.push_back(index3);
						tmp.push_back(index4);
						impropers.push_back(tmp);
						

						i++;
					}
				}

				// Next line/section
				continue;
			}

			
		}

		

	}

	return true;
}
