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

#include "LinearProgrammingOptimization.h"

#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "MslTools.h"

using namespace MslTools;


LinearProgrammingOptimization::LinearProgrammingOptimization(){

	// Set defaults
	selfEnergy              = NULL;
	pairEnergy              = NULL;
	totalNumRotamers        = 0;
	totalNumPositions       = 0;
	numConstraints          = 0;
	numDecisionVariables    = 0;
	responsibleForEnergyTableMemory = false;
	verbose = false;
	
}
LinearProgrammingOptimization::~LinearProgrammingOptimization(){


	
}

void LinearProgrammingOptimization::readEnergyTable(string _filename){
	/*
	  ENERGY TABLE FILE FORMAT:
	  
	  # COMMENT LINES
	  # SELF TERMS
	  #  POSITION_INDEX ROTAMER_INDEX SELF_ENERGY
	  0 0 -0.5
	  0 1 -0.3
	  1 0  0.1
	  ..

	  # PAIR TERMS
	  #  POSITION_INDEX ROTAMER_INDEX POSITION_INDEX ROTAMER_INDEX PAIR_ENERGY
	  0 0 1 0 -0.2
	  0 0 1 1  0.5
	  ...
	  
	*/

	// This object is now responsible for the energy table memory.
	responsibleForEnergyTableMemory = true;

	selfEnergy = new vector<vector<double> >();
	pairEnergy = new vector<vector<vector<vector<double> > > >();

	ifstream fin;
	string line;
	try {
		fin.open(_filename.c_str());
		if (fin.fail()){
			cerr << "ERROR 8904 in LinearProgrammingOptimization::readEnergyTable opening file "<<_filename<<endl;
			exit(8904);
		}
		
		bool firstPair = true;
		while (!fin.fail()){
			getline(fin,line);
			if (line[0] == '#') continue;
			if (line.length() <= 1) continue;



			vector<string> toks = tokenize(line);

			
			// self Energy Line has 3 numbers on it
			if (toks.size() == 3){
				int posIndex = toInt(toks[0]);
				if (selfEnergy->size() < posIndex+1){
					vector<double> tmp;
					selfEnergy->push_back(tmp);
				}

				(*selfEnergy)[posIndex].push_back(toDouble(toks[2]));
				
			}



			// pair Energy Line has 5 numbers on it
			/*
			  TODO:
			     Make sure pairEnergy is not a full sized table.
			      it should only be a triangle - lower triangular.
			      resize(i+1) for pairEnergy[i]?
			 */
			if (toks.size() == 5){

				if (firstPair){
					pairEnergy->resize(selfEnergy->size());
					// For each position resize to
					// num rotamers
					for (uint i = 0; i < selfEnergy->size();i++){
						(*pairEnergy)[i].resize((*selfEnergy)[i].size());

						// For each rotamer,
						// resize to num positions
						for (uint j = 0; j < (*selfEnergy)[i].size();j++){
							(*pairEnergy)[i][j].resize((*selfEnergy).size());

							// For each position resize to num rotamers
							for (uint k = 0; k < (*selfEnergy).size();k++){
								(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

							}
						}

					}
				}
				firstPair = false;
				
				// Add to pairEnergy object. Indices pairEnergy[POS1][ROT1][POS2][ROT2] = ENERGY;
				(*pairEnergy)[toInt(toks[0])][toInt(toks[1])][toInt(toks[2])][toInt(toks[3])] = toDouble(toks[4]);

				// Add symmetric entries into the pairEnergy table ****NO NEED****
				// pairEnergy[toInt(toks[2])][toInt(toks[3])][toInt(toks[0])][toInt(toks[1])] = toDouble(toks[4]);
				
			}
		}
	} catch (...){
		cerr << "ERROR 7809 in LinearProgrammingOptimization::readEnergyTable in file: "<<_filename<<endl;
		exit(7809);
	}

}


void LinearProgrammingOptimization::analyzeEnergyTable(){
	cout << "Done reading energy file, now anaylize properites... "<<endl;


	// PAIRTYPE:
	// Define each pair of positions as POS or NEG.
	//   POS = at least one rotamer pair between two positions is positive.
	//   NEG = not POS

	// rotamerState,rotamerPairs
	// Define an ordered list of rotamer indexes:  
	//   1. self-rotamer variable (rotamer state)
	//   2. pair-rotamer variable (rotamer pair)
	

	// Elimination of rotamers...
	//  1. All rotamers between 2 positions is 0 eliminate all constraints between those two.
	//  2. Input Mask?

	// Keep track of number of constraints
	numConstraints = 0;

	// PairType 
	pairType.resize(pairEnergy->size());

	fprintf(stdout, "\tSelfEnergy   : %10d\n",(int)selfEnergy->size());
	fprintf(stdout, "\tPairEnergy   : %10d\n",(int)pairEnergy->size());

	// For each pair of rotamers..
	for (uint pos1 = 0; pos1 < pairEnergy->size();pos1++){

		// Resize pairType[pos1]
		pairType[pos1].resize(pairEnergy->size());


		cout << "PairE_I:"<<pos1<<" "<<(*pairEnergy)[pos1].size()<<endl;
		for (uint pos2 = 0; pos2 < pairEnergy->size();pos2++){

			cout << "PairE_J:"<<pos2<<" "<<(*pairEnergy)[pos2].size()<<endl;


			// Initialize positive to false, meaning we 
			bool positiveFlag = false;
			for (uint rot1 = 0; rot1 < (*pairEnergy)[pos1].size();rot1++){			

				// addToMask, skip here? numConstraints, rotamerState, rotamerPair
				if (inputMasks.size() > pos1 && inputMasks[pos1].size() > rot1 && !inputMasks[pos1][rot1]){
					continue;
				}
				
				if (pos1 != pos2){
					numConstraints++;
				}
				// Keep track of each self-rotamer only the first time over pos2.
				if (pos2 == 0){


					vector<int> tmp;	
					tmp.push_back(pos1);
					tmp.push_back(rot1);
					rotamerState.push_back(tmp);

				}
				for (uint rot2 = 0; rot2 < (*pairEnergy)[pos2].size();rot2++){			
					

					//cout << "PairENERGY: "<<pos1<<","<<rot1<<","<<pos2<<","<<rot2<<endl;
					if (pos1 < pos2 && (*pairEnergy)[pos1][rot1][pos2][rot2] > 0){
						positiveFlag = true;
					}

					if (pos1 > pos2 && (*pairEnergy)[pos2][rot2][pos1][rot1] > 0){
						positiveFlag = true;
					}

					if (pos1 == pos2) continue;

					if (pos1 < pos2){
						vector<int> tmp;
						tmp.push_back(pos1);
						tmp.push_back(rot1);
						tmp.push_back(pos2);
						tmp.push_back(rot2);
						rotamerPair.push_back(tmp);
					}
				}

			}

			pairType[pos1][pos2] = positiveFlag;
			if (verbose) {
				cout << "PairType: "<<pos1<<","<<pos2<<" = "<<positiveFlag<<endl;
			}

		}
	}



	// Store total num positions
	totalNumPositions = selfEnergy->size();

	// Compute total number of rotatmers, store 'self' rotamer pairs in columnTable
	totalNumRotamers = 0;
	for (uint i = 0; i < totalNumPositions;i++){
		totalNumRotamers += (*selfEnergy)[i].size();
	}

	numConstraints += totalNumPositions;

	numDecisionVariables = rotamerState.size() + rotamerPair.size();

	cout << "Totals: "<<endl;
	fprintf(stdout, "\tNumPositions : %10d\n",totalNumPositions);
	fprintf(stdout, "\tNumRotamers  : %10d\n",totalNumRotamers);
	fprintf(stdout, "\tRotamerStates: %10d\n",(int)rotamerState.size());
	fprintf(stdout, "\tRotamerPairs : %10d\n",(int)rotamerPair.size());
	fprintf(stdout, "\tPairType     : %10d\n",(int)pairType.size());
	fprintf(stdout, "\tSelfEnergy   : %10d\n",(int)selfEnergy->size());
	fprintf(stdout, "\tPairEnergy   : %10d\n",(int)pairEnergy->size());
	fprintf(stdout, "\tConstraints  : %10d\n",numConstraints);
	fprintf(stdout, "\tDecisions    : %10d\n",numDecisionVariables);

}
void LinearProgrammingOptimization::solveLP(){


	glp_simplex(lp,NULL);
	glp_intopt(lp,NULL);
	double Z = glp_get_obj_val(lp);
	double Zint = glp_mip_obj_val(lp);


	masks.clear();
	masks.resize(totalNumPositions);
	rotamerSelection.clear();
	rotamerSelection.resize(totalNumPositions);
	for (uint i = 0; i < totalNumPositions;i++){
		masks[i].resize(rotamerState[i].size());
	}
	for (uint i = 1; i <= rotamerState.size()+rotamerPair.size();i++){
		double x = glp_get_col_prim(lp, i);
		double xint = glp_mip_col_val(lp, i);
		string s = glp_get_col_name(lp, i);
		int kind = glp_get_col_kind(lp, i);

		printf("%10s,%4d  x[%4d] = %8.3f,%8.3f; ", s.c_str(), kind, i,x,xint);

		// 
		if (i <= rotamerState.size()){
			masks[rotamerState[i-1][0]][rotamerState[i-1][1]] = xint;
			if (xint){
				rotamerSelection[rotamerState[i-1][0]] = rotamerState[i-1][1];
			}
			printf(" RS");
		}
		printf("\n");

		// Break out of loop if verbose flag is not set
		if (!verbose && i >= rotamerState.size()) break;

		
	}
	printf("\nZ = %8.3f; %8.3f \n", Z,Zint);
	glp_delete_prob(lp);
	

}



void LinearProgrammingOptimization::createLP(){

	if (selfEnergy->size() == 0 || pairEnergy->size() == 0){
		cerr << "ERROR 2436 LinearProgrammingOptimization::createLPproblem selfEnergy or pairEnergy lists are NULL or 0 in size: "<<selfEnergy->size()<<" "<<pairEnergy->size()<<endl;
		exit(2436);
	}
	


	// Create an LP problem object
	lp = glp_create_prob();

	glp_set_prob_name(lp, "rotamerOptimization");

	// Set direction to MINIMIZE (not MAXIMIZE)
	glp_set_obj_dir(lp,GLP_MIN);



	/*
	      Add Constraints/ROWS
	*/

	cout << "Add constraints/rows\n";

	glp_add_rows(lp, numConstraints);
	map<int,vector<int> > indexConstraints;

	// Set up positional constraints
	cout << "\tPositional Constriants\n";
	stringstream consName;
	int index = 1;
	for (uint i = 0; i < totalNumPositions; i++){

		vector<int> tmp;
		tmp.push_back(i);
		indexConstraints[index-1] = tmp;

		consName.str("");
 		consName << "Pos-"<<i;
		glp_set_row_name(lp, index, consName.str().c_str());
		glp_set_row_bnds(lp, index, GLP_FX, 1, 1);

		if (verbose){
			cout << "\t\t "<<consName.str()<<" ("<<index<<")"<<endl;
		}
		index++;
	}


	// Set up postionX_positionY+rotamerZ constraints	
	cout << "\tPos-PosRot Constriants\n";
	for (uint pos1 = 0; pos1 < pairEnergy->size();pos1++){
		for (uint pos2 = 0; pos2 < pairEnergy->size();pos2++){
			if (pos1 == pos2) continue;

			for (uint rot2 = 0; rot2 < (*pairEnergy)[pos2].size();rot2++){

				vector<int> tmp;
				tmp.push_back(pos1);
				tmp.push_back(pos2);
				tmp.push_back(rot2);
				indexConstraints[index-1] = tmp;

				consName.str("");
				consName << "Pair-"<<pos1<<"_"<<pos2<<"-"<<rot2;
				glp_set_row_name(lp,index,consName.str().c_str());

				bool result = pairType[pos1][pos2];
				if (pos2 > pos1){
					result = pairType[pos2][pos1];
				}
				if (result){
					glp_set_row_bnds(lp,index,GLP_FX,0,0);
					if (verbose){
						cout <<"\t\t"<<consName.str()<<" ("<<index<<") BOUNDED TO 0."<<endl;
					}
				} else {
					glp_set_row_bnds(lp,index,GLP_UP,0,0);

					if (verbose){
						cout <<"\t\t"<<consName.str()<<" ("<<index<<") BOUNDED UPTO 0."<<endl;
					}
				}

				index++;


			}
		}
	}

	
	/*
	     Add decision variables (one per rotamer pair).
	          rotamer state 
		  rotamer pair
	*/

	cout << "Add decision variables\n";
	stringstream colName;
	glp_add_cols(lp, rotamerState.size()+rotamerPair.size());


	// Add rotamer selected
	index = 1;
	cout << "\tRotamer States\n";
	for (uint i = 0 ;i < rotamerState.size();i++){

		colName.str("");
		colName << "RotState-"<<rotamerState[i][0]<<"-"<<rotamerState[i][1];
		glp_set_col_name(lp,index,colName.str().c_str());
		glp_set_col_kind(lp,index,GLP_IV);
		glp_set_col_bnds(lp,index,GLP_DB,0,1);
		glp_set_obj_coef(lp,index++,(*selfEnergy)[rotamerState[i][0]][rotamerState[i][1]]);

		if (verbose){
			cout << "\t\tState "<<i<<"("<<index-1<<") : "<<(*selfEnergy)[rotamerState[i][0]][rotamerState[i][1]]<<endl;
		}
	}


	// Add pair/edge selected
	cout << "\tRotamer Pairs\n";
	for (uint i = 0; i < rotamerPair.size();i++){

		colName.str("");
		colName << "RotPair-"<<rotamerPair[i][0]<<"-"<<rotamerPair[i][1]<<"_"<<rotamerPair[i][2]<<"-"<<rotamerPair[i][3];
		glp_set_col_name(lp,index,colName.str().c_str());
		glp_set_col_kind(lp,index,GLP_IV);
		glp_set_col_bnds(lp,index,GLP_DB,0,1);


		if (rotamerPair[i][0] < rotamerPair[i][2]) {
			glp_set_obj_coef(lp,index,(*pairEnergy)[rotamerPair[i][0]][rotamerPair[i][1]][rotamerPair[i][2]][rotamerPair[i][3]]);

			if (verbose){
				cout << "\t\t"<<colName.str()<<" "<<i<<"("<<index<<") : "<<(*pairEnergy)[rotamerPair[i][0]][rotamerPair[i][1]][rotamerPair[i][2]][rotamerPair[i][3]]<<endl;
			}
		} else {
			glp_set_obj_coef(lp,index,(*pairEnergy)[rotamerPair[i][2]][rotamerPair[i][3]][rotamerPair[i][0]][rotamerPair[i][1]]);

			if (verbose){
				cout << "\t\t"<<colName.str()<<" "<<i<<"("<<index<<") : "<<(*pairEnergy)[rotamerPair[i][2]][rotamerPair[i][3]][rotamerPair[i][0]][rotamerPair[i][1]]<<endl;
			}
		}

		index++;

	}

	//cout << "Num Integer Columns: "<<glp_get_num_int(lp)<<endl;


	// For each constraint/row
	cout << "Compute non-zero coefficent matrix \n";
	int nonZeroCount = 0;
	vector<vector<double> > coeff;
	for (uint i = 0; i < numConstraints;i++){
		for (uint j = 0; j < numDecisionVariables;j++){
			if (i < totalNumPositions){

				if (j < rotamerState.size() && rotamerState[j][0] == indexConstraints[i][0]){
					nonZeroCount++;
					vector<double> tmp;
					tmp.push_back(i+1);
					tmp.push_back(j+1);
					tmp.push_back(1);
					coeff.push_back(tmp);
					
				}
			} else {
				if (j >= rotamerState.size() && rotamerPair[j-rotamerState.size()][0] == indexConstraints[i][0] &&
				                                rotamerPair[j-rotamerState.size()][2] == indexConstraints[i][1] &&
				                                rotamerPair[j-rotamerState.size()][3] == indexConstraints[i][2]){
					vector<double> tmp;
					tmp.push_back(i+1);
					tmp.push_back(j+1);
					tmp.push_back(1);
					coeff.push_back(tmp);
					nonZeroCount++;
				} else {

					if (j >= rotamerState.size() && rotamerPair[j-rotamerState.size()][2] == indexConstraints[i][0] &&
					    rotamerPair[j-rotamerState.size()][0] == indexConstraints[i][1] &&
					    rotamerPair[j-rotamerState.size()][1] == indexConstraints[i][2]){
						vector<double> tmp;
						tmp.push_back(i+1);
						tmp.push_back(j+1);
						tmp.push_back(1);
						coeff.push_back(tmp);
						nonZeroCount++;
					}
				}

				// Inside rotamerState columns and column Position2 == row Position2 and column Rotamer2 == row Rotamer 2.
				if (j < rotamerState.size() && rotamerState[j][0] == indexConstraints[i][1] && rotamerState[j][1] == indexConstraints[i][2]){
					vector<double> tmp;
					tmp.push_back(i+1);
					tmp.push_back(j+1);
					tmp.push_back(-1);
					coeff.push_back(tmp);
					nonZeroCount++;
				}
			}
		}
	}

	cout << "Create coefficent matrix ("<<(nonZeroCount+1)<<")\n";
	int row[(nonZeroCount+1)];
	int col[(nonZeroCount+1)];
	double val[(nonZeroCount+1)];

	index = 1;

	for (uint i = 1; i <= nonZeroCount;i++){
		row[i] = (int)coeff[i-1][0];
		col[i] = (int)coeff[i-1][1];
		val[i] = coeff[i-1][2];
	}




	cout << "Load Matrix"<<endl;
	// Load coefficient matrix from constraints
	glp_load_matrix(lp, nonZeroCount, row,col,val);

}


void LinearProgrammingOptimization::printMe(bool _selfOnly){

	fprintf(stdout,"Self terms:\n");
	for (uint i = 0; i < (*selfEnergy).size();i++){
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){

			fprintf(stdout, "    %4d %4d %8.3f", i, j, (*selfEnergy)[i][j]);
			if (masks[i][j]) {
				fprintf(stdout, " **** ");
			}
			fprintf(stdout,"\n");


		}

	}
	if (_selfOnly) return;

	fprintf(stdout,"Pair terms:\n");
	for (uint i = 0; i < (*pairEnergy).size();i++){
		for (uint j = 0; j < (*pairEnergy)[i].size();j++){
			//for (uint k = i+1 ; k < (*pairEnergy)[i][j].size();k++){	
			for (uint k = 0 ; k < (*pairEnergy)[i][j].size();k++){	
				if (i == k){
					continue;
				}
				for (uint l = 0 ; l < (*pairEnergy)[i][j][k].size();l++){	
					fprintf(stdout, "    %4d %4d %4d %4d %8.3f", i, j, k, l, (*pairEnergy)[i][j][k][l]);

					if (masks[i][j] && masks[k][l]) {
						fprintf(stdout, " **** ");
					}
					fprintf(stdout,"\n");
				}
			}
		}

	}
}

void LinearProgrammingOptimization::addEnergyTable(vector<vector<double> > &_selfEnergy, vector<vector<vector<vector<double> > > > &_pairEnergy){

	/*
	  Pair Energy Table Layout:
	  
	  Example 3 positions [ 2 1 2 ] = # of rotamers at each position


	  ====== LowerTriangular =====
	  Indices:

	  0 

	  1 -> 0

	       |

	       0 -> 0   1


	 2 -> 0             1
	    
	      |             |
	      
	      0 -> 0 1      0 -> 0 1
	       
	      1 -> 0        1 -> 0
	  

	 so:
	  pairEnergy[0].size() = 0;
	  pairEnergy[1].size() = 1;
	  pairEnergy[1][0].size() = 1;
	  pairEnergy[1][0][0].size() = selfEnergy[0].size();
	  
	  pairEnergy[2].size() = 2;
	  pairEnergy[2][0].size() = 2;
	  pairEnergy[2][0][0].size() = selfEnergy[0].size();
	  pairEnergy[2][0][1].size() = selfEnergy[1].size();
	  pairEnergy[2][1][0].size() = selfEnergy[0].size();
	  pairEnergy[2][1][1].size() = selfEnergy[1].size();



	  ====== UpperTriangular =====
	  Indices:

	  0 -> 0          1            
	            
               |          |            
		    	  	     	       
	       0  	  0            
	       		  	     	       
               1 -> 0	  1 -> 0       
			                       
               2 -> 0 1	  2 -> 0 1     


	  1 -> 0

	       |

	       0
	       
	       1

	       2 -> 0 1


	 2 -> 0             1
	    


	  

	 so:
	  pairEnergy[0].size() = selfEnergy[0].size();
	  pairEnergy[0][0].size() = selfEnergy.size();
	  pairEnergy[0][0][0].size() = 0;
	  pairEnergy[0][0][1].size() = selfEnergy[1].size();
	  pairEnergy[0][0][2].size() = selfEnergy[2].size();
	  pairEnergy[0][1][0].size() = 0;
	  pairEnergy[0][1][1].size() = selfEnergy[1].size();
	  pairEnergy[0][1][2].size() = selfEnergy[2].size();
	  

	  pairEnergy[1].size() = selfEnergy[1].size();
	  pairEnergy[1][0].size() = selfEnergy.size();
	  pairEnergy[1][0][0].size() = 0;
	  pairEnergy[1][0][1].size() = 0;
	  pairEnergy[1][0][2].size() = selfEnergy[2].size();
	  
	  
	  pairEnergy[2].size() = 2;
	  pairEnergy[2][0].size() = 0;
	  pairEnergy[2][1].size() = 0;


	  
	 */

	selfEnergy = &_selfEnergy;


	pairEnergy = new vector<vector<vector<vector<double> > > >();
	pairEnergy->resize(selfEnergy->size());

	// For each position resize to
	// num rotamers
	for (uint i = 0; i < selfEnergy->size();i++){
		(*pairEnergy)[i].resize((*selfEnergy)[i].size());

		// For each rotamer,
		// resize to num positions
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){
			(*pairEnergy)[i][j].resize((*selfEnergy).size());

			// For each position resize to num rotamers
			for (uint k = 0; k < (*selfEnergy).size();k++){
				(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

			}
		}

	}


	// Convert to upper triangular
	for (uint pos1 = 0 ; pos1 < _pairEnergy.size();pos1++){
		for (uint pos2 = pos1+1 ; pos2 < _pairEnergy.size();pos2++){
			for (uint rot2 = 0 ; rot2 < _pairEnergy[pos2].size();rot2++){
				for (uint rot1 = 0 ; rot1 < _pairEnergy[pos2][rot2][pos1].size();rot1++){
					(*pairEnergy)[pos1][rot1][pos2][rot2] =  _pairEnergy[pos2][rot2][pos1][rot1];
				}
			}
		}
	}


}


vector<vector<bool> > LinearProgrammingOptimization::getMask(){
	return masks;
	
}

double LinearProgrammingOptimization::getTotalEnergy(){
	
	double energy = 0.0;
	for (uint i = 0 ;i < totalNumPositions;i++){
		
		energy += (*selfEnergy)[i][rotamerSelection[i]];
		for (uint j = i+1; j < totalNumPositions;j++){

			energy += (*pairEnergy)[j][rotamerSelection[j]][i][rotamerSelection[i]];

		}
	}

	return energy;
}


string LinearProgrammingOptimization::getRotString(){
	string result = "";

	for (uint i = 0; i < totalNumPositions;i++){
			result += MslTools::intToString(rotamerSelection[i])+":";
	}

	return result;
}
