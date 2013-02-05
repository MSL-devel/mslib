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

#include "LinearProgrammingOptimization.h"

#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "MslTools.h"

using namespace MSL;
using namespace std;



LinearProgrammingOptimization::LinearProgrammingOptimization(){

	// Set defaults
	selfEnergy              = NULL;
	pairEnergy              = NULL;
	lp = NULL;
	totalNumRotamers        = 0;
	totalNumPositions       = 0;
	numConstraints          = 0;
	numDecisionVariables    = 0;
	responsibleForEnergyTableMemory = false;
	inputMasks.clear();
	verbose = false;
	
}
LinearProgrammingOptimization::~LinearProgrammingOptimization(){
	deletePointers();
}


void LinearProgrammingOptimization::deletePointers(){
	if(lp) {
		glp_delete_prob(lp);
		lp = NULL;
	}
	deleteEnergyTables();
}

void LinearProgrammingOptimization::deleteEnergyTables(){
	if(responsibleForEnergyTableMemory) {
		if(selfEnergy) {
			delete selfEnergy;
			selfEnergy = NULL;
		}
		if(pairEnergy) {
			delete pairEnergy;
			pairEnergy = NULL;
		}
	}
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
	deleteEnergyTables();
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



			vector<string> toks = MslTools::tokenize(line);

			
			// self Energy Line has 3 numbers on it
			if (toks.size() == 3){
				int posIndex = MslTools::toInt(toks[0]);
				if (selfEnergy->size() < posIndex+1){
					vector<double> tmp;
					selfEnergy->push_back(tmp);
				}

				(*selfEnergy)[posIndex].push_back(MslTools::toDouble(toks[2]));
				
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
					for (uint i = 0; i < pairEnergy->size();i++){
						(*pairEnergy)[i].resize((*selfEnergy)[i].size());

						// For each rotamer,
						// resize to num positions
						for (uint j = 0; j < (*pairEnergy)[i].size();j++){
							(*pairEnergy)[i][j].resize(i);

							// For each position resize to num rotamers
							for (uint k = 0; k < (*pairEnergy)[i][j].size();k++){
								(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

							}
						}

					}
				}
				firstPair = false;
				// Add to pairEnergy object. Indices pairEnergy[POS1][ROT1][POS2][ROT2] = ENERGY;
				int p1 = MslTools::toInt(toks[0]);
				int r1 = MslTools::toInt(toks[1]);
				int p2 = MslTools::toInt(toks[2]);
				int r2 = MslTools::toInt(toks[3]);
				if(p1 > p2) {
					(*pairEnergy)[p1][r1][p2][r2] = MslTools::toDouble(toks[4]);
				} else if(p2 > p1) {
					(*pairEnergy)[p2][r2][p1][r1] = MslTools::toDouble(toks[4]);
				} else {
					// p1 == p2
					cerr << "WARNING 12343: IGNORING LINE: " << line << endl;
				}	
				
			}
		}
	} catch (...){
		cerr << "ERROR 7809 in LinearProgrammingOptimization::readEnergyTable in file: "<<_filename<<endl;
		exit(7809);
	}

}


void LinearProgrammingOptimization::analyzeEnergyTable(){
	if(verbose) {
		cout << "Done reading energy file, now anaylize properites... "<<endl;
	}


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


	// PairType 
	pairType.resize(pairEnergy->size());

	if(verbose) {
		fprintf(stdout, "\tSelfEnergy   : %10d\n",(int)selfEnergy->size());
		fprintf(stdout, "\tPairEnergy   : %10d\n",(int)pairEnergy->size());
	}


	// Load the rotamers of pos0 seperately here
	for (uint rot1 = 0; rot1 < (*selfEnergy)[0].size();rot1++){
		// if the rotamer is not masked out
		if(inputMasks.size() == 0 || inputMasks[0][rot1]) {
			vector<int> tmp;	
			tmp.push_back(0);
			tmp.push_back(rot1);
			rotamerState.push_back(tmp);
		}
	}

	// For each pair of rotamers..
	for (uint pos1 = 0; pos1 < pairEnergy->size();pos1++){

		// Resize pairType[pos1]
		pairType[pos1].resize(pairEnergy->size());


		if(verbose) {
			cout << "PairE_I:"<<pos1<<" "<<(*selfEnergy)[pos1].size()<<endl;
		}
		for (uint pos2 = 0; pos2 < pos1;pos2++){

			if(verbose) {
				cout << "PairE_J:"<<pos2<<" "<<(*selfEnergy)[pos2].size()<<endl;
			}


			// Initialize positive to false, meaning we 
			bool positiveFlag = false;
			for (uint rot1 = 0; rot1 < (*pairEnergy)[pos1].size();rot1++){			

				// addToMask, skip here? numConstraints, rotamerState, rotamerPair
				if (inputMasks.size() > pos1 && inputMasks[pos1].size() > rot1 && !inputMasks[pos1][rot1]){
					continue;
				}
				
				// Keep track of each self-rotamer only the first time over pos2.
				if (pos2 == 0){
					vector<int> tmp;	
					tmp.push_back(pos1);
					tmp.push_back(rot1);
					rotamerState.push_back(tmp);
				}

				for (uint rot2 = 0; rot2 < (*pairEnergy)[pos1][rot1][pos2].size();rot2++){			
					//cout << "PairENERGY: "<<pos1<<","<<rot1<<","<<pos2<<","<<rot2<<endl;
					if (inputMasks.size() > pos2 && inputMasks[pos2].size() > rot2 && !inputMasks[pos2][rot2]){
						continue;
					}
					if ((*pairEnergy)[pos1][rot1][pos2][rot2] > 0){
						positiveFlag = true;
					}

					vector<int> tmp;
					tmp.push_back(pos1);
					tmp.push_back(rot1);
					tmp.push_back(pos2);
					tmp.push_back(rot2);
					rotamerPair.push_back(tmp);
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
		for(uint j = 0; j < (*selfEnergy)[i].size(); j++ ) {
			if(inputMasks.size() == 0 || inputMasks[i][j]) {
				totalNumRotamers++;
			}
		}
	}

	// Keep track of number of constraints
	numConstraints = (totalNumPositions - 1) * totalNumRotamers + totalNumPositions;

	numDecisionVariables = rotamerState.size() + rotamerPair.size(); 
	if(verbose) {
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

}

int LinearProgrammingOptimization::writeCPLEXFile(string _filename) {
	return glp_write_lp(lp,NULL,_filename.c_str());
}

void LinearProgrammingOptimization::solveLP(){

	time_t start,end;
	time(&start);
	int simplexResult = glp_simplex(lp,NULL);
	time(&end);
	double Z = glp_get_obj_val(lp);
	if(verbose) {
		cout << "Simplex Result " << simplexResult << endl;
		cout << "Simplex Time " << difftime(end,start) << endl;
		cout << "Simplex Objective Value " << Z << endl;
	}

	rotamerSelection.clear();
	rotamerSelection.resize(totalNumPositions);
	map<int,double> bestRotWeight; // keep track of the weight of the best rotamer for each position and choose the one with the maximum weight
	for (uint i = 1; i <= rotamerState.size()+rotamerPair.size();i++){
		double x = glp_get_col_prim(lp, i);
		string s = glp_get_col_name(lp, i);
		int kind = glp_get_col_kind(lp, i);

		if(verbose) {
			printf("%10s,%4d  x[%4d] = %8.3f; ", s.c_str(), kind, i,x);
		}

		// 
		if (i <= rotamerState.size()){
			// The self Energy variables
			if(bestRotWeight.find(rotamerState[i-1][0]) == bestRotWeight.end() || bestRotWeight[rotamerState[i-1][0]] < x ) {
				rotamerSelection[rotamerState[i-1][0]] = rotamerState[i-1][1];
				bestRotWeight[rotamerState[i-1][0]] = x;
			}

			if(verbose) {
				printf(" RS");
			}
		}
		if(verbose) {
			printf("\n");
		}

		// Break out of loop if verbose flag is not set
		if (!verbose && i >= rotamerState.size()) break;
	}
}
void LinearProgrammingOptimization::solveMIP(){

	time_t start,end;
	time(&start);
	int intoptResult = glp_intopt(lp,NULL);
	time(&end);
	double Zint = glp_mip_obj_val(lp);
	if(verbose) {
		cout << "Intopt Result " << intoptResult << endl;
		cout << "Intopt Time " << difftime(end,start) << endl;
		cout << "Intopt Objective Value " << Zint << endl;
	}



	rotamerSelection.clear();
	rotamerSelection.resize(totalNumPositions);
	for (uint i = 1; i <= rotamerState.size()+rotamerPair.size();i++){
		double xint = glp_mip_col_val(lp, i);
		string s = glp_get_col_name(lp, i);
		int kind = glp_get_col_kind(lp, i);

		if(verbose) {
			printf("%10s,%4d  x[%4d] = %8.3f; ", s.c_str(), kind, i,xint);
		}

		// 
		if (i <= rotamerState.size()){
			// The self Energy variables
			if (xint){
				rotamerSelection[rotamerState[i-1][0]] = rotamerState[i-1][1];
			}
			if(verbose) {
				printf(" RS");
			}
		}
		if(verbose) {
			printf("\n");
		}

		// Break out of loop if verbose flag is not set
		if (!verbose && i >= rotamerState.size()) break;
	}
}



void LinearProgrammingOptimization::createLP(){

	if (selfEnergy->size() == 0 || pairEnergy->size() == 0){
		cerr << "ERROR 2436 LinearProgrammingOptimization::createLPproblem selfEnergy or pairEnergy lists are NULL or 0 in size: "<<selfEnergy->size()<<" "<<pairEnergy->size()<<endl;
		exit(2436);
	}
	


	// Create an LP problem object, dealloc any memory already allocated
	if(lp) {
		glp_delete_prob(lp);
	}
	lp = glp_create_prob();

	glp_set_prob_name(lp, "rotamerOptimization");

	// Set direction to MINIMIZE (not MAXIMIZE)
	glp_set_obj_dir(lp,GLP_MIN);



	/*
	      Add Constraints/ROWS
	*/

	if(verbose) {
		cout << "Add constraints/rows\n";
	}

	glp_add_rows(lp, numConstraints);
	map<int,vector<int> > indexConstraints;

	// Set up positional constraints
	if(verbose) {
		cout << "\tPositional Constriants\n";
	}
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
	if(verbose) {
		cout << "\tPos-PosRot Constriants\n";
	}
	for (uint pos1 = 0; pos1 < selfEnergy->size();pos1++){
		for (uint pos2 = 0; pos2 < selfEnergy->size();pos2++){
			if (pos1 == pos2) continue;

			for (uint rot2 = 0; rot2 < (*selfEnergy)[pos2].size();rot2++){
				if(inputMasks.size() && !inputMasks[pos2][rot2]) {
					continue;
				}

				vector<int> tmp;
				tmp.push_back(pos1);
				tmp.push_back(pos2);
				tmp.push_back(rot2);
				indexConstraints[index-1] = tmp;

				consName.str("");
				consName << "Pair-"<<pos1<<"_"<<pos2<<"-"<<rot2;
				glp_set_row_name(lp,index,consName.str().c_str());
				glp_set_row_bnds(lp,index,GLP_FX,0,0);
				if (verbose){
					cout <<"\t\t"<<consName.str()<<" ("<<index<<") FIXED TO 0."<<endl;
				}
				
				/*
				// If we get rid of the E(u,v) = 0 variables then some constraints become upper bounded
				bool result = false;
				if (pos2 > pos1){
					result = pairType[pos2][pos1];
				} else {
					result = pairType[pos1][pos2];
				}
				if (result){
					glp_set_row_bnds(lp,index,GLP_FX,0,0);
					if (verbose){
						cout <<"\t\t"<<consName.str()<<" ("<<index<<") FIXED TO 0."<<endl;
					}
				} else {
					glp_set_row_bnds(lp,index,GLP_UP,0,0);

					if (verbose){
						cout <<"\t\t"<<consName.str()<<" ("<<index<<") BOUNDED UPTO 0."<<endl;
					}
				}
				*/
				
				index++;


			}
		}
	}

	
	/*
	     Add decision variables (one per rotamer pair).
	          rotamer state 
		  rotamer pair
	*/

	if(verbose) {
		cout << "Add decision variables\n";
	}
	stringstream colName;
	glp_add_cols(lp, rotamerState.size()+rotamerPair.size());


	// Add rotamer selected
	index = 1;
	if(verbose) {
		cout << "\tRotamer States\n";
	}
	for (uint i = 0 ;i < rotamerState.size();i++){

		colName.str("");
		colName << "RotState-"<<rotamerState[i][0]<<"-"<<rotamerState[i][1];
		glp_set_col_name(lp,index,colName.str().c_str());
		glp_set_col_bnds(lp,index,GLP_DB,0.0,1.0); 
		glp_set_col_kind(lp,index,GLP_IV); // The simplex and interior point solvers ignore it anyway
		glp_set_obj_coef(lp,index++,(*selfEnergy)[rotamerState[i][0]][rotamerState[i][1]]);

		if (verbose){
			cout << "\t\tState "<<i<<"("<<index-1<<") : "<<(*selfEnergy)[rotamerState[i][0]][rotamerState[i][1]]<<endl;
		}
	}


	// Add pair/edge selected
	if(verbose) {
		cout << "\tRotamer Pairs\n";
	}
	for (uint i = 0; i < rotamerPair.size();i++){

		colName.str("");
		colName << "RotPair-"<<rotamerPair[i][0]<<"-"<<rotamerPair[i][1]<<"_"<<rotamerPair[i][2]<<"-"<<rotamerPair[i][3];
		glp_set_col_name(lp,index,colName.str().c_str());
		glp_set_col_bnds(lp,index,GLP_DB,0.0,1.0); 
		glp_set_col_kind(lp,index,GLP_IV); 

		if (rotamerPair[i][0] > rotamerPair[i][2]) {
			glp_set_obj_coef(lp,index,(*pairEnergy)[rotamerPair[i][0]][rotamerPair[i][1]][rotamerPair[i][2]][rotamerPair[i][3]]);

			if (verbose){
				cout << "\t\t"<<colName.str()<<" "<<i<<"("<<index<<") : "<<(*pairEnergy)[rotamerPair[i][0]][rotamerPair[i][1]][rotamerPair[i][2]][rotamerPair[i][3]]<<endl;
			}
		} else {
			// This is the greater than case. We should not have an equal case
			glp_set_obj_coef(lp,index,(*pairEnergy)[rotamerPair[i][2]][rotamerPair[i][3]][rotamerPair[i][0]][rotamerPair[i][1]]);

			if (verbose){
				cout << "\t\t"<<colName.str()<<" "<<i<<"("<<index<<") : "<<(*pairEnergy)[rotamerPair[i][2]][rotamerPair[i][3]][rotamerPair[i][0]][rotamerPair[i][1]]<<endl;
			}
		}

		index++;

	}

	//cout << "Num Integer Columns: "<<glp_get_num_int(lp)<<endl;


	// For each constraint/row
	if(verbose) {
		cout << "Compute non-zero coefficent matrix \n";
	}
	int nonZeroCount = 0;
	vector<vector<double> > coeff;
	for (uint i = 0; i < numConstraints;i++){
		/*
		if (i < totalNumPositions){
			cout << indexConstraints[i][0] << endl;
		} else {
			cout << indexConstraints[i][0] << " " << indexConstraints[i][1] << " " << indexConstraints[i][2] << endl ;
		}
		*/
		for (uint j = 0; j < numDecisionVariables;j++){
			if (i < totalNumPositions){
				// all the self rotamer constraints
				if (j < rotamerState.size() && rotamerState[j][0] == indexConstraints[i][0]){
					// only for rotamerState.size() number of decision variables
					nonZeroCount++;
					vector<double> tmp;
					tmp.push_back(i+1);
					tmp.push_back(j+1);
					tmp.push_back(1);
					coeff.push_back(tmp);
					//cout << "CONSTRAINT " << i+1 << "," << j+1 << ",1" << endl;
				}
			} else {
				// all the pair constraints
				if (j >= rotamerState.size() && rotamerPair[j-rotamerState.size()][0] == indexConstraints[i][0] &&
				                                rotamerPair[j-rotamerState.size()][2] == indexConstraints[i][1] &&
				                                rotamerPair[j-rotamerState.size()][3] == indexConstraints[i][2]){
					vector<double> tmp;
					tmp.push_back(i+1);
					tmp.push_back(j+1);
					tmp.push_back(1);
					coeff.push_back(tmp);
					//cout << "CONSTRAINT " << i+1 << "," << j+1 << ",1" << endl;
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
						//cout << "CONSTRAINT " << i+1 << "," << j+1 << ",1" << endl;
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
					//cout << "CONSTRAINT " << i+1 << "," << j+1 << ",-1" << endl;
				}
			}
		}
	}

	if(verbose) {
		cout << "Create coefficent matrix ("<<(nonZeroCount+1)<<")\n";
	}
	int row[(nonZeroCount+1)];
	int col[(nonZeroCount+1)];
	double val[(nonZeroCount+1)];

	index = 1;

	for (uint i = 1; i <= nonZeroCount;i++){
		row[i] = (int)coeff[i-1][0];
		col[i] = (int)coeff[i-1][1];
		val[i] = coeff[i-1][2];
	}




	if(verbose) {
		cout << "Load Matrix"<<endl;
	}
	// Load coefficient matrix from constraints
	glp_load_matrix(lp, nonZeroCount, row,col,val);

}


void LinearProgrammingOptimization::printMe(bool _selfOnly){

	fprintf(stdout,"Self terms:\n");
	for (uint i = 0; i < (*selfEnergy).size();i++){
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){

			fprintf(stdout, "    %4d %4d %8.3f", i, j, (*selfEnergy)[i][j]);
			if (rotamerSelection[i] == j) {
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

					if (rotamerSelection[i] == j && rotamerSelection[k] == l) {
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

	deleteEnergyTables();
	selfEnergy = &_selfEnergy;
	pairEnergy = &_pairEnergy;

}

vector<vector<bool> > LinearProgrammingOptimization::getMask() {
	vector<vector<bool> > aliveRotamers;
	aliveRotamers.resize(totalNumPositions);
	for(int i = 0; i < aliveRotamers.size(); i++) {
		aliveRotamers[i].resize((*selfEnergy)[i].size(),false);
	}
	for(int i = 0; i < totalNumPositions; i++) {
		aliveRotamers[i][rotamerSelection[i]] = true;
	} 
	return aliveRotamers;
}



double LinearProgrammingOptimization::getTotalEnergy(){
	
	double energy = 0.0;
	for (uint i = 0 ;i < totalNumPositions;i++){
		
		energy += (*selfEnergy)[i][rotamerSelection[i]];
		for (uint j = 0; j < i;j++){

			energy += (*pairEnergy)[i][rotamerSelection[i]][j][rotamerSelection[j]];

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

vector<unsigned int>& LinearProgrammingOptimization::getSolution(bool _runMIP) {
	analyzeEnergyTable();
	createLP();
	// The LP needs to be run even if we plan to run MIP
	solveLP();
	if( _runMIP) {
		solveMIP();
	}
	return rotamerSelection;
}
