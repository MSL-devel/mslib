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

#include "Predicate.h"
#include "MslTools.h"

using namespace MSL;
using namespace std;


Predicate::Predicate(){
	init();
}

Predicate::~Predicate(){
}

void Predicate::init(){
	op = "";
	operands.clear();
	fullText = "";
	results.clear();
}

void Predicate::parseAsInFixText(vector<string> &_validOperators){

	/*
	  Better way to do it:
	   read each char, add to a token
	   if whitespace, then check if token is a valid operator.
	   if it is a valid operator, then 
	          set the op value
		  set remaining chars(-op chars) from token to operand1
		  set flag operatorFound = true
	 */


	//int operandCount = 0;
	//bool inOperand = false;
	string token = "";
	op = "";

	bool operatorFound = false;
	// Parse each line of the text
	for (uint i = 0 ; i < fullText.size();i++){

		//cout << "Character: "<<fullText.substr(i,1)<<endl;
		if (i > 0 && fullText.substr(i,1) == " "){
			
			operatorFound = false;
			for (uint i = 0; i < _validOperators.size();i++){
				int opLength = _validOperators[i].size();

				if (opLength >= token.size()) continue;
				if (MslTools::toUpper(token.substr(token.size()-opLength,opLength)) == _validOperators[i]){
					op = _validOperators[i];
					if (trim(token.substr(0,token.size()-opLength-1)) != ""){
							operands.push_back(token.substr(0,token.size()-opLength-1)); // operand1
					}
					operatorFound = true;
					token = "";
					break;
				}
			}

			if (operatorFound){
				continue;
			}
		}

		
		token += fullText.substr(i,1);
		
	}

	if (trim(token) != ""){
		operands.push_back(token); // operand2
	}



	// Test for data filled ?
	//  op
	//  operand1
	//  operand2
	//cout << toString();


}
void Predicate::parseAsPostFixText(vector<string> &_validOperators){

	/*
	int operandCount = 0;
	bool inOperand = false;
	string token = "";
	for (uint i = 0 ; i < fullText.size();i++){


		if (fullText[i] == '\''){

			// Switch inOperand flag
			if (inOperand){

				// Add text to approriate operand
				if (operand1 == ""){
					operand1 = token;
				} else {
					if (operand2 == ""){
						operand2 = token;
					}
				}


				operandCount++;
				inOperand = false;
				token = "";
			} else {
				inOperand = true;

			}	

			continue;
		}


		token += fullText[i];
	}

	token = trim(token);
	//token = transform(token.begin(),token.end(),token.begin(),toupper);

	// Check which operator
	op = "";
	for (uint i = 0; i < _validOperators.size();i++){
		if (token == _validOperators[i]){
			op = _validOperators[i];
			break;
		}
	}

	if (op == ""){
		cerr << "ERROR 2341 Predicate::parsePostFixText invalid operator: "<<token<<endl;
		exit(2341);
	}
	*/

}


string Predicate::toString(){
	stringstream sout;
	if (operands.size() == 2){
		//sout << "(Operand,Operator, Operand) = ("<<operands[0]<<","<<op<<","<<operands[1]<<")"<<endl;
		sout << "("<<operands[0]<<","<<op<<","<<operands[1]<<")";
	} else if (operands.size() == 1){
		//sout << "(Operand,Operator, TREE) = ("<<operands[0]<<","<<op<<",TREE)"<<endl;
		sout << "("<<operands[0]<<","<<op<<",TREE)";
	} else {
		//sout << "(TREE, Operator, TREE) = (TREE,"<<op<<",TREE)"<<endl;
		sout << "(TREE,"<<op<<",TREE)";
	}
	return sout.str();
}
string Predicate::trim(string str){
  // trim leading whitespace
  string::size_type  notwhite = str.find_first_not_of(" \t\n\r");
  str.erase(0,notwhite);

  // trim trailing whitespace
  notwhite = str.find_last_not_of(" \t\n\r"); 
  str.erase(notwhite+1); 

  return str;
}

vector<string> Predicate::tokenizeWord(string _input, string _delimiter){

  vector<string> results;

  int start  = _input.find_first_not_of(_delimiter);
  int end    = 0;
  int loop   = 0;
  while ( start != std::string::npos){
    if (loop == 0) start = 0;
    end    = _input.find(_delimiter, start);
    results.push_back(_input.substr(start, end-start));
    start  = _input.find_first_not_of(_delimiter,end);
    loop++;
  }

  return results;
}
