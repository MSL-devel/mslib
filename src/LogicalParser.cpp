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

#include "LogicalParser.h"
#include "MslTools.h"



#include <sstream>
#include <queue>

using namespace MSL;
using namespace std;


LogicalParser::LogicalParser(){
	logicStatementInFix = "";
	logicStatementPostFix = "";
	treeRoot          = NULL;
	validOperators.push_back("AND");
	validOperators.push_back("OR");
	validOperators.push_back("XOR");
	debug = false;

}

LogicalParser::~LogicalParser(){

	// DONE    TODO : Need to delete Tree from treeRoot... hmmm... 
	delete treeRoot;
}

void LogicalParser::setLogicStatementInFix(string _infix){
	logicStatementInFix = _infix;

	logicStatementInFix  = MslTools::trim(logicStatementInFix);
	if (logicStatementInFix[0] != '('){
		logicStatementInFix = "(" + logicStatementInFix + ")";
	}
 
	
}

void LogicalParser::parse(){

	// Current Tree pointer
	Tree<Predicate> *current;

	// Remove whitespace from back and front of statement
	logicStatementInFix = MslTools::toUpper(MslTools::trim(logicStatementInFix));
	//cout << "Input Statement: "<<logicStatementInFix;

	// TODO: A memory leak?  .. This means you can't use a LogicalParser object 2 times. NOT BUG, BUT POOR BEHAVIOR.
	if (treeRoot != NULL){
		cerr << "ERROR 1234 LogicalParser::parse() treeRoot is not NULL.  Setting to NULL will generate a memory leak.\n";
		exit(1234);
	}

	// Initialize tree
	treeRoot = NULL;  
	current  = NULL;


	// Search string for indices of operators
	map<int,bool> validOpIndices;
	for (uint j = 0; j < validOperators.size();j++){
		stringstream ss;
		ss << " "<<validOperators[j]<<" ";
		int start  = logicStatementInFix.find_first_not_of(ss.str());
		int loop = 0;
		while (start != std::string::npos){
			if (loop == 0) start = 0;
			int pos    = logicStatementInFix.find(ss.str(), start);
			if (pos >= 0){
				validOpIndices[pos] = true;
				if (debug){
					cout << "Operator("<<ss.str()<<") index is: "<<pos<<" start: "<<start<<endl;
				}
			}
			start      = logicStatementInFix.find_first_not_of(ss.str(), pos);
			loop++;
		}
	}

//	if (validOpIndices.size() == 0){
//		cerr << "ERROR 1234 LogicalParser::parse no valid OPERATORS in statement: "<<logicStatementInFix<<endl;
//		exit(1234);
//	}

	// Counters for our next loop
	int operatorCount         = 0;
	int openParenthesisCount  = 0;
	int closeParenthesisCount = 0;
	string correctFormat =  "";
	openParenthesisCount = closeParenthesisCount =0;

	if (debug){
		cout << "Fix this please: "<<logicStatementInFix<<endl;
	}


	// Check for problems with string before we do anything
	//  TO MANY BUGS... REQUIRE GOOD LOGICAL STATEMENT!


// 	for (uint i = 0 ; i < logicStatementInFix.size();i++){

// 		// Keep track of # open,close parens. At end of loop, we can die if counts are not equal.
// 		cout << "I["<<i<<"]: "<<logicStatementInFix[i]<<endl;
// 		if (logicStatementInFix[i] == '(') openParenthesisCount++;
// 		if (logicStatementInFix[i] == ')') closeParenthesisCount++;

// 		// Increment operator count if we have an operator index  (i).
// 		//    Also check for being in between ( ) . In that case don't add to count.
// 		if (validOpIndices[i] && openParenthesisCount != closeParenthesisCount+1) {

// 			operatorCount++; 

// 			if (debug){
// 				cout << "Index: "<<i<<" -> Increment OpCount("<<operatorCount<<"), have validOP and there are "<<openParenthesisCount<<" '(' and "<<closeParenthesisCount<<" ')'"<<endl;
// 			}
// 		}

// 		// PROBLEM: found 2 operators within '( )' term.  Add in '(' and ')'.
// 		if (operatorCount > 1){
// 			correctFormat = "(" + correctFormat + ") ";
// 			operatorCount = 1;
// 			continue;
// 		}


// 		correctFormat += logicStatementInFix[i];

		
// 	}

// 	// Exit if parenthesis count does not match (i.e. poorly formed selection string!)
// 	if (openParenthesisCount != closeParenthesisCount){
// 		cerr << "ERROR 3434 Parenthesis count are not equal : (open,close) ("<<openParenthesisCount<<","<<closeParenthesisCount<<") for statment: "<<logicStatementInFix<<endl;
// 		exit(3434);
// 	}


// 	//cout << "\tOpCount: "<<operatorCount<<" "<<openParenthesisCount<<" Corrected Statement: '"<<correctFormat<<"'"<<endl;

// 	// Reset logic statement
// 	logicStatementInFix = "(" + correctFormat + ")";



	// Special Case, no operators only have single '(' and ')'
// 	if (operatorCount == 0){
// 		std::replace(logicStatementInFix.begin(),logicStatementInFix.end(),'(',' ');
// 		std::replace(logicStatementInFix.begin(),logicStatementInFix.end(),')',' ');
// 		logicStatementInFix = "("+logicStatementInFix+")";
		
// 	}

	if (debug){
		cout << "LOGIC STATEMENT: "<<logicStatementInFix<<" "<<operatorCount<<endl;
	}
	// For each character
	for (uint i = 0 ; i < logicStatementInFix.size();i++){

		// Break loop if current pointer is NULL
		if (i > 0 && current == NULL) break;

		//cout << "Character: "<<logicStatementInFix[i]<<endl;

		// For each '(' create a new branch of the tree
		if (logicStatementInFix[i] == '('){

			// Create the tree root if tree root doesn't exist
			if (treeRoot == NULL) {
				treeRoot = new Tree<Predicate>();
				current = treeRoot;
				predicateList.push_back(current->getData());

				if (debug){
					cout << "tree root created: ->"<<current->getData()->getText()<<"<-"<<endl;
				}
			} else {

				// Add a subtree to the current tree
				current->addSubTree(new Tree<Predicate>(current));


				if (debug){
					cout << "Subtree added: ->"<<current->getData()->getText()<<"<-"<<endl;
				}
				// If no text yet, SubTree is 0.  Else Subtree is 1.
				if (current->getData()->getText() == ""){
					//cout << "\tBranch, SubTree1 (or Predicate1)"<<endl;
					current = (*current)[0]; // subTree 0 
				} else {
					//cout << "\tBranch, SubTree2 (or Predicate2) "<<endl;
					current = (*current)[1]; // subTree 1 
				}
				predicateList.push_back(current->getData());

			}
			continue;
		}
			

		// Ending logic term , parse the text into operand, operand, operator and move back up the tree
		if (logicStatementInFix[i] == ')'){
			if (debug){
				cout << "Moving back up the tree"<<endl;
			}
			(current->getData())->parseAsInFixText(validOperators);
			current = current->getParent();
			continue;
		}		


		// Add characters to growing text for current logic term
		(current->getData())->addToText(logicStatementInFix.substr(i,1));
		

	}


}


/*
  Logic examples

 (A and B or (C and D) or E)

 ( ((A and B) or (C and D)) or E )


 ((('A' 'B' and) ('C' 'D' and) or) 'E' or)

       Will be:
                              OR
			     /  \
			    OR  E
			   /  \
		       AND    AND
 		       / \    / \
		      A   B  C   D


 ( 'E' (('A' 'B' and) ('C' 'D' and) or) or)
 
 */


/*
  Print out by traversing breadth-first..
 */
void LogicalParser::printLogicTree(){

	treeRoot->printTree(100);
}

bool LogicalParser::eval(KeyLookup &_aLookupObject){

	for (uint i = 0; i < predicateList.size();i++){
		predicateList[i]->clearResults();
	}
	
	// No recusive algorithm for a single node...(treeRoot only)
	if (predicateList.size() == 1){


		evalOperand(_aLookupObject, *treeRoot->getData(),0);

		if (treeRoot->getData()->getOperator() != ""){
			evalOperand(_aLookupObject, *treeRoot->getData(),1);
		}


		if (treeRoot->getData()->getOperator() == "AND"){
			return treeRoot->getData()->getResult(0) && treeRoot->getData()->getResult(1);
		} 

		if (treeRoot->getData()->getOperator() == "OR"){
			return treeRoot->getData()->getResult(0) || treeRoot->getData()->getResult(1);
		}


		return treeRoot->getData()->getResult(0);
	}



	// A recursive algorithm...
	return recursiveEval(_aLookupObject, treeRoot);
}

bool LogicalParser::recursiveEval(KeyLookup &_aLookupObject, Tree<Predicate> *_node){
	
	Predicate *cur  = _node->getData();
	int numSubTrees = _node->getNumSubtrees();
	int numResults;
	numResults = cur->getNumResults();

	bool recursiveResult = false;
	//cout << "Current: "<<_node->getData()->toString()<<" subTrees: "<<numSubTrees<<" results: "<<cur->getNumResults()<<endl;
	if (numSubTrees == 2 && cur->getNumResults() == 0){
		recursiveResult = recursiveEval(_aLookupObject, (*_node)[0]);
		cur->addResult(recursiveResult);
	}

	if (numSubTrees == 2 && cur->getNumResults() == 1){
		int result = bailOutEarly(*cur);

		if (result == 0 || result == 1){
			return result;
		}

		recursiveResult = recursiveEval(_aLookupObject, (*_node)[1]);
		cur->addResult(recursiveResult);
	}


	if (numSubTrees == 1 && cur->getNumResults() == 0){
		evalOperand(_aLookupObject, *cur,0);
	}

	if (numSubTrees == 1 && cur->getNumResults() == 1){

		int result = bailOutEarly(*cur);
		

		if (result == 0 || result == 1){
			return result;
		}

		recursiveResult = recursiveEval(_aLookupObject, (*_node)[0]);
		cur->addResult(recursiveResult);
	}

	
	if (numSubTrees == 0 && cur->getNumResults() == 0){
		evalOperand(_aLookupObject, *cur,0);
	}

	if (numSubTrees == 0 && cur->getNumResults() == 1){
		int result = bailOutEarly(*cur);
		

		if (result == 0 || result == 1){
			return result;
		}


		
		evalOperand(_aLookupObject, *cur,1);
	}

	// Final case
	if (cur->getNumResults() == 2){
		
		// This is no good. there should be no hard coding of
		// operators.
		if (cur->getOperator() == "AND") {
			return cur->getResult(0) && cur->getResult(1);
		} else if (cur->getOperator() == "OR") {
			return cur->getResult(0) || cur->getResult(1);
		} else {
			cerr << "OPERATOR: "<<cur->getOperator()<<" was not found."<<endl;
		}
	}
	
	// IS THIS RIGHT? OR SHOULD IT BE FALSE?
	return true;
}


int LogicalParser::bailOutEarly(Predicate &_predObj){

	if (_predObj.getOperator() == "OR" && _predObj.getResult(0)) return 1;

	if (_predObj.getOperator() == "AND" && !_predObj.getResult(0)) return 0;

	return -2;
}
void LogicalParser::evalOperand(KeyLookup &_aLookupObject,Predicate &_predObj, int _operand){


	vector<string> toks;
	toks = MslTools::tokenize(MslTools::trim(_predObj.getOperand(_operand))," ");

	string keyword = toks[0];

	bool flipResult = false;
	if (toks[0] == "NOT"){
		flipResult = true;
		keyword = toks[1];
	}

	string keywordType = _aLookupObject.isValidKeyword(keyword);
	bool result = false;
	if ( keywordType != ""){


		string value   = toks[1];
		if (flipResult) {
			value   = toks[2];
		}


		vector<string> multipleValues = MslTools::tokenize(value, "+");
		for (uint i = 0; i < multipleValues.size();i++){

			
			if (keywordType == "string"){
				if (_aLookupObject.getString(keyword) == multipleValues[i]){
					result  = true;
				}
			}

			if (keywordType == "real"){
				double diff = _aLookupObject.getReal(keyword)- MslTools::toDouble(multipleValues[i]);
				if ( diff < 0.00001 && diff > -0.00001){
					result  = true;
				}
			}

			if (keywordType == "int"){

				vector<string> range = MslTools::tokenize(multipleValues[i],"-");
				if (range.size() == 1 && _aLookupObject.getInt(keyword) == MslTools::toInt(multipleValues[i])){
					result  = true;
				}

				if (range.size() == 2 && _aLookupObject.getInt(keyword) >= MslTools::toInt(range[0]) && _aLookupObject.getInt(keyword) <= MslTools::toInt(range[1])){
					result = true;
				}
				
			}

			if (keywordType == "bool"){

				if ( _aLookupObject.getBool(keyword) == MslTools::toInt(multipleValues[i])){
					result = true;
				}
			}

			if (keywordType == "queryBool"){
			  if ( _aLookupObject.getQueryBool(keyword,multipleValues[i])){
			       result = true;
			     }
			  
			}

			if (result) break;
		}


	} else {

		if (toks.size() == 1) {

			// Try using this keyword as a name in selectionFlags hash.
			if (_aLookupObject.getSelectionFlag(keyword)){
				result = true;
			} else {
				result = false;

			}
		} else {

			cerr << "Keyword "<<keyword<<" is invalid."<<endl;
		}


	}


	if (flipResult){
		result = !result;
	}


	_predObj.addResult(result); 

	/*
	if (result){
		cout << "Operand: "<<_predObj.getOperand(_operand)<< " was "<<result<<endl;
	}
	*/


}




vector<string> LogicalParser::tokenize(string input) {
	input = MslTools::trim(input);
	vector<string> matches;
	vector<string>::iterator k;

	bool open = false;
	for (int i=0; i< input.size(); i++) {
		if (input.substr(i,1) == " " || input.substr(i,1) == "\t" || input.substr(i,1) == "\n") {
			open = false;
			continue;
		} else {
			if (open) {
				*k += input[i];
			} else {
				open = true;
				matches.push_back(input.substr(i,1));
				k = matches.end() - 1;
			}
		}
	}
	return matches;
}

