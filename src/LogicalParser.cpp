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

#include "LogicalParser.h"
#include "MslTools.h"



#include <sstream>
#include <queue>
#include <stack>

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("LogicalParser");

LogicalParser::LogicalParser(){
	logicStatementInFix = "";
	logicStatementPostFix = "";
	treeRoot          = NULL;
	validOperators.push_back("AND");
	validOperators.push_back("XOR");
	validOperators.push_back("OR");

	
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
/*
  Parse:
  1 condition:
  (condition)

  2 conditions:
  (condition or condition)

  3 conditions:
  (condition or condition and conidition)
  ((condition and conidition) or condition)
  (condition or (condition and conidition))

  4 conditions:
  (condition or condition and condition or condition)
  ((condition and condition) or condition or condition)
  (condition or (condition and condition) or condition)
  (condition or condition or (condition and condition)) 

 */
void LogicalParser::parse(){

	// Current Tree pointer
	Tree<Predicate> *current;

	// Remove whitespace from back and front of statement
	logicStatementInFix = MslTools::toUpper(MslTools::trim(logicStatementInFix));
	MSLOUT.stream() << "Input Statement: "<<logicStatementInFix<<endl;

	// TODO: A memory leak?  .. This means you can't use a LogicalParser object 2 times. NOT BUG, BUT POOR BEHAVIOR.
	if (treeRoot != NULL){
		cerr << "ERROR 1234 LogicalParser::parse() treeRoot is not NULL.  Setting to NULL will generate a memory leak.\n";
		exit(1234);
	}

	// Initialize tree
	treeRoot = NULL;  
	current  = NULL;

	stack<string> orderTree;
	string currentStr = "";
	for (uint i = 0 ; i < logicStatementInFix.size();i++){
	  
	  // Add to stack if open paren
	  if (logicStatementInFix[i] == '('){

	    if (currentStr != ""){
	      orderTree.push(currentStr);
	      currentStr = "";
	    }

	    orderTree.push(logicStatementInFix.substr(i,1));
	    continue;
	  }

	  // Add to tree if close paren
	  if (logicStatementInFix[i] == ')'){

	    createPredicate(orderTree, currentStr, &treeRoot, &current);

	    // add predicate back to tree
	    
	    continue;
	  }


	  if (logicStatementInFix[i] == ' '){
	    orderTree.push(currentStr);
	    currentStr = "";
	    continue;
	  }

	  currentStr += logicStatementInFix.substr(i,1);


	} // For each charater in logic statement
	
	// Must flush out stack here...
	while (!orderTree.empty()){
	  createPredicate(orderTree, currentStr, &treeRoot, &current);
	}

  
}
void LogicalParser::createPredicate(stack<string> &orderTree, string &currentStr, Tree<Predicate> **treeRoot, Tree<Predicate> **current ){

            if (currentStr != ""){
	      orderTree.push(currentStr);
	    }
	    currentStr = "";
	    
	    // Unstack until '(' or operatorCount == 2
	    string predicateString = "";
	    int operatorCount = 0;
	    while (!orderTree.empty()){
	      string topStr = orderTree.top();
	      orderTree.pop();
	      
	      // Check if this is an operator...
	      bool validOperator = false;
	      for (uint v = 0; v < validOperators.size();v++){
		//int opLength = validOperators[v].size();
		//if (opLength >= topStr.size()) continue;

		if (MslTools::toUpper(topStr) == validOperators[v]){
		  validOperator = true;
		  break;
		}
	      }
	      //cout << "Considering topStr: "<<topStr<<" validOp: "<<validOperator<<" opCount: "<<operatorCount<<endl;
	      if (validOperator){
		operatorCount++;
		/*

		  // Here is where ordering of operators would go...

		localStack_Operands.push(predicateString);
		predicateString = "";

		// If this operator has lower preference than the top of the stack, invert them (plus the operand stack)
		if (localStack_Operators.top()

		localStack_Operators.push(validOperator);
		*/
	      }

	      if (operatorCount == 2){
		orderTree.push(topStr);
		break;
	      }

	      if (topStr == "("){
		break;
	      }


	      if (predicateString == ""){
		predicateString = topStr;
	      } else {
		predicateString = topStr + " " + predicateString;
	      }

	    } // While stack not empty

	    predicateString += " ";
	    MSLOUT.stream() << "Creating predicate from "<<predicateString<<"."<<endl;

	    // We are now on a single level , order the operators
	    stack<string> operators;
	    stack<string> operands;
	    string currentItem = "";
	    for (uint i = 0; i < predicateString.size();i++){
	      
	    }
	    // Add each one to this subtree
	    


	    // Parse lastString  ... assume at most 1 operand, 2 operators
	    Predicate *p = new Predicate();
	    p->setText(predicateString);
	    p->parseAsInFixText(validOperators);
	    *treeRoot = new Tree<Predicate>();
	    (*treeRoot)->setData(*p);
	    delete(p);

	    if ( (*current) != NULL){
	      (*treeRoot)->addSubTree(*current);
	      (*current)->setParent(*treeRoot);
	    }

	    *current = *treeRoot;
	    predicateList.push_back((*current)->getData());
}
void LogicalParser::parse_almost(){

	// Current Tree pointer
	Tree<Predicate> *current;

	// Remove whitespace from back and front of statement
	logicStatementInFix = MslTools::toUpper(MslTools::trim(logicStatementInFix));
	MSLOUT.stream() << "Input Statement: "<<logicStatementInFix;

	// TODO: A memory leak?  .. This means you can't use a LogicalParser object 2 times. NOT BUG, BUT POOR BEHAVIOR.
	if (treeRoot != NULL){
		cerr << "ERROR 1234 LogicalParser::parse() treeRoot is not NULL.  Setting to NULL will generate a memory leak.\n";
		exit(1234);
	}

	// Initialize tree
	treeRoot = NULL;  
	current  = NULL;


	// Search string for indices of operators
	map<int,int> validOpIndices;
	int lastOpIndex =0;

	for (uint j = 0; j < validOperators.size();j++){
		stringstream ss;
		ss << " "<<validOperators[j]<<" ";
		int start  = logicStatementInFix.find_first_not_of(ss.str());
		int loop = 0;
		while (start != std::string::npos){
			if (loop == 0) start = 0;
			int pos    = logicStatementInFix.find(ss.str(), start);
			if (pos >= 0){
			        validOpIndices[pos] = validOperators[j].size();
				if (pos > lastOpIndex){
				  lastOpIndex = pos;
				}
				if (debug){
					MSLOUT.stream() << "Operator("<<ss.str()<<") index is: "<<pos<<" start: "<<start<<endl;
				}
			}
			start      = logicStatementInFix.find_first_not_of(ss.str(), pos);
			loop++;
		}
	}


	MSLOUT.stream() << "LAST OP INDEX: "<<lastOpIndex<<endl;
	int numberOfSeqOperators = 0;
	vector<pair<int,int> > op_positions_in_subtext;
	for (uint i = 0 ; i < logicStatementInFix.size();i++){
	  
	        // Parenthesis-defined End of a logic subtree
		if (logicStatementInFix[i] == ')'){
		  MSLOUT.stream() << "I: "<<i<<" Ending logic statement parse it1: "<<current->getData()->toString()<<endl;
		  (current->getData())->parseAsInFixText(validOperators);
		  MSLOUT.stream() << "I: "<<i<<" Ending logic statement parse it2: "<<current->getData()->toString()<<endl;
		  if (current == treeRoot) break;
		  current = current->getParent();

		  numberOfSeqOperators=0;
		  op_positions_in_subtext.clear();
		  //treeRoot->printTree();

		  // New subtree because of sequential seq operations disrupted by parenthesis (name CA or (name CB and chain B) or resn thr)
		  //if (current->getNumSubtrees() > 0 && i+1 < lastOpIndex){
		  if (i+1 <= lastOpIndex){
		    MSLOUT.stream() << "ADDING ADDITIONAL SUBTREE"<<endl;
		    current->addSubTree(new Tree<Predicate>(current));
		    current = (*current)[current->getNumSubtrees()-1]; // new subtree
		  }

		  continue;
		}
	  
	  	// Parenthesis-defined New logic subtree
		if (logicStatementInFix[i] == '('){

		         if (treeRoot == NULL){
			        MSLOUT.stream() << "Creating root of tree"<<endl;
				treeRoot = new Tree<Predicate>();
				current = treeRoot;
				predicateList.push_back(current->getData());
			 } else {


			   // Add a subtree to the current tree
			   MSLOUT.stream() << "Add subtree"<<endl;
			   current->addSubTree(new Tree<Predicate>(current));

			   if (numberOfSeqOperators == 1){
			     string subStatement = (current->getData())->getText();			     
			     (current->getData())->setText(subStatement.substr(0,op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+2));
			     (current->getData())->parseAsInFixText(validOperators);


			     current = (*current)[0]; // subTree 1 
			     numberOfSeqOperators = 0;
			     op_positions_in_subtext.clear();
			     //treeRoot->printTree();
			   } else {

			     // If no text yet, SubTree is 0.  Else Subtree is 1.
			     if (current->getData()->getText() == ""){
			       MSLOUT.stream() << "\tBranch, SubTree1 (or Predicate1)"<<endl;
			       current = (*current)[0]; // subTree 0 
			     } else {
			       MSLOUT.stream() << "\tBranch, SubTree2 (or Predicate2) "<<endl;
			       current = (*current)[1]; // subTree 1 
			     }
			   }
			   predicateList.push_back(current->getData());
			   numberOfSeqOperators = 0;
			   op_positions_in_subtext.clear();
			   //treeRoot->printTree();		 
			 }

			continue; 

		}

		map<int,int>::iterator it;
		it = validOpIndices.find(i);
		if (it != validOpIndices.end()) { 
		  numberOfSeqOperators++;
		  op_positions_in_subtext.push_back(pair<int,int>((current->getData())->getText().size(),it->second));
		}


		// New subtree because of sequential seq operations
		if (numberOfSeqOperators == 2) {

		    string subStatement = (current->getData())->getText();

		    MSLOUT.stream() << "Substatement: "<<subStatement<<" indices: "<<op_positions_in_subtext[0].first<<" "<<op_positions_in_subtext[0].second<<" "<<subStatement.size()<<" next: "<<subStatement.substr(op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+1)<<"."<<" substr: "<<subStatement.substr(0,op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+2)<<"."<<endl;
		    // Use first operand as text
		    //(current->getData())->setText(subStatement.substr(0,op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+1));
		    (current->getData())->setText(subStatement.substr(0,op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+2));
		    (current->getData())->parseAsInFixText(validOperators);

		    //(current->getData())->addOperand(subStatement.substr(0,op_positions_in_subtext[0].first));
		    //(current->getData())->setOperator(subStatement.substr(op_positions_in_subtext[0].first,op_positions_in_subtext[0].second));

		    // Add a subtree to the current tree
 		    current->addSubTree(new Tree<Predicate>(current));
		    current = (*current)[current->getNumSubtrees()-1]; // new subtree
		    predicateList.push_back(current->getData());
		    current->getData()->setText(subStatement.substr(op_positions_in_subtext[0].first+op_positions_in_subtext[0].second+1));
		    (current->getData())->addToText(logicStatementInFix.substr(i,1));

		    numberOfSeqOperators = 1;
		    op_positions_in_subtext.erase(op_positions_in_subtext.begin(),op_positions_in_subtext.begin());
		    //treeRoot->printTree();
		    continue; 
		}




		

		// Add characters to growing text for current logic term
		(current->getData())->addToText(logicStatementInFix.substr(i,1));
		MSLOUT.stream() << "Current: "<<current->getData()->getText()<<"."<<endl;
		
	}
	
}
void LogicalParser::parse2(){

	// Current Tree pointer
	Tree<Predicate> *current;

	// Remove whitespace from back and front of statement
	logicStatementInFix = MslTools::toUpper(MslTools::trim(logicStatementInFix));
	//MSLOUT.stream() << "Input Statement: "<<logicStatementInFix;

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
					MSLOUT.stream() << "Operator("<<ss.str()<<") index is: "<<pos<<" start: "<<start<<endl;
				}
			}
			start      = logicStatementInFix.find_first_not_of(ss.str(), pos);
			loop++;
		}
	}

	/*
	  This doesn't work with the simple case 'all'
	if (validOpIndices.size() == 0){
		cerr << "ERROR 1234 LogicalParser::parse no valid OPEATORS in statement: "<<logicStatementInFix<<endl;
		exit(1234);
	}
	*/

	// Counters for our next loop
	int operatorCount         = 0;
	int openParenthesisCount  = 0;
	int closeParenthesisCount = 0;
	string correctFormat =  "";
	openParenthesisCount = closeParenthesisCount =0;

	if (debug){
		MSLOUT.stream() << "Fix this please: "<<logicStatementInFix<<endl;
	}


	// Check for problems with string before we do anything
	//  TO MANY BUGS... REQUIRE GOOD LOGICAL STATEMENT!


	for (uint i = 0 ; i < logicStatementInFix.size();i++){

		// Keep track of # open,close parens. At end of loop, we can die if counts are not equal.
		MSLOUT.stream() << "I["<<i<<"]: "<<logicStatementInFix[i]<<endl;
		if (logicStatementInFix[i] == '(') openParenthesisCount++;
		if (logicStatementInFix[i] == ')') closeParenthesisCount++;

		// Increment operator count if we have an operator index  (i).
		//    Also check for being in between ( ) . In that case don't add to count.
		if (validOpIndices[i] && openParenthesisCount == closeParenthesisCount+1) {

			operatorCount++; 
			if (debug){
				MSLOUT.stream() << "Index: "<<i<<" -> Increment OpCount("<<operatorCount<<"), have validOP and there are "<<openParenthesisCount<<" '(' and "<<closeParenthesisCount<<" ')'"<<endl;
			}
		}

		// PROBLEM: found 2 operators within '( )' term.  Add in '(' and ')'.
		if (operatorCount > 1){
			correctFormat = "(" + correctFormat + ") ";
			operatorCount--;
			continue;
		}


		correctFormat += logicStatementInFix[i];

		
	}

	// Exit if parenthesis count does not match (i.e. poorly formed selection string!)
	if (openParenthesisCount != closeParenthesisCount){
		cerr << "ERROR 3434 Parenthesis count are not equal : (open,close) ("<<openParenthesisCount<<","<<closeParenthesisCount<<") for statment: "<<logicStatementInFix<<endl;
		exit(3434);
	}


	//MSLOUT.stream() << "\tOpCount: "<<operatorCount<<" "<<openParenthesisCount<<" Corrected Statement: '"<<correctFormat<<"'"<<endl;

	// Reset logic statement
	logicStatementInFix = "(" + correctFormat + ")";



	//Special Case, no operators only have single '(' and ')'
	if (operatorCount == 0){
		std::replace(logicStatementInFix.begin(),logicStatementInFix.end(),'(',' ');
		std::replace(logicStatementInFix.begin(),logicStatementInFix.end(),')',' ');
		logicStatementInFix = "("+logicStatementInFix+")";
		
	}

	if (debug){
		MSLOUT.stream() << "LOGIC STATEMENT: "<<logicStatementInFix<<" "<<operatorCount<<endl;
	}

	// For each character
	for (uint i = 0 ; i < logicStatementInFix.size();i++){

		// Break loop if current pointer is NULL
		if (i > 0 && current == NULL) break;

		//MSLOUT.stream() << "Character: "<<logicStatementInFix[i]<<endl;

		// For each '(' create a new branch of the tree
		if (logicStatementInFix[i] == '('){

			// Create the tree root if tree root doesn't exist
			if (treeRoot == NULL) {
				treeRoot = new Tree<Predicate>();
				current = treeRoot;
				predicateList.push_back(current->getData());

				if (debug){
					MSLOUT.stream() << "tree root created: ->"<<current->getData()->getText()<<"<-"<<endl;
				}
			} else {

				// Add a subtree to the current tree
				current->addSubTree(new Tree<Predicate>(current));


				if (debug){
					MSLOUT.stream() << "Subtree added: ->"<<current->getData()->getText()<<"<-"<<endl;
				}
				// If no text yet, SubTree is 0.  Else Subtree is 1.
				if (current->getData()->getText() == ""){
					//MSLOUT.stream() << "\tBranch, SubTree1 (or Predicate1)"<<endl;
					current = (*current)[0]; // subTree 0 
				} else {
					//MSLOUT.stream() << "\tBranch, SubTree2 (or Predicate2) "<<endl;
					current = (*current)[1]; // subTree 1 
				}
				predicateList.push_back(current->getData());

			}
			continue;
		}
			

		// Ending logic term , parse the text into operand, operand, operator and move back up the tree
		if (logicStatementInFix[i] == ')'){
			if (debug){
				MSLOUT.stream() << "Moving back up the tree"<<endl;
			}
			(current->getData())->parseAsInFixText(validOperators);
			current = current->getParent();
			continue;
		}		


		// Add characters to growing text for current logic term
		(current->getData())->addToText(logicStatementInFix.substr(i,1));
		

	}


	//printLogicTree();
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

   treeRoot->printTree();
  //treeRoot->printTreeDFS(treeRoot);
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

	bool recursiveResult = false;
	//MSLOUT.stream() << "Current: "<<_node->getData()->toString()<<" subTrees: "<<numSubTrees<<" results: "<<cur->getNumResults()<<endl;
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
				double diff = _aLookupObject.getReal(keyword)-MslTools::toDouble(multipleValues[i]);

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
		MSLOUT.stream() << "Operand: "<<_predObj.getOperand(_operand)<< " was "<<result<<endl;
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

