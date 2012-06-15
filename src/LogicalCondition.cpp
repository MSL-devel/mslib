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

#include "LogicalCondition.h"

using namespace std;
using namespace MSL;

LogicalCondition::LogicalCondition() {
	setup();
}

LogicalCondition::~LogicalCondition() {
	deletePointers();
}

void LogicalCondition::deletePointers() {
	if (pNextCondition != NULL) {
		pNextCondition->deletePointers();
		delete pNextCondition;
		pNextCondition = NULL;
	}
	if (pSubCondition != NULL) {
		pSubCondition->deletePointers();
		delete pSubCondition;
		pSubCondition = NULL;
	}
}

void LogicalCondition::setup() {
	pNextCondition = NULL;
	pSubCondition = NULL;
	pCurrentCondition = this;
	pParentCondition = NULL;
	pPrevCondition = NULL;
	pHeadCondition = this;
	pBranchHeadCondition = this;
	logicalOperator = "";
	operatorCode = 0;
	localValue = true;
	overallValue = true;
	localDone = false;
	treeDone = true;
	computed = false;
	validOperators.push_back("AND");
	validOperators.push_back("OR");
	validOperators.push_back("XOR");
	NOT_flag = false;
}

void LogicalCondition::copy(const LogicalCondition & _cond) {
	reset();

	if (_cond.pSubCondition != NULL) {
		pSubCondition = new LogicalCondition(*_cond.pSubCondition);
		pSubCondition->setParent(this);
		pSubCondition->setHead(pHeadCondition);
		pSubCondition->setBranchHead(pBranchHeadCondition);
	}
	if (_cond.pNextCondition != NULL) {
		pNextCondition = new LogicalCondition(*_cond.pNextCondition);
		pNextCondition->setPrev(this);
		pNextCondition->setHead(pHeadCondition);
		pSubCondition->setBranchHead(pBranchHeadCondition);
	}
	logicalOperator = _cond.logicalOperator;
	operatorCode = _cond.operatorCode;
	localValue = _cond.localValue;
	overallValue = _cond.overallValue;
	localDone = _cond.localDone;
	treeDone = _cond.treeDone;
	computed = _cond.computed;
	NOT_flag = _cond.NOT_flag;

}

void LogicalCondition::reset() {
	deletePointers();
	pCurrentCondition = this;
	logicalOperator = "";
	operatorCode = 0;
	localValue = true;
	overallValue = true;
	localDone = false;
	treeDone = true;
	computed = false;
	NOT_flag = false;
}

void LogicalCondition::operator=(const LogicalCondition & _cond) {
	copy(_cond);
}

vector<string> LogicalCondition::parseBrackets(string _input) {
	/*********************************************************
	 *  Turns
	 *    (RESI 7 AND (NAME CA OR NAME CB)) OR RESN ALA
	 *
	 *  into the following vector
	 *
	 *   0 (
	 *   1 RESI
	 *   2 7
	 *   3 AND
	 *   4 (
	 *   5 NAME
	 *   6 CA
	 *   7 OR
	 *   8 NAME
	 *   9 CB
	 *  10 )
	 *  11 )
	 *  12 OR
	 *  13 RESN
	 *  14 ALA
	 * 
	 *  Note: if the string starts and end with unnecessary ( ) these will be 
	 *        removed 
	 *        (RESI 7 AND (NAME CA OR CB)) -> RESI 7 AND (NAME CA OR CB)
	 *********************************************************/
	// Remove whitespace from back and front of statement
	_input =MslTools::toUpper(MslTools::trim(_input));
	vector<string> tokens = MslTools::tokenizeAndTrim(_input);
	for (unsigned int i=0; i<tokens.size(); i++) {
		size_t opos = tokens[i].find("(",0);
		size_t cpos = tokens[i].find(")",0);
		string braketType = "";
		bool found = false;
		size_t pos = 0;
		if (opos != std::string::npos) {
			braketType = "(";
			found = true;
			pos = opos;
		} else if (cpos != std::string::npos) {
			braketType = ")";
			found = true;
			pos = cpos;
		}
		if (found) {
			if (tokens[i] == braketType) {
				continue;
			}
			string before = tokens[i].substr(0, pos);
			string after = tokens[i].substr(pos+1, tokens[i].size() - pos - 1);
			tokens[i] = braketType;
			if (before != "") {
				tokens.insert(tokens.begin()+i, before);
				i++;
			}
			if (after != "") {
				tokens.insert(tokens.begin()+i+1, after);
			}
		}
	}
	clearExternalBrakets(tokens);
	return tokens;
}

void LogicalCondition::clearExternalBrakets(vector<string> & _tokens) {
	// remove useless brakets that wrap around the entire statement
	unsigned int counter = 1;
	while (_tokens.size() >=2 && _tokens[0] == "(" && _tokens.back() == ")") {
		bool erase = true;
		for (unsigned int i=1; i<_tokens.size()-1; i++) {
			if (_tokens[i] == "(") {
				counter++;
			} else if (_tokens[i] == ")") {
				counter--;
				if (counter == 0) {
					erase = false;
					break;
				}
			}
		}
		if (erase) {
			_tokens.erase(_tokens.end()-1);
			_tokens.erase(_tokens.begin());
		} else {
			break;
		}
	}
}

bool LogicalCondition::splitLogic(vector<string> & _tokenizedLogic, vector<string> & _leftSideTokens, vector<string> & _rightSideTokens, string & _logicalOperator, bool & _isSimpleLeft, bool & _isSimpleRight) {
	/**************************************************************
	 *  Split the tokenizedLogic in two parts, finding the main opertor,
	 *  and checks if the other parts have operators in it, in
	 *  which case they need to be further processed (they are not
	 *  simple
	 **************************************************************/

	_isSimpleLeft = true;
	_isSimpleRight = true;

	bool foundSplicingOperator = false;
	unsigned int braketLevel = 0;
	for (int i=0; i<_tokenizedLogic.size(); i++) {
		bool foundOperator = false;
		for (unsigned j=0; j<validOperators.size(); j++) {
			if (_tokenizedLogic[i] == validOperators[j]) {
			       // it is an operator
				foundOperator = true;
				break;
			}
		}
		if (foundOperator) {
			if (braketLevel == 0 && !foundSplicingOperator) {
				// it is the main operator
				foundSplicingOperator = true;
				_logicalOperator = _tokenizedLogic[i];
				continue;
			}
			if (foundSplicingOperator) {
				// found an extra operator on the right side, pattern is not simple
				_isSimpleRight = false;
			} else {
				// found an extra operator on the left side, pattern is not simple
				_isSimpleLeft = false;
			}
		}
		// check if it is a braket and adjust the braket level
		if (_tokenizedLogic[i] == "(") {
			// make sure the braket comes after an operator (AND OR XOR)
			bool afterValidOperator = false;
			if (i != 0) {
				for (unsigned j=0; j<validOperators.size(); j++) {
					if (_tokenizedLogic[i-1] == validOperators[j]) {
						afterValidOperator = true;
						break;
					}
					if (_tokenizedLogic[i-1] == "NOT" || _tokenizedLogic[i-1] == "(") {
						// it could also be NOT or another braket
						afterValidOperator = true;
					}
				}
				if (!afterValidOperator) {
					cerr << "WARNING 48226: brakets not preceed by operator in bool LogicalCondition::splitLogic(vector<string> & _tokenizedLogic, vector<string> & _leftSideTokens, vector<string> & _rightSideTokens, string & _logicalOperator, bool & _isSimpleLeft, bool & _isSimpleRight)" << endl;
					return false;
				}
			}
			braketLevel++;
				
		} else if (_tokenizedLogic[i] == ")") {
			if (braketLevel == 0) {
				cerr << "WARNING 48231: Invalit brakets in string in bool LogicalCondition::splitLogic(vector<string> & _tokenizedLogic, vector<string> & _leftSideTokens, vector<string> & _rightSideTokens, string & _logicalOperator, bool & _isSimpleLeft, bool & _isSimpleRight)" << endl;
				return false;
			}
			// make sure the braket is followed by an operator (AND OR XOR but not NOT)
			bool beforeValidOperator = false;
			if (i != _tokenizedLogic.size()-1) {
				for (unsigned j=0; j<validOperators.size(); j++) {
					if (_tokenizedLogic[i+1] == validOperators[j]) {
						beforeValidOperator = true;
						break;
					}
				}
				if (_tokenizedLogic[i+1] == ")") {
					// it could also be another braket
					beforeValidOperator = true;
				}
				if (!beforeValidOperator) {
					cerr << "WARNING 48226: brakets not followed by operator in bool LogicalCondition::splitLogic(vector<string> & _tokenizedLogic, vector<string> & _leftSideTokens, vector<string> & _rightSideTokens, string & _logicalOperator, bool & _isSimpleLeft, bool & _isSimpleRight)" << endl;
					return false;
				}
			}
			braketLevel--;
		}
		if (foundSplicingOperator) {
			// add to the right side tokenizedLogic
			_rightSideTokens.push_back(_tokenizedLogic[i]);
		} else {
			// add to the left side tokenizedLogic
			_leftSideTokens.push_back(_tokenizedLogic[i]);
		}
	}
	//clearExternalBrakets(_leftSideTokens);
	//clearExternalBrakets(_rightSideTokens);
	if (braketLevel != 0 || (foundSplicingOperator && (_leftSideTokens.size() == 0 || _rightSideTokens.size() == 0))) {
		// if we have unmatched brakets or an operator both next to an emtpy condition ("AND RESI 37" or "NAME CA AND")
		return false;
	}
	return true;
}

bool LogicalCondition::setLogic(vector<string> _logicTokens) {
	/*******************************************************
	 *   This function creates a linked chains of conditions.
	 *   If a condition is complex (brakets) then these are
	 *   processed as sub-conditions (themselves linked
	 *   chain of conditions, that can also be branched)
	 *
	 *   "RESI 37 AND (NAME CB OR (RESN VAL AND NAME CG1) AND CHAIN B"
	 *
	 *    RESI 37 --AND--> * --AND--> CHAIN B
	 *                     |
	 *                     NAME CB --OR--> *
	 *                                     |
	 *                                     RESN VAL --AND--> NAME CG1
	 *******************************************************/
	// cleanup
	reset();
	
	vector<string> leftSideTokens;
	vector<string> rightSideTokens;
	string logicalOp;
	bool isSimpleLeft;
	bool isSimpleRight;
	if (!splitLogic(_logicTokens, leftSideTokens, rightSideTokens, logicalOp, isSimpleLeft, isSimpleRight)) {
		// failed, cleanup
		deletePointers();
		localDone = true;
		treeDone = true;
		conditionTokens.clear();
		return false;
	}
	if (leftSideTokens[0] == "NOT") {
		NOT_flag = true;
		leftSideTokens.erase(leftSideTokens.begin());
	}
	if (isSimpleLeft) {
		// simple conditions, it has its tokens (i.e. "NAME" "CA")
		conditionTokens = leftSideTokens;
	} else {
		clearExternalBrakets(leftSideTokens);
		// complex condition (i.e. "NAME" "CA" "AND" "RESI" "37"), it will be parsed by a subcondition
		pSubCondition = new LogicalCondition;
		pSubCondition->setParent(this);
		pSubCondition->setHead(pHeadCondition);
		pSubCondition->setBranchHead(pSubCondition);
		if (!pSubCondition->setLogic(leftSideTokens)) {
			return false;
		}
	}
	if (rightSideTokens.size() > 0) {
		// we have a logical condition after this
		pNextCondition = new LogicalCondition;
		pNextCondition->setPrev(this);
		pNextCondition->setHead(pHeadCondition);
		pNextCondition->setBranchHead(pBranchHeadCondition);
		if (!pNextCondition->setLogic(rightSideTokens)) {
			return false;
		}
		pNextCondition->setLogicalOperator(logicalOp); // AND/OR/XOR
		treeDone = false;
	}
	/*******************************************************
	 *  TODO:
	 *  Here code could be added to move around the LogicalConditions   
	 *  in the linked chain so that the least heavy are executed
	 *  first, for speed
	 *******************************************************/
	return true;
}

vector<string> LogicalCondition::getLogicalCondition() {
	/*******************************************************
	 *  This object returns a sequence of logical questions
	 *  (tokenized, i.e. ("NAME" "CA"), going through the
	 *  subconditions and the chain of conditions, keeping
	 *  track of what to return next.
	 * 
	 *  The pCurrentCondition pointer keeps track of where the logical
	 *  question was coming from, so that a subsequent call to
	 *  setLogicalConditionValue(true/false) can be placed in the
	 *  right place.
	 * 
	 *  We can find out if we are done by calling logicComplete()
	 *
	 *  At the end we will be able to get a true/false for 
	 *  the overall logic with a call to getOverallBooleanState()
	 *******************************************************/
	vector<string> out;
	if (localDone && treeDone) {
		// we are already done
		return out;
	}
	if (!localDone) {
		// the local condition is not done
		if (pSubCondition != NULL) {
			// this condition has a subcondition, we'll pass the ball
			pCurrentCondition = pSubCondition;
			out = pSubCondition->getLogicalCondition();
			if (pSubCondition->logicComplete()) {
				// sub condition is done
				localDone = true;
			}
			return out;
		} else {
			if (pNextCondition == NULL) {
				// end of a chain, tell the head that the chain is done
				setTreeDone();
			}
			pCurrentCondition = this;
			localDone = true;
			// return the tokens, we are done with this condition
			return conditionTokens;
		}
	} else {
		// this condition is done, let's call the next condition (if there is one)
		if (pNextCondition != NULL) {
			if (!pNextCondition->logicComplete()) {
				pCurrentCondition = pNextCondition;
				return pNextCondition->getLogicalCondition();
			}
		}
	}
	// something went wrong
	localDone = true;
	treeDone = true;
	return out;
}


//void LogicalCondition::setLogicalConditionValue(bool _value, bool _phony) {
void LogicalCondition::setLogicalConditionValue(bool _value) {
	/************************************************************
	 *  After asking a question with getLogicalCondition() we input the
	 *  answer in the LogicalCondition object, it will be used to
	 *  compute the overall state of the complex condition
	 ************************************************************/
	if (pCurrentCondition == this) {
		localValue = _value ^ NOT_flag;
		overallValue = localValue;
	} else {
		pCurrentCondition->setLogicalConditionValue(_value);
	}
}

void LogicalCondition::setTreeDone() {
	/******************************************
	 *   Propagate backward in the chain to let
	 *   the beginning condition know that the
	 *   chain is done
	 ******************************************/
	if (pPrevCondition != NULL) {
		pPrevCondition->setTreeDone();
	} else {
		treeDone = true;
		if (pParentCondition != NULL) {
			pParentCondition->localDone = true;
			if (pParentCondition->pNextCondition == NULL) {
				/******************************************
				 * it this is a sub-condition and the pParentCondition is 
				 * an end of chain, let's initiate the call in the
				 * pParentCondition 
				 ******************************************/
				pParentCondition->setTreeDone();
			}
		}
	}
}

string LogicalCondition::printLogicalConditions() const {
	ostringstream ss;
	if (logicalOperator != "") {
		ss << " " << logicalOperator << " ";
	}
	if (NOT_flag) {
		ss << "NOT ";
	}
	if (pSubCondition != NULL) {
		ss << "(" << pSubCondition->printLogicalConditions() << ")";
	} else {
		ss << MslTools::joinLines(conditionTokens);
	}
	if (pNextCondition != NULL) {
		ss << pNextCondition->printLogicalConditions();
	}
	return ss.str();
}

string LogicalCondition::printLogicalConditionsWithValues() const {
	ostringstream ss;
	if (logicalOperator != "") {
		ss << " " << logicalOperator << " ";
	}
	if (NOT_flag) {
		ss << "NOT ";
	}
	if (pSubCondition != NULL) {
		ss << "(" << pSubCondition->printLogicalConditionsWithValues() << ") [" << localValue << "]";
	} else {
		ss << MslTools::joinLines(conditionTokens) << " [" << localValue << "]";
	}
	if (pNextCondition != NULL) {
		ss << pNextCondition->printLogicalConditionsWithValues();
	}
	if (pPrevCondition == NULL && pParentCondition == NULL) {
		ss << " = ";
		if (logicComplete()) {
			ss << "[" << overallValue << "]";
		} else {
			ss << "[?]" << endl;
		}
	}

	return ss.str();
}


bool LogicalCondition::getOverallBooleanState() {
	if (!logicComplete()) {
		cerr << "WARNING: requested overall boolean status of the complex LogicalCondition but the logic is not complete, in bool LogicalCondition::getOverallBooleanState()" << endl;
	}
	if (pSubCondition != NULL) {
		// resolve the overall value of a sub-condition branch and set it as local
		localValue = pSubCondition->getOverallBooleanState() ^ NOT_flag;
		overallValue = localValue;
	}
	if (pNextCondition != NULL) {
		// we have a next condition, propagate the logic forward.  This goes in 3 stages
		// to respect the priority of the AND -> XOR -> OR operators.
		// the number 1 indicates that it will resolve the AND operations
		pNextCondition->propagateLogicForward(1);
	}

	return overallValue;
}

void LogicalCondition::propagateLogicForward(unsigned int _stage) {
	/*************************************************************
	 * This subroutine is called in 3 stages, that resolve the 
	 * AND, XOR and OR relationship in their order of priority, which is
	 * kept by the _stage variable
	 *  1 AND
	 *  2 XOR
	 *  3 OR
	 * If the current LogicalCondition is computed or the stage does not
	 * correspond to the right operator, the cycle skips to the
	 * next one.  Otherwise, the value of the condition and that of the first
	 * previous condition that has not been computed yet are combined, that value
	 * is stored in the previous condition and this is marked as computed. 
	 *
	 * Note: the code assume that pPrevCondition is defined
	 *************************************************************/
	if (_stage == 1 && pSubCondition != NULL) {
		// resolve the overall value of a sub-condition branch and set it as local
		localValue = pSubCondition->getOverallBooleanState() ^ NOT_flag;
		overallValue = localValue;
	}
	if (!computed) {
		if (operatorCode == 1 && _stage == 1) {
			// AND
			bool newValue = pPrevCondition->getOverallValue() & overallValue;
			pPrevCondition->setOverallValue(newValue);
			computed = true;
		} else if (operatorCode == 2 && _stage == 2) {
			// XOR
			bool newValue = pPrevCondition->getOverallValue() ^ overallValue;
			pPrevCondition->setOverallValue(newValue);
			computed = true;
		} else if (operatorCode == 3 && _stage == 3) {
			// OR
			bool newValue = pPrevCondition->getOverallValue() | overallValue;
			pPrevCondition->setOverallValue(newValue);
			computed = true;
		}
	}
	if (pNextCondition != NULL) {
		pNextCondition->propagateLogicForward(_stage);
	} else {
		// we have the overall results at the end of the chain: send it back to the head LogicalCondition
		if (_stage < 3) {
			_stage++;
			pBranchHeadCondition->pNextCondition->propagateLogicForward(_stage);
		} else {
			// we are done
			if (pBranchHeadCondition->pParentCondition != NULL) {
				// propagate up the result to the branch parent
				pBranchHeadCondition->pParentCondition->overallValue = overallValue;
			}
		}
	}
}

bool LogicalCondition::getOverallValue() const {
	// recursively get the overallValue from
	// the first non-computed previous condition
	// in the chain
	if (!computed) {
		return overallValue;
	} else {
		// assuming pPrevCondition is not NULL
		return pPrevCondition->getOverallValue();
	}
}
void LogicalCondition::setOverallValue(bool _val) {
	// recursively set the overallValue from
	// the first non-computed previous condition
	// in the chain
	if (!computed) {
		overallValue = _val;
	} else {
		// assuming pPrevCondition is not NULL
		pPrevCondition->setOverallValue(_val);
	}
}

