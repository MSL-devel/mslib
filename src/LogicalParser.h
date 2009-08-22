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

#ifndef LOGICALPARSER_H
#define LOGICALPARSER_H

#include <iostream>
#include <string>
#include <math.h>
using namespace std;


#include "Predicate.h"
#include "Tree.h"
#include "Selectable.h"


class LogicalParser {
	public:
		LogicalParser(); 
		LogicalParser(string &_str) : logicStatementInFix(_str) {};
		~LogicalParser();

		inline string getLogicStatementInFix()            { return logicStatementInFix; }
		void setLogicStatementInFix(string _str);

		inline string getLogicStatementPostFix()          { return logicStatementPostFix; }
		inline void setLogicStatementPostFix(string _str) { logicStatementPostFix = _str; }

		// add operators
		inline void addOperator(string &_op)              { validOperators.push_back(_op); }
		inline vector<string>& getOperators()             { return validOperators; }

		// Create a postfix from infix and parse
		void parse();

		// Print the logic tree
		void printLogicTree();

		// Evaluate a selectable object here or in selection?
		bool eval(KeyLookup &_aLookupObject);

		// Temporary functions until tools class is created
		vector<string> tokenize(string input);		

		
		void setDebugFlag(bool _flag);
		bool getDebugFlag();

	private:

		// Evaluation sub-functions
		bool recursiveEval(KeyLookup &_aLookupObject, Tree<Predicate> *_node);
		int bailOutEarly(Predicate &predObj);
		void evalOperand(KeyLookup &_aLookupObject,Predicate &_predObj, int operand);

		// Input strings
		string logicStatementInFix;
		string logicStatementPostFix;

		// Tree in flat form
		vector<Predicate *> predicateList;

		// Root of predicate tree (decision tree).
		Tree<Predicate>     *treeRoot;

		// Valid operators for this logical parser:
		vector<string> validOperators;

		bool debug;

};

// INLINE
inline void LogicalParser::setDebugFlag(bool _flag) { debug = _flag; }
inline bool LogicalParser::getDebugFlag() { return debug; }
#endif
