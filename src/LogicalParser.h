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

#ifndef LOGICALPARSER_H
#define LOGICALPARSER_H

#include <iostream>
#include <string>
#include <math.h>
#include <stack>

#include "Predicate.h"
#include "Tree.h"
#include "Selectable.h"


namespace MSL { 
class LogicalParser {
	public:
		LogicalParser(); 
		LogicalParser(std::string &_str) : logicStatementInFix(_str) {};
		~LogicalParser();

		inline std::string getLogicStatementInFix()            { return logicStatementInFix; }
		void setLogicStatementInFix(std::string _str);

		inline std::string getLogicStatementPostFix()          { return logicStatementPostFix; }
		inline void setLogicStatementPostFix(std::string _str) { logicStatementPostFix = _str; }

		// add operators
		inline void addOperator(std::string &_op)              { validOperators.push_back(_op); }
		inline std::vector<std::string>& getOperators()             { return validOperators; }

		// Create a postfix from infix and parse
		void parse();
		void parse_almost();
		void parse2();



		// Print the logic tree
		void printLogicTree();

		// Evaluate a selectable object here or in selection?
		bool eval(KeyLookup &_aLookupObject);

		// Temporary functions until tools class is created
		std::vector<std::string> tokenize(std::string input);		

		
		void setDebugFlag(bool _flag);
		bool getDebugFlag();

	private:

		// Evaluation sub-functions
		bool recursiveEval(KeyLookup &_aLookupObject, Tree<Predicate> *_node);
		int bailOutEarly(Predicate &predObj);
		void evalOperand(KeyLookup &_aLookupObject,Predicate &_predObj, int operand);

		void createPredicate(stack<string> &orderTree, string &currentStr, Tree<Predicate> **treeRoot, Tree<Predicate> **current );

		// Input std::strings
		std::string logicStatementInFix;
		std::string logicStatementPostFix;

		// Tree in flat form
		std::vector<Predicate *> predicateList;

		// Root of predicate tree (decision tree).
		Tree<Predicate>     *treeRoot;

		// Valid operators for this logical parser:
		std::vector<std::string> validOperators;

		bool debug;

};

// INLINE
inline void LogicalParser::setDebugFlag(bool _flag) { debug = _flag; }
inline bool LogicalParser::getDebugFlag() { return debug; }
}

#endif
