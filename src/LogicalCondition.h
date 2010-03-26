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

#ifndef LOGICALCONDITION_H
#define LOGICALCONDITION_H

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "MslTools.h"

namespace MSL {
	/**********************************************************************
	 *  This object is a logical tree (a set of possibly nested condtions).
	 *  It also takes values for the conditions and calculates and overall
	 *  state (true or false)
	 *
	 *  The structure of the tree is the following:
	 *   - a head node followed by a chain of linked nodes
	 *   - any node could have insead of a simple condition (i.e. "NAME CA")
	 *     a complex condition ("NAME CA AND RESI 37)", in which case the node
	 *     will point to a sub condition and its value will be the value of the
	 *     subcondition
	 *
	 *  Examples:
	 *     Simple:
	 *     "NAME CA OR NAME CB"
	 *
	 *         [NAME CA]->[OR NAME CB]
	 *
	 *     One sub-tree:
	 *     "NAME CA AND NOT (RESN ILE OR RESN LEU) AND CHAIN B"
	 *
	 *         [NAME CA]->[AND NOT *]->[AND CHAIN B]
	 *                             |
	 *                            [RESN ILE]->[OR RESN LEU]
	 *
	 *
	 *     Multiple sub-trees
	 *     "NAME CA AND ((RESI 37 AND CHAIN A) OR (RESI 41 AND CHAIN B)) AND NOT RESN ILE"
	 *
	 *         [NAME CA]->[AND *]->[AND NOT RESN ILE]
	 *                         |
	 *                        [*]->[AND *]
	 *                         |        |
	 *                         |       [RESI 41]->[AND CHAIN B]
	 *                         |
	 *                        [RESI 37]->[AND CHAIN A]
	 *
	 *   The logic is set via a string with setLogic(string _logic);
	 *   Each condition is returned with subsequent calls to getLogicalCondition()
	 *   For example, in the 2nd example, it will return "NAME CA" "RESN ILE" "RESN LEU"
	 *   and finally "CHAIN B".
	 *
	 *   The value (obtained by evaluating the object (i.e. Atom), are input back in the
	 *   LogicalCondition
	 *
	 *   When all questions are asked and answered logicComplete() will return true.
	 *
	 *   The overall value of the logic will be obtained with getOverallBooleanState()
	 *
	 *   Example:
	 *     setLogic("NAME CA AND NOT (RESN ILE OR RESN LEU) AND CHAIN B");
	 *     getLogicalCondition()  -> "NAME CA"
	 *     setLogicalConditionValue(true);
	 *     getLogicalCondition()  -> "RESN ILE"
	 *     setLogicalConditionValue(false);
	 *     getLogicalCondition()  -> "RESN LEU"
	 *     setLogicalConditionValue(false);
	 *     getLogicalCondition()  -> "CHAIN B"
	 *     setLogicalConditionValue(true);
	 *     logicComplete() -> true
	 *     getOverallBooleanState() -> true
	 *    TODO
	 *
	 *  1 This object supports a mechanism to
	 *    skip unnecessary checks but it is disabled 
	 *    (shortCutLogic_flag set to false), it needs 
	 *    more debugging
	 *
	 *  2 Also, after that is done, add support for swapping
	 *    the order or "heavy" and light "branches" so
	 *    that the first get potentially skipped
	 * 
	 *  3 Add the support for a NOT operator that
	 *    inverts the logic
	 *
	 *  4 The AND/NOT/XOR should not be looked up as
	 *    strings but converted to int (enum) to
	 *    speedup matching
	 **********************************************************************/

	class LogicalCondition {
		public:
			LogicalCondition();
			~LogicalCondition();

			bool setLogic(std::string _logic);
			void setLogicalOperator(std::string _operator);

			void operator=(const LogicalCondition & _cond);

			std::vector<std::string> getLogicalCondition();

			void setLogicalConditionValue(bool _value, bool _phony=false);
			bool getOverallBooleanState();
			bool logicComplete() const;

			std::string printLogicalConditions() const;
			std::string printLogicalConditionsWithValues() const;

			void restartQuery(); // forget the previously inputed logical values and restart querying the logic
			void reset(); // remove the logic completely

		private:
			void setup();
			void copy(const LogicalCondition & _cond);
			void deletePointers();
			std::vector<std::string> parseBrackets(std::string _input);
			bool setLogic(std::vector<std::string> _logic);
			bool splitLogic(std::vector<std::string> & _tokenizedLogic, std::vector<std::string> & _leftSideTokens, std::vector<std::string> & _rightSideTokens, std::string & _logicalOperator, bool & _isSimpleLeft, bool & _isSimpleRight);
			void setParent(LogicalCondition * _pParentCondition);
			void setPrev(LogicalCondition * _pPrevCondition);
			void setHead(LogicalCondition * _pHeadCondition);
			void setBranchHead(LogicalCondition * _pBranchHeadCondition);
			void setTreeDone();
			void clearExternalBrakets(std::vector<std::string> & _tokens);
			void propagateLogicForward();
			void headCalculatesOverallBooleanState();
			void invertCondition(bool _invert);

			std::vector<std::string> conditionTokens;
			std::string logicalOperator;
			unsigned int operatorCode;
			LogicalCondition * pNextCondition;
			LogicalCondition * pSubCondition;
			LogicalCondition * pCurrentCondition;
			LogicalCondition * pParentCondition;
			LogicalCondition * pPrevCondition;
			LogicalCondition * pHeadCondition;
			LogicalCondition * pBranchHeadCondition;
			bool localValue;
			bool localValueIsPhony;
			bool overallValue;
			bool localDone;
			bool treeDone;
		//	bool computed;
			bool NOT_flag; // if NOT was given the object will have to invert the logic
			bool shortCutLogic_flag;
			std::vector<std::string> validOperators;
	};

	inline bool LogicalCondition::logicComplete() const {
		return localDone && treeDone;
	}
	inline bool LogicalCondition::setLogic(std::string _logic) {
		_logic = MslTools::trim(MslTools::toUpper(_logic));
		std::vector<std::string> tokens = parseBrackets(_logic);
		return setLogic(tokens);
	}
	inline void LogicalCondition::setParent(LogicalCondition * _pParentCondition) {
		pParentCondition = _pParentCondition;
	}
	inline void LogicalCondition::setPrev(LogicalCondition * _pPrevCondition) {
		pPrevCondition = _pPrevCondition;
	}
	inline void LogicalCondition::setHead(LogicalCondition * _pHeadCondition) {
		pHeadCondition = _pHeadCondition;
	}
	inline void LogicalCondition::setBranchHead(LogicalCondition * _pBranchHeadCondition) {
		pBranchHeadCondition = _pBranchHeadCondition;
	}
	inline void LogicalCondition::setLogicalOperator(std::string _operator) {
		logicalOperator = _operator;
		if (_operator == "AND") {
			operatorCode = 1;
		} else if (_operator == "OR") {
			operatorCode = 2;
		} else if (_operator == "XOR") {
			operatorCode = 3;
		} else {
			std::cerr << "WARNING 23829: invalid operator " << _operator << " in inline void LogicalCondition::setLogicalOperator(std::string _operator)" << std::endl;
			operatorCode = 0;
		}
	}
	inline void LogicalCondition::restartQuery() {
		pCurrentCondition = this;
		localValue = true;
		localValueIsPhony = true;
		localDone = false;
		if (pNextCondition != NULL) {
			treeDone = false;
			pNextCondition->restartQuery();
		} else {
			treeDone = true;
		}
		if (pSubCondition != NULL) {
			pSubCondition->restartQuery();
		}
	//	computed = false;
	}
			
};
#endif

