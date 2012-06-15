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
	 *
	 *
	 *  NOTE: the priority of the operators follows the same rules of C++ and
	 *        other languages (as well as PyMOL).  AND is evaluated first, then
	 *        XOR and then OR.
	 *
	 *        For example:  "NAME CA OR RESI 7 AND CHAIN B XOR RESN ILE" is the same
	 *        of "NAME CA OR ((RESI 7 AND CHAIN B) XOR RESN ILE)"
	 *
	 *    TODO
	 *
	 *  1 Add the support for a mechanism to  skip unnecessary checks. 
	 *   For example "NAME CA AND RESI 7" should not ask for the 
	 *   residue number is the atom name is not CA.  "NAME CA OR RESI 7"
	 *   should not ask for the residue number is the name is CA.
	 *
	 *  2 Also, after that is done, add support for swapping
	 *    the order or "heavy" and light "branches" so
	 *    that the first get potentially skipped
	 * 
	 **********************************************************************/

	class LogicalCondition {
		public:
			LogicalCondition();
			~LogicalCondition();

			bool setLogic(std::string _logic);
			void setLogicalOperator(std::string _operator);

			void operator=(const LogicalCondition & _cond);

			std::vector<std::string> getLogicalCondition();

			//void setLogicalConditionValue(bool _value, bool _phony=false);
			void setLogicalConditionValue(bool _value);
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
//<<<<<<< .mine
			void propagateLogicForward(unsigned int _stage);
//=======
//			void propagateLogicForward();
//			void headCalculatesOverallBooleanState();
//>>>>>>> .r450
			void invertCondition(bool _invert);
			bool getOverallValue() const;
			void setOverallValue(bool _val);

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
			bool overallValue;
			bool localDone;
			bool treeDone;
//<<<<<<< .mine
			bool computed;
			bool NOT_flag; // if NOT was given the object will have to invert the logic
//=======
//		//	bool computed;
//			bool NOT_flag; // if NOT was given the object will have to invert the logic
//			bool shortCutLogic_flag;
//>>>>>>> .r450
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
//<<<<<<< .mine
		if (_operator == "AND") {
			operatorCode = 1;
		} else if (_operator == "XOR") {
			operatorCode = 2;
		} else if (_operator == "OR") {
			operatorCode = 3;
		} else {
			std::cerr << "WARNING 23829: invalid operator " << _operator << " in inline void LogicalCondition::setLogicalOperator(std::string _operator)" << std::endl;
			operatorCode = 0;
		}
//=======
//		if (_operator == "AND") {
//			operatorCode = 1;
//		} else if (_operator == "OR") {
//			operatorCode = 2;
//		} else if (_operator == "XOR") {
//			operatorCode = 3;
//		} else {
//			std::cerr << "WARNING 23829: invalid operator " << _operator << " in inline void LogicalCondition::setLogicalOperator(std::string _operator)" << std::endl;
//			operatorCode = 0;
//		}
//>>>>>>> .r450
	}
	inline void LogicalCondition::restartQuery() {
		pCurrentCondition = this;
		localValue = true;
		//localValueIsPhony = true;
		localDone = false;
		computed = false;
		if (pNextCondition != NULL) {
			treeDone = false;
			pNextCondition->restartQuery();
		} else {
			treeDone = true;
		}
		if (pSubCondition != NULL) {
			pSubCondition->restartQuery();
		}
	}

			
};
#endif

