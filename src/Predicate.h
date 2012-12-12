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

#ifndef PREDICATE_H
#define PREDICATE_H


#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>



namespace MSL { 
class Predicate {
	public:
		Predicate();
		~Predicate();

		void init();
		void setOperator(std::string _op) { op = _op;}
		std::string  getOperator() { return op;}

		void addOperand(std::string _op1) { operands.push_back(_op1); }
		std::string getOperand(int i) { return operands[i];} // needs error checking

		void addToText(std::string _str) { fullText += _str; }
		void setText(std::string _str) { fullText = _str; }
		std::string getText() { return fullText; }

		void parseAsPostFixText(std::vector<std::string> &_validOperators);
		void parseAsInFixText(std::vector<std::string> &_validOperators);

		bool getResult(int i) { return results[i]; } // needs error checking
		void addResult(bool _result) { results.push_back(_result); }

		int getNumResults() { return results.size(); }

		void clearResults() { results.clear(); }
		std::string toString();
		std::string trim(std::string str);
		std::vector<std::string> tokenizeWord(std::string _input, std::string _delimiter);

		friend std::ostream & operator<<(std::ostream &_os, Predicate &_pred)  {_os << _pred.toString(); return _os;};

	private:
		
		std::string op;           // an operator
		std::vector<std::string> operands;
		std::vector<bool> results;
		std::string fullText;     // text containing all three: op, operand1, operand2


};
}

#endif
