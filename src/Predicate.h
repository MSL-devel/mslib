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

#ifndef PREDICATE_H
#define PREDICATE_H


#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>


using namespace std;

class Predicate {
	public:
		Predicate();
		~Predicate();

		void init();
		void setOperator(string op);
		string  getOperator() { return op;}

		void addOperand(string _op1) { operands.push_back(_op1); }
		string getOperand(int i) { return operands[i];} // needs error checking

		void addToText(string _str) { fullText += _str; }
		string getText() { return fullText; }

		void parseAsPostFixText(vector<string> &_validOperators);
		void parseAsInFixText(vector<string> &_validOperators);

		bool getResult(int i) { return results[i]; } // needs error checking
		void addResult(bool _result) { results.push_back(_result); }

		int getNumResults() { return results.size(); }

		void clearResults() { results.clear(); }
		string toString();
		string trim(string str);
		vector<string> tokenizeWord(string _input, string _delimiter);

		friend ostream & operator<<(ostream &_os, Predicate &_pred)  {_os << _pred.toString(); return _os;};

	private:
		
		string op;           // an operator
		vector<string> operands;
		vector<bool> results;
		string fullText;     // text containing all three: op, operand1, operand2


};
#endif
