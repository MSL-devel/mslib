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

#include <iostream>

#include "MslTools.h"

using namespace std;

using namespace MSL;


int main() {

	cout << "Test string trim(const string & _str);" << endl;

	string a = "     pizza   pizza pizza                  ";
	cout << "a = >" << a << "<, trim(a) = >" << MslTools::trim(a) << "<" << endl;

	a = "pizza   pizza pizza                  ";
	cout << "a = >" << a << "<, trim(a) = >" << MslTools::trim(a) << "<" << endl;

	a = "     pizza   pizza pizza";
	cout << "a = >" << a << "<, trim(a) = >" << MslTools::trim(a) << "<" << endl;

	a = "a     pizza   pizza pizza                  a";
	cout << "a = >" << a << "<, trim(a) = >" << MslTools::trim(a) << "<" << endl;

	a = "                           ";
	cout << "a = >" << a << "<, trim(a) = >" << MslTools::trim(a) << "<" << endl;

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test double toDouble(const string & _string, const string & _msg=(string)"");" << endl;

	a = "8.521454365555";
	double b = 0.0;
	try {
		b = MslTools::toDouble(a);
		cout << "a = " << a << ", b = " << b << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an double" << endl;
	}

	a = "-9";
	try {
		b = MslTools::toDouble(a);
		cout << "a = " << a << ", b = " << b << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an double" << endl;
	}

	a = "-.5898842";
	try {
		b = MslTools::toDouble(a);
		cout << "a = " << a << ", b = " << b << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an double" << endl;
	}

	a = "pizza";
	try {
		b = MslTools::toDouble(a);
		cout << "a = " << a << ", b = " << b << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an double" << endl;
	}

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test int toInt(const string & _string, const string & _msg=(string)"")" << endl;

	a = "8.521454365555";
	int c = 0;
	try {
		c = MslTools::toInt(a);
		cout << "a = " << a << ", c = " << c << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an int" << endl;
	}

	a = "-9";
	try {
		c = MslTools::toInt(a);
		cout << "a = " << a << ", c = " << c << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an int" << endl;
	}

	a = "-.5898842";
	try {
		c = MslTools::toInt(a);
		cout << "a = " << a << ", c = " << c << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an int" << endl;
	}

	a = "pizza";
	try {
		c = MslTools::toInt(a);
		cout << "a = " << a << ", c = " << c << endl;
	} catch (exception &e) {
		cout << "a = " << a << ", ERROR: is not an int" << endl;
	}


	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string toString(const int & _i, const string & _msg=(string)"")" << endl;

	c = 765;
	a = MslTools::intToString(c);
	cout << "c = " << c << ", a = " << a << endl;

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test void splitIntAndString(const string _input, int & _intResult, string & _stringResult)" << endl;

	a = "85C";	
	string d;
	MslTools::splitIntAndString(a, c, d);
	cout << "a = " << a << ", c = " << c << ", d = " << d << endl;

	a = "  85 C ";	
	MslTools::splitIntAndString(a, c, d);
	cout << "a = " << a << ", c = " << c << ", d = " << d << endl;

	a = "85";	
	MslTools::splitIntAndString(a, c, d);
	cout << "a = " << a << ", c = " << c << ", d = " << d << endl;

	a = "C";	
	MslTools::splitIntAndString(a, c, d);
	cout << "a = " << a << ", c = " << c << ", d = " << d << endl;


	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string toUpper(const string & _input)" << endl;
	cout << endl;

	a = "This is a mixed cased Sentence";	
	d = MslTools::toUpper(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "tHIS ALSO IS A MIXED CASED sENTENCE";	
	d = MslTools::toUpper(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "this also is a lower cased sentence";	
	d = MslTools::toUpper(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "THIS ALSO IS A UPPER CASED SENTENCE";	
	d = MslTools::toUpper(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "This sentence has number 1 2 3 4 and punctuation # $ ^ &";	
	d = MslTools::toUpper(a);
	cout << "a = " << a << ", d = " << d << endl;


	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string toLower(const string & _input)" << endl;
	cout << endl;

	a = "This is a mixed cased Sentence";	
	d = MslTools::toLower(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "tHIS ALSO IS A MIXED CASED sENTENCE";	
	d = MslTools::toLower(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "this also is a lower cased sentence";	
	d = MslTools::toLower(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "THIS ALSO IS A UPPER CASED SENTENCE";	
	d = MslTools::toLower(a);
	cout << "a = " << a << ", d = " << d << endl;

	a = "This sentence has number 1 2 3 4 and punctuation # $ ^ &";	
	d = MslTools::toLower(a);
	cout << "a = " << a << ", d = " << d << endl;

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test vector<string> tokenize(const string & _input, const string & _delimiter)" << endl;

	cout << endl;
	a = "This is a random sentence";
	d = " ";
	cout << "Tokenize >" << a << "< with >" << d << "<" << endl;
	vector<string> e = MslTools::tokenize(a, d);
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	a = "This is a random sentence";
	d = "s";
	cout << "Tokenize >" << a << "< with >" << d << "<" << endl;
	e = MslTools::tokenize(a, d);
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	a = "This is a random sentence";
	d = "en";
	cout << "Tokenize >" << a << "< with >" << d << "<" << endl;
	e = MslTools::tokenize(a, d);
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string joinLines(const vector<string> & _input)" << endl;

	cout << endl;
	e.clear();
	e.push_back("Line 1");
	e.push_back("Line 2");
	e.push_back("Line 3");
	e.push_back("Line 4");
	a = MslTools::joinLines(e);
	cout << "Join:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << a << endl;


	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string joinLines(const string & _spacer, const vector<string> & _input)" << endl;

	cout << endl;
	d = "-*-";
	a = MslTools::joinLines(e, d);
	cout << "Join with space \"" << d << "\":" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << a << endl;


	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
        /*
	cout << "Test vector<string> joinBackslashedLines(const vector<string> & _input)" << endl;

	cout << endl;
	e.clear();
	e.push_back("Line 1");
	e.push_back("Line 2\\");
	e.push_back("Line 3     \\");
	e.push_back("Line 4");
	e.push_back("Line 5  \\ ");
	e.push_back("Line 6");
        */
        vector<string> f;
        /*
        vector<string> f = MslTools::joinBackslashedLines(e);
	cout << "Join:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << "Result:" << endl;
	for (vector<string>::iterator k=f.begin(); k!=f.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	e.clear();
	e.push_back("Line 1\\");
	e.push_back("Line 2\\");
	e.push_back("Line 3\\");
	e.push_back("Line 4\\");
	e.push_back("Line 5\\");
	e.push_back("Line 6\\");
	f = MslTools::joinBackslashedLines(e);
	cout << "Join:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << "Result:" << endl;
	for (vector<string>::iterator k=f.begin(); k!=f.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	e.clear();
	e.push_back("Line 1");
	e.push_back("Line 2");
	e.push_back("Line 3");
	e.push_back("Line 4");
	e.push_back("Line 5");
	e.push_back("Line 6\\");
	f = MslTools::joinBackslashedLines(e);
	cout << "Join:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << "Result:" << endl;
	for (vector<string>::iterator k=f.begin(); k!=f.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
         */
	cout << "Test vector<string> removeEmptyLines(const vector<string> & _input)" << endl;

	cout << endl;
	e.clear();
	e.push_back("Line 1");
	e.push_back("Line 2");
	e.push_back("");
	e.push_back("");
	e.push_back("Line 3");
	e.push_back("Line 4");
	e.push_back("Line 5");
	e.push_back("             ");
	e.push_back("Line 6");
	e.push_back("");
	f = MslTools::removeEmptyLines(e);
	cout << "Join:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << "Result:" << endl;
	for (vector<string>::iterator k=f.begin(); k!=f.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test string uncomment(const string & _input, const string & _commentString=\"#\")" << endl;

	a = "This is a value  # this is a comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;

	a = "This sentence has a \"# within quotes\" and it is not a comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;

	a = "This sentence has it \"outside\" the quotes and it is a comment # here is the comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;

	a = "This sentence has a \'# within single quotes\' and it is not a comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;

	a = "This sentence has it \'outside\' single quotes and it is a comment # here is the comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;

	a = "# everything is a comment";
	d = MslTools::uncomment(a);
	cout << "Uncomment \"" << a << "\" -> \"" << d << "\"" << endl;


	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test vector<string> uncomment(const vector<string> & _input, const string & _commentString=\"#\")" << endl;
	
	e.clear();
	e.push_back("This is a value  # this is a comment");
	e.push_back("This sentence has a \"# within quotes\" and it is not a comment");
	e.push_back("This sentence has it \"outside\" the quotes and it is a comment # here is the comment");
	e.push_back("This sentence has a \'# within single quotes\' and it is not a comment");
	e.push_back("This sentence has it \'outside\' single quotes and it is a comment # here is the comment");
	e.push_back("# everything is a comment");

	f = MslTools::uncomment(e);
	cout << "Uncomment:" << endl;
	for (vector<string>::iterator k=e.begin(); k!=e.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}
	cout << "Result:" << endl;
	for (vector<string>::iterator k=f.begin(); k!=f.end(); k++) {
		cout << ">" << *k << "<" <<  endl;
	}

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;
	cout << "Test double smartRound(double coord, double gridSize)\n";

	cout << "smartRound(4.31, 0.5) = "<< MslTools::smartRound(4.31,0.5)<<"\n";
	cout << "smartRound(4.31, 1.0) = "<< MslTools::smartRound(4.31,1.0)<<"\n";
	cout << "smartRound(4.31, 10.0) = "<<MslTools::smartRound(4.31,10.0)<<"\n";
	cout << "smartRound(5.31, 10.0) = "<<MslTools::smartRound(5.31,10.0)<<"\n";
	cout << "smartRound(4.61, 0.5) = "<< MslTools::smartRound(4.61,0.5)<<"\n";
	cout << "smartRound(4.61, 0.2) = "<< MslTools::smartRound(4.61,0.2)<<"\n";
	cout << "smartRound(4.61, 0.1) = "<< MslTools::smartRound(4.61,0.1)<<"\n";
	cout << "smartRound(4.65, 0.1) = "<< MslTools::smartRound(4.65,0.1)<<"\n";
	cout << "smartRound(4.66, 0.1) = "<< MslTools::smartRound(4.66,0.1)<<"\n";
	cout << "smartRound(90, 180) = "<< MslTools::smartRound(90,180)<<"\n";
	cout << "smartRound(179, 180) = "<< MslTools::smartRound(179,180)<<"\n";
	cout << "smartRound(270, 180) = "<< MslTools::smartRound(270,180)<<"\n";





	return 0;
}
