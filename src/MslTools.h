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

#ifndef MSLTOOLS_H_
#define MSLTOOLS_H_

// STL includes
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <limits>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <time.h>

// MSL includes
#include "MslExceptions.h"
#include "Real.h"
#include "release.h"

#ifdef __GSL__
#include "RandomNumberGenerator.h"
#endif

// BOOST Includes
#ifdef __BOOST__
#include <boost/regex.hpp>
#endif

using namespace std;


namespace MslTools {




	/*
             ******************************************
	     *          CONSTANTS
	     ******************************************
	*/

	const float  floatMax  = std::numeric_limits<float>::max();
	const double doubleMax = std::numeric_limits<double>::max();
	const int    intMax    = std::numeric_limits<int>::max();

	


	/*
             ******************************************
	     *          STRING FUNCTIONS
	     ******************************************
	*/

	string trim(const string & _str, const string & _trimString=" \t\n\r");
	vector<string> trim(const vector<string> & _str, const string & _trimString=" \t\n\r");
	
	// CONVERSIONS
	double toDouble(const string & _string, const string & _msg=(string)"");
        Real toReal(const string & _string, const string & _msg=(string)"");
	int toInt(const string & _string, const string & _msg=(string)"");
	string intToString(const int & _i, const string & _msg=(string)"");
	string doubleToString(const double & _d, const string & _msg=(string)"");

	// CHECK CONTENT OF STRING
	bool isDigitChars(string _input); // 0-9
	bool isAlphaNumericChars(string _input); // 0-9 A-Z a-z
	bool isAlphaChars(string _input); // A-Z a-z
	bool isWhiteSpaces(string _input); // " " \t \n \r

	 // split a string made by digits followed by letters into an int and a string (i.e. "734B" -> 734 and "B")
	void splitIntAndString(const string & _input, int & _intResult, string & _stringResult);

	string toUpper(const string & _input);
	string toLower(const string & _input);

	// JOIN AND SPLIT
	vector<string> tokenize(const string & _input, const string & _delimiter=" ", bool _allowEmtpy=false);
	vector<vector<string> > dualLevelTokenizer(const string & _input, const string & _delimiter=" ", const string & _leftBraket="[", const string _rightBraket="]");
	vector<string> extractBraketed(const string & _input, const string & _leftBraket="[", const string _rightBraket="]");
	//vector<string> tokenizeWord(const string & _input, const string & _delimiter=" ", const bool & _trim=false);
	vector<string> tokenizeAndTrim(const string & _input, const string & _delimiter=" ", bool _allowEmtpy=false, const string & _trimString=" \t\n\r");

	string joinLines(const vector<string> & _input, const string & _spacer=" ");

	vector<string> joinConnectedLines(const vector<string> & _input, string _marker="\\", string _spacer=" ");
	//vector<string> joinBackslashedLines(const vector<string> & _input, const string & _spacer=" ");

	vector<string> removeEmptyLines(const vector<string> & _input);
	string uncomment(const string & _input, const string & _commentString="#");
	vector<string> uncomment(const vector<string> & _input, const string & _commentString="#");
	

	// Extract the name of a file from a full path , subtract extension as well.
	bool readTextFile(vector<string> & _container, const string & _filename);

	string getMSLversion();


	// RegEx Functions  (only works with compile __BOOST__ flag on, otherwise returns false immediately)
	//   Remember escaped characters need to be double-escaped in expression variable so '\s' is '\\s' 
	bool regex(string lineToMatch, string expression, vector<string> &matches);
		
	/*
             ******************************************
	     *          MATH FUNCTIONS
	     ******************************************
	*/

	double correlate(vector<double> _data1, vector<double> _data2);
	double smartRound(double coord, double gridSize);
	double mod(double x, double y);
	double setPrecision(double _d, unsigned int _significantDigits);

	Real round(Real value);

	
	/*
             ******************************************
	     *          PATH FUNCTIONS
	     ******************************************
	*/
	string pathRoot(string _path); // /the/path/theFile.ext -> /the/path/theFile
	string pathHead(string _path); // /the/path/theFile.ext -> /the/path
	string pathExtension(string _path); // /the/path/theFile.ext -> ext
	string pathTail(string _path); // /the/path/theFile.ext -> theFile.ext
	string getFileName(string fullpath);

	string createDir(string _name);
	string outputFileNameParser(string _name);
	bool   mkNestedDir(string _dir, mode_t _mode);
	string getRandomAlphaNumString(unsigned int _size, bool _alphaOnly=false);
	unsigned int getRandomInt(unsigned int _max);

	/*
              ******************************************
              * Tools to convert different AA codes.
	      ******************************************
	*/
	static std::map<string, string> threeToOneLetter;
	static std::map<string, string> oneToThreeLetter;
	void loadAAConversionTables();
	string getOneLetterCode(string threeLetterCode);
	string getThreeLetterCode(string oneLetterCode);



	/*
              ******************************************
              * STL Helper functions (sort predicates,etc..)
	      ******************************************
        */

	bool sortPairIntDoubleAscending(const pair<int,double> &left, const pair<int,double> &right);
	bool sortPairIntDoubleDecending(const pair<int,double> &left, const pair<int,double> &right);

	/*
              ******************************************
              * Help sorting residues by resnum and icode
	      ******************************************
        */
	bool sortByResnumIcodeAscending(int _resnum1, string _icode1, int _resnum2, string _icode2);

        template <class A, class B> inline bool mapHasKey(std::map<A, B> &_mapObj, A _key) {
            return( _mapObj.find(_key) != _mapObj.end());
        };


	/*
              ******************************************
              * Sorting a vector of doubles
	      ******************************************
	*/
	void quickSort(vector<double> & _vec);
	void quickSort(vector<double> & _vec, unsigned int _left, unsigned int _right);
	double partition(vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex);
	/*
              ******************************************
              * Sorting a vector of doubles and keeping track of
	      * the original index in another vector (_index
	      * should initially be _index[i] := i)
	      ******************************************
	*/
	void quickSortWithIndex(vector<double> & _vec, vector<unsigned int> &_index);
	double partitionWithIndex(vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex, vector<unsigned int> &_index);
	void quickSortWithIndex(vector<double> & _vec, unsigned int _left, unsigned int _right, vector<unsigned int> &_index);


	/*
             ******************************************
	     *          COLOR FUNCTIONS
	     ******************************************
	*/
	vector<double> getRGB(vector<double> &_startRGB, vector <double> &_endRGB, double _minValue, double _maxValue, double _value);
	void rgb2hsv(vector<double> &_rgb, vector<double> &_hsv);
	void hsv2rgb(vector<double> &_hsv, vector<double> &_rgb);
	

};

// inlined functions

#endif
