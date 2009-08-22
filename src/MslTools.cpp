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

#include "MslTools.h"
#include "MslExceptions.h"
#include "RandomNumberGenerator.h"

#include <sys/stat.h>
#include <errno.h>
#include "math.h"
#include <iostream>

string MslTools::trim(const string & _str, const string &_trimString){
	string out = _str;

	// trim leading whitespace
	string::size_type  notwhite = out.find_first_not_of(_trimString);
	out.erase(0,notwhite);

	// trim trailing whitespace
	notwhite = out.find_last_not_of(_trimString); 
	out.erase(notwhite+1); 

	return out;
}

vector<string> MslTools::trim(const vector<string> & _str,const string &_trimString){
	vector<string> out;
	for (vector<string>::const_iterator k=_str.begin(); k!=_str.end(); k++) {
		out.push_back(trim(*k,_trimString));
	}
	return out;
}

double MslTools::toDouble(const string & _string, const string & _msg) {
	istringstream ss(_string);
	double d = 0.0;
	ss >> d;

	if (!ss) throw ConvertDoubleException(" string "+ _string +" is not a double. " + _msg);

	return d;
}

Real MslTools::toReal(const string & _string, const string & _msg) {
	istringstream ss(_string);
	Real d = 0.0;
	ss >> d;

	if (!ss) throw ConvertDoubleException(" string "+ _string +" is not a Real. " + _msg);

	return d;
}

int MslTools::toInt(const string & _string, const string & _msg) {

	istringstream ss(_string);
	int i = 0;
	ss >> i;

	if (!ss) throw ConvertIntException(" string "+ _string +" is not an int. "+ _msg);

	return i;
}

string MslTools::intToString(const int & _i, const string & _msg){
	stringstream ss;
	ss << _i;

	return ss.str();
	
}

string MslTools::doubleToString(const double & _d, const string & _msg){
	stringstream ss;
	ss << _d;

	return ss.str();
	
}

string MslTools::toUpper(const string & _input) {
	string out = _input;
	for(int i=0; i<out.length(); i++) {
		out[i] = toupper(out[i]);
	}
	return out;
}

string MslTools::toLower(const string & _input) {
	string out = _input;
	for(int i=0; i<out.length(); i++) {
		out[i] = tolower(out[i]);
	}
	return out;
}

vector<string> MslTools::tokenize(const string & _input, const string & _delimiter, bool _allowEmtpy){
	vector<string> results;
	
	if (_allowEmtpy) {
		size_t prePos = 0;
		size_t pos  = _input.find(_delimiter);
        unsigned int delimiterSize = _delimiter.size();
		string left = _input, right;

		while (pos != std::string::npos) {
            results.push_back(left.substr(prePos, pos));
            if( pos + delimiterSize <= left.size() )
                left = left.substr(pos + delimiterSize, left.size() );
            else
                left = "";
			pos  = left.find(_delimiter);
		}

		results.push_back(left);
	} else {
		int start  = _input.find_first_not_of(_delimiter);
		int end    = 0;
		string cur = _input;


		while (start != std::string::npos){
			end    = _input.find_first_of(_delimiter, start);
			results.push_back(_input.substr(start, end-start));
			start  = _input.find_first_not_of(_delimiter, end);
		}
	}
	return results;

}

vector<string> MslTools::tokenizeAndTrim(const string & _input, const string & _delimiter, bool _allowEmtpy, const string & _trimString) {
	vector<string> results = tokenize(_input, _delimiter, _allowEmtpy);
	return trim(results, _trimString);
}

/*
vector<string> MslTools::tokenizeWord(const string & _input, const string & _delimiter, const bool & _trim){
	vector<string> results;

	int start  = _input.find_first_not_of(_delimiter);
	int end    = 0;
	int loop   = 0;
	while ( start != std::string::npos){
		if (loop == 0) start = 0;
		end    = _input.find(_delimiter, start);

		if (_trim){
			results.push_back(trim(_input.substr(start, end-start)));
		} else {
			results.push_back(_input.substr(start, end-start));			
		}
		start  = _input.find_first_not_of(_delimiter,end);
		loop++;
	}

	return results;
}
*/


vector<vector<string> > MslTools::dualLevelTokenizer(const string & _input, const string & _delimiter, const string & _leftBraket, const string _rightBraket) {
	vector<vector<string> > out;

	unsigned int cycle = 0;
	vector<string> split = extractBraketed(_input, _leftBraket, _rightBraket);
	while (true) {
		vector<string> tokens0 = tokenize(split[0], _delimiter);
		for (vector<string>::iterator k=tokens0.begin(); k!=tokens0.end(); k++) {
			out.push_back(vector<string>(1, *k));
		}
		vector<string> tokens1 = tokenize(split[1], _delimiter);
		if (tokens1.size() > 0) {
			out.push_back(tokens1);
		}
		split = extractBraketed(split[2], _leftBraket, _rightBraket);
		cycle++;
		if (trim(split[2]).length() == 0) {
			vector<string> tokens0 = tokenize(split[0], _delimiter);
			for (vector<string>::iterator k=tokens0.begin(); k!=tokens0.end(); k++) {
				out.push_back(vector<string>(1, *k));
			}
			vector<string> tokens1 = tokenize(split[1], _delimiter);
			if (tokens1.size() > 0) {
				out.push_back(tokens1);
			}
			break;
		}
	}
	return out;
}

vector<string> MslTools::extractBraketed(const string & _input, const string & _leftBraket, const string _rightBraket) {
	vector<string> out;
	size_t posL  = _input.find(_leftBraket);
	size_t posR  = _input.find(_rightBraket, posL);
	if (posL != std::string::npos && posR != std::string::npos) {
		out.push_back(_input.substr(0, posL));
		out.push_back(_input.substr(posL+1, posR-posL-1));
		out.push_back(_input.substr(posR+1, _input.length()-posR-1));
	} else {
		out.push_back(_input);
		out.push_back("");
		out.push_back("");
	}
	return out;
}


void MslTools::splitIntAndString(const string & _input, int & _intResult, string & _stringResult) {
	_intResult = 0;
	_stringResult = "";

	string trimmed = trim(_input);

	int scaleFactor = 1;
	int startIndex  = 0;
	if (trimmed[0] == '-'){
		scaleFactor = -1;
		startIndex = 1;
	}
	for (int i=startIndex; i<trimmed.size(); i++) {

		int asciiCode = trimmed[i];

		if (asciiCode >= 48 && asciiCode <= 57) {
			// it is a digit
			int digit = trimmed[i] - '0';
			_intResult *= 10;
			_intResult += digit;
		} else {
			_stringResult = trim(trimmed.substr(i, trimmed.size() - i));
			_intResult *= scaleFactor;
			return;
		}
	}

	_intResult *= scaleFactor;
}

string MslTools::joinLines(const vector<string> & _input, const string & _spacer) {
	string out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		if (k != _input.begin()) {
			out += _spacer;
		}
		out += *k;
	}
	return out;
}

vector<string> MslTools::joinConnectedLines(const vector<string> & _input, string _marker, string _spacer) {
	/******************************************************
	 *
	 *  Connect lines if they have a special terminating character
	 *  (the defauls is a backlash)
	 *
	 *  Example 
	 *
	 *    This is a first line
	 *    This is a second line
	 *    This is a third line, backslashed,\
	 *    and this is the continuation of the third line
	 *    This is a fourth line
	 *
	 *  gives:
	 *    This is a first line
	 *    This is a second line
	 *    This is a third line, backslashed, and this is the continuation of the third line
	 *    This is a fourth line
	 *******************************************************/

	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		out.push_back("");
		bool first = true;
		while (k->size() > 0 && k->substr(k->size()-1, 1) == _marker) {
			if (first) {
				first = false;
			} else {
				out.back() += _spacer;
			}
			out.back() += k->substr(0, k->size()-1);
			k++;
			if (k==_input.end()) {
				return out;
			}
		}
		if (first) {
			first = false;
		} else {
			out.back() += _spacer;
		}
		out.back() += *k;
	}
	return out;
}


string MslTools::uncomment(const string & _input, const string & _commentString) {

	size_t foundComment = _input.find(_commentString);
	if (foundComment == string::npos) {
		// no comment string, we are done
		return _input;
	}

	bool openSingleQuote = false;
	bool openDoubleQuote = false;
	for (unsigned int i=0; i< _input.size() - _commentString.size() + 1; i++) {
		// truncate a line at the # character
		if (!openSingleQuote && !openDoubleQuote && _input.substr(i, _commentString.size()) == _commentString) {
			return _input.substr(0,i);
		}
		// but do not consider if the # is within quotes
		if (!openSingleQuote && _input[i] == '"') {
			openDoubleQuote = !openDoubleQuote;
		}
		if (!openDoubleQuote && _input[i] == '\'') {
			openSingleQuote = !openSingleQuote;
		}
	}
	return _input;
}

vector<string> MslTools::uncomment(const vector<string> & _input, const string & _commentString) {
	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		out.push_back(uncomment(*k, _commentString));
	}
	return out;
}

vector<string> MslTools::removeEmptyLines(const vector<string> & _input) {
	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		// remove empty lines
		if (trim(*k).size() != 0) {
			out.push_back(*k);
		}
	}
	return out;
}



string MslTools::getFileName(string fullpath){

  string name = fullpath;

  // Assume linux paths '/' , sorry windows users!
  vector<string> paths = tokenize(name, "/");


  if (paths.size() > 0){
    name = paths[paths.size()-1];
  }

  int cut;
  if (name.find_last_of(".") > 1000 || name.find_last_of(".") <= 0){
	  cut = name.length();
  } else {
	  cut = name.length() - name.find_last_of(".");
  }

  name.erase(name.length() - cut);

  return name;
}




double MslTools::correlate(vector<double> _data1, vector<double> _data2){

	// GSL Way... (GSL 1.10+ required)
	/*
	double d1[_data1.size()];
	double d2[_data1.size()];
	for (uint i = 0; i < _data1.size();i++){
		d1[i] = _data1[i];
		d2[i] = _data2[i];
	}

	// This gets the Pearson correlation coefiicient (the 1's are strides: skipping number in array)
	//  double r =  gsl_stats_correlation(d1,1,d2,1,_data1.size());  

	*/         


	// Pearson product-moment correlation coefficient 
	//                               n sum(x y) - sum(x) sum(y)
	//        r(x, y) := ---------------------------------------------------
	//                                        2                           2
	//                   sqrt(n sum(x x) - sum (x)) sqrt(n sum(y y) - sumy )


	double sumX;
	double sumY;
	double sumXY;
	double sumXX;
	double sumYY;
	sumX = sumY = sumXY = sumXX = sumYY = 0.0;
	for (uint i = 0 ; i < _data1.size();i++){

		sumX += _data1[i];
		sumY += _data2[i];
		sumXY += (_data1[i] *_data2[i]);
		sumXX += (_data1[i] *_data1[i]);
		sumYY += (_data2[i] *_data2[i]);
		
	}

	int n = _data1.size();
	double r =  (n * sumXY - sumX*sumY) / (sqrt(n * sumXX - (sumX*sumX)) * sqrt(n * sumYY - (sumY*sumY)));
	return r;
}



bool MslTools::readTextFile(vector<string> & _container, const string & _filename) {

	_container.clear();

	ifstream fs;
	fs.open(_filename.c_str());
	if (fs.fail()) {
		return false;
	}

	while (true) {
		string line;
		getline(fs,line);

		if (fs.fail()) {
			// no more lines to read from file, quit the while loop
			break;
		}

		_container.push_back(line);

	}
	fs.close();

	return true;

}



string MslTools::pathRoot(string _path) {
	// /the/path/theFile.ext -> /the/path/theFile
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path;
		} else if (_path.substr(i,1) == ".") {
			return _path.substr(0,i);
		}
	}
	return _path;
}

string MslTools::pathHead(string _path) {
	// /the/path/theFile.ext -> /the/path
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path.substr(0,i);
		}
	}
	return _path;
}

string MslTools::pathExtension(string _path) {
	// /the/path/theFile.ext -> ext
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return "";
		} else if (_path.substr(i,1) == ".") {
			return _path.substr(i+1,_path.size()-i-1);
		}
	}
	return "";
}

string MslTools::pathTail(string _path) {
	// /the/path/theFile.ext -> theFile.ext
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path.substr(i+1,_path.size()-i-1);
		}
	}
	return _path;
}




double MslTools::smartRound(double coord, double gridSize){

  if (gridSize > 1) {
    int cint = int(coord*100);
    int gint = int(gridSize *100);
    double remainder = gint - (cint % gint);
    if (remainder / (double)gint < 0.5)
      return int(coord + remainder/100);

    return int((coord*100 - (gridSize*100 - remainder))/100);

  } else {

    double frac_coord       = coord - int(coord);
    int shift_frac_coord = int(100 * frac_coord);
    double correction = int(shift_frac_coord % int(gridSize *100));
    double frac_result = (shift_frac_coord - correction)/100;
    if (correction/(int (gridSize*100)) >= 0.5){
      frac_result += gridSize;
      if (frac_result > 1) frac_result = int(frac_result);	
    }

    return int(coord) + (frac_result);
  }
}


double MslTools::mod(double x, double y){
	return x-y *floor(x/y);
}


Real MslTools::round(Real value) {
    Real half = (Real) 0.5f;
    int val = (int)( (value > (Real)0.0f)?(value + half):(value - half) );
    return( (Real)val );
}


string MslTools::createDir(string _name){

	// First make sure name is unique
	string prevOutputdir = _name;
	string outputdir = outputFileNameParser(_name);
	while (true) {
		if (prevOutputdir == outputdir) {
			// the outputdir doesn't have randomly generated parts
			break;
		}
		ifstream test_fs;
		test_fs.open(outputdir.c_str());
		if (test_fs.is_open()) {
			// directory exists, try again
			prevOutputdir = outputdir;
			outputdir = outputFileNameParser(outputdir);
			test_fs.close();
		} else {
			test_fs.close();
			break;
		}
	}

	// Now make the directory
	mkNestedDir(outputdir, 0755); 

	return outputdir;
}

string MslTools::outputFileNameParser(string _name) {
	vector<unsigned int> startSub;
	vector<unsigned int> endSub;
	vector<unsigned int> size;
	vector<string> type;
	for (unsigned int i=0; i<_name.size(); i++) {
		if (_name.substr(i, 1) == "%") {
			unsigned int start = i;
			i++;
			bool invalid = false;
			string digits;
			while (i<_name.size() && _name.substr(i,1) != "s" && _name.substr(i,1) != "r" && _name.substr(i,1) != "t" && _name.substr(i,1) != "d") {
				if (_name[i] > 47 && _name[i] < 58) {
					digits += _name.substr(i,1);
					i++;
				} else {
					invalid = true;
					break;
				}
			}
			if (!invalid) {
				if (digits.size() == 0) {
					digits = "1";
				}
				unsigned int number = toInt(digits);
				if ( _name.substr(i,1) != "s" || _name.substr(i,1) != "r" || _name.substr(i,1) != "t" || _name.substr(i,1) != "d") {
					startSub.push_back(start);
					endSub.push_back(i);
					size.push_back(number);
					type.push_back(_name.substr(i,1));
				}
			}
		}
	}

	time_t currentTime = time(NULL);; 
	tm *curr = localtime(&currentTime);
	//time(&currentTime);
	stringstream ss;
	ss << (int)currentTime;
	string timeString = ss.str();
	char c [1000];
	sprintf(c, "%04d%02d%02d", curr->tm_year+1900, curr->tm_mon+1, curr->tm_mday);
	string dateStamp = c;
	sprintf(c, "%02d%02d%02d", curr->tm_hour, curr->tm_min, curr->tm_sec);
	string timeStamp = c;

	for (int i=startSub.size()-1; i>=0; i--) {
		if (type[i] == "s") {
			while (timeString.size() < size[i]) {
				timeString = "0" + timeString;
			}
			if (timeString.size() > size[i]) {
				timeString = timeString.substr(timeString.size() - size[i], size[i]);
			}
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, timeString);
		} else if (type[i] == "d") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, dateStamp);
		} else if (type[i] == "t") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, timeStamp);
		} else if (type[i] == "r") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, getRandomAlphaNumString(size[i],true));
		} 
	}
	return _name;

}

bool MslTools::mkNestedDir(string _dir, mode_t _mode) {

	_dir = trim(_dir);
	vector<string> directories;
	vector<string>::iterator k;

	bool absolutePath = false;
	string prevPath = "";
	if (_dir.substr(0,1) == "/") {
		absolutePath = true;
		prevPath = "/";
		_dir = _dir.substr(1,_dir.size()-1);
	}
	
	bool open = false;
	for (int i=0; i< _dir.size(); i++) {
		if (_dir.substr(i,1) == "/") {
			open = false;
			continue;
		} else {
			if (open) {
				*k += _dir[i];
			} else {
				open = true;
				directories.push_back(_dir.substr(i,1));
				k = directories.end() - 1;
			}
		}
	}

	for (unsigned int i=0; i<directories.size(); i++) {
		// don't make . and ..
		if (directories[i] != "." || directories[i] != "..") {
			string dirToMake = prevPath + directories[i];
			
			if (mkdir(dirToMake.c_str(), _mode) != 0) {
				// failed because it did already exist?
				if (errno != EEXIST) {
					// no, we failed
					return false;
				}
			}
			prevPath += directories[i] + (string)"/";
		}
	}
	return true;
}

string MslTools::getRandomAlphaNumString(unsigned int _size, bool _alphaOnly) {
	vector<string> characters;
	characters.push_back("a");
	characters.push_back("b");
	characters.push_back("c");
	characters.push_back("d");
	characters.push_back("e");
	characters.push_back("f");
	characters.push_back("g");
	characters.push_back("h");
	characters.push_back("i");
	characters.push_back("j");
	characters.push_back("k");
	characters.push_back("l");
	characters.push_back("m");
	characters.push_back("n");
	characters.push_back("o");
	characters.push_back("p");
	characters.push_back("q");
	characters.push_back("r");
	characters.push_back("s");
	characters.push_back("t");
	characters.push_back("u");
	characters.push_back("v");
	characters.push_back("w");
	characters.push_back("x");
	characters.push_back("y");
	characters.push_back("z");
	characters.push_back("A");
	characters.push_back("B");
	characters.push_back("C");
	characters.push_back("D");
	characters.push_back("E");
	characters.push_back("F");
	characters.push_back("G");
	characters.push_back("H");
	characters.push_back("I");
	characters.push_back("J");
	characters.push_back("K");
	characters.push_back("L");
	characters.push_back("M");
	characters.push_back("N");
	characters.push_back("O");
	characters.push_back("P");
	characters.push_back("Q");
	characters.push_back("R");
	characters.push_back("S");
	characters.push_back("T");
	characters.push_back("U");
	characters.push_back("V");
	characters.push_back("W");
	characters.push_back("X");
	characters.push_back("Y");
	characters.push_back("Z");
	characters.push_back("0");
	characters.push_back("1");
	characters.push_back("2");
	characters.push_back("3");
	characters.push_back("4");
	characters.push_back("5");
	characters.push_back("6");
	characters.push_back("7");
	characters.push_back("8");
	characters.push_back("9");
	string out;

	int size = characters.size();
	if (_alphaOnly){
		size = 52; //  53 ?
	}

	RandomNumberGenerator rng;
	rng.setRNGType("knuthran2");
	rng.setRNGTimeBasedSeed();

	for (unsigned int i=0; i<_size; i++) {
		int randomN = (int)(rng.getRandomInt() % size);
		out += characters[randomN];
	}
	
	return out;

}

void MslTools::loadAAConversionTables() {
	threeToOneLetter["ALA"] = "A";
	oneToThreeLetter["A"]   = "ALA";

	threeToOneLetter["CYS"] = "C";
	oneToThreeLetter["C"]   = "CYS";

	threeToOneLetter["ASP"] = "D";
	oneToThreeLetter["D"]   = "ASP";

	threeToOneLetter["GLU"] = "E";
	oneToThreeLetter["E"]   = "GLU";

	threeToOneLetter["PHE"] = "F";
	oneToThreeLetter["F"]   = "PHE";

	threeToOneLetter["GLY"] = "G";
	oneToThreeLetter["G"]   = "GLY";

	threeToOneLetter["HIS"] = "H";
	oneToThreeLetter["H"]   = "HIS";

	threeToOneLetter["ILE"] = "I";
	oneToThreeLetter["I"]   = "ILE";

	threeToOneLetter["LYS"] = "K";
	oneToThreeLetter["K"]   = "LYS";

	threeToOneLetter["LEU"] = "L";
	oneToThreeLetter["L"]   = "LEU";

	threeToOneLetter["MET"] = "M";
	oneToThreeLetter["M"]   = "MET";

	threeToOneLetter["ASN"] = "N";
	oneToThreeLetter["N"]   = "ASN";

	threeToOneLetter["PRO"] = "P";
	oneToThreeLetter["P"]   = "PRO";

	threeToOneLetter["GLN"] = "Q";
	oneToThreeLetter["Q"]   = "GLN";

	threeToOneLetter["ARG"] = "R";
	oneToThreeLetter["R"]   = "ARG";

	threeToOneLetter["SER"] = "S";
	oneToThreeLetter["S"]   = "SER";

	threeToOneLetter["THR"] = "T";
	oneToThreeLetter["T"]   = "THR";

	threeToOneLetter["VAL"] = "V";
	oneToThreeLetter["V"]   = "VAL";

	threeToOneLetter["TRP"] = "W";
	oneToThreeLetter["W"]   = "TRP";

	threeToOneLetter["TYR"] = "Y";
	oneToThreeLetter["Y"]   = "TYR";
}
string MslTools::getOneLetterCode(string threeLetterCode) { 
    string oneLetterCode = "X";
    if (threeToOneLetter.size()<= 0) 
        loadAAConversionTables();
    threeLetterCode = toUpper(threeLetterCode);
    map<string, string>::iterator entry = threeToOneLetter.find(threeLetterCode);
    if (entry != threeToOneLetter.end())
        oneLetterCode = entry->second;
    return oneLetterCode;
}

string MslTools::getThreeLetterCode(string oneLetterCode) { 
    string threeLetterCode = "XXX";
    if (oneToThreeLetter.size()<= 0) 
        loadAAConversionTables();
    oneLetterCode = toUpper(oneLetterCode);
    map<string, string>::iterator entry = oneToThreeLetter.find(oneLetterCode);
    if (entry != oneToThreeLetter.end())
        threeLetterCode = entry->second;
    return threeLetterCode;
}


bool MslTools::sortPairIntDoubleAscending(const pair<int,double> &left, const pair<int,double> &right){
	return left.second < right.second;
}
bool MslTools::sortPairIntDoubleDecending(const pair<int,double> &left, const pair<int,double> &right){
	return left.second > right.second;
}



double MslTools::setPrecision(double _d, unsigned int _significantDigits) {
	return (double)int(_d * pow(10.0, _significantDigits) + 0.5) / pow(10.0, _significantDigits);
}




bool MslTools::sortByResnumIcodeAscending(int _resnum1, string _icode1, int _resnum2, string _icode2) {
	if (_resnum1 < _resnum2 || (_resnum1 == _resnum2 && _icode1.c_str()[0] < _icode1.c_str()[1])) {
		return true;
	}
	return false;
}

