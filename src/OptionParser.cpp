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

#include "OptionParser.h"
#include "release.h"
#include "MslOut.h"

using namespace MSL;
using namespace std;

static MslOut MSLOUT("OptionParser");

OptionParser::OptionParser() {
	errorFlag = false;
	commandName = "";
	confFile = "";
	translateEnv_flag = true;
}

OptionParser::OptionParser(const OptionParser & _OP) {
	copy(_OP);
}

OptionParser::~OptionParser() {
}

void OptionParser::operator=(const OptionParser & _OP) {
	copy(_OP);
}

void OptionParser::copy(const OptionParser & _OP) {
	opts = _OP.opts;
	vals = _OP.vals;
	source = _OP.source;
	freeArgs = _OP.freeArgs;

	shortOptEquivalent = _OP.shortOptEquivalent;

	commandName = _OP.commandName;

	confFile = _OP.confFile;

	required = _OP.required;
	allowed = _OP.allowed;
	mutualExclusive = _OP.mutualExclusive;
	dependentOn = _OP.dependentOn;
	interDependent = _OP.interDependent;
	linked = _OP.linked;
	oneRequired = _OP.oneRequired;

	missing = _OP.missing;
	disallowed = _OP.disallowed;
	ambiguous = _OP.ambiguous;
	disallowedTogether = _OP.disallowedTogether;
	missingDependency = _OP.missingDependency;
	missingInterdipendenty = _OP.missingInterdipendenty;
	missingLinked = _OP.missingLinked;
	missingOneRequired = _OP.missingOneRequired;

	defaultArguments = _OP.defaultArguments;
	
	errorFlag = _OP.errorFlag;
	translateEnv_flag = _OP.translateEnv_flag;	
}


void OptionParser::readArgv(int theArgc, char *theArgv[]) {

	/***************************************************
	 *  Parse the command line arguments
	 ***************************************************/

	errorFlag = false;

	/***************************************************
	 *  save the command (argv[0])
	 ***************************************************/
	vector<string> commandLineArg;
	for (int i=0; i<theArgc; i++) {
		if (i==0) {
			commandName = theArgv[i];
		} else {
			commandLineArg.push_back(theArgv[i]);
		}
	}

	for (int i=0; i<commandLineArg.size(); i++) {
		/***************************************************
		 *  a long argument has 2 dashes and multiple characters (i.e --help) 
		 ***************************************************/
		if (commandLineArg[i].size() > 3 && commandLineArg[i].substr(0,2) == "--") {
			string s = commandLineArg[i].substr(2,commandLineArg[i].size()-2);

			string equalVal;
			for (unsigned int j=1; j<s.size(); j++) {
				if (s.substr(j, 1) == "=") {
					equalVal = s.substr(j+1,s.size()-j-1);
					s = s.substr(0, j);
				}
			}

			bool found = false;
			int foundAt = -1;
			
			/***************************************************
			 *  If the same argument was aready given ignore it
			 ***************************************************/
			for (int j=0; j<opts.size(); j++) {
				if (s == opts[j]) {
					found = true;
					foundAt = j;
					break;
				}
			}

			if (!found) {
				opts.push_back(s);
				vals.push_back(vector<string>());
				source.push_back(vector<string>());
				foundAt = opts.size() - 1;
			}
			vector<string> valvector;
			if (equalVal != "") {
				valvector.push_back(equalVal);
			}
			while (i<commandLineArg.size()-1 && commandLineArg[i+1].substr(0,1) != "-") {
				i++;
				valvector.push_back(commandLineArg[i]);
			}
			string v = MslTools::joinLines(valvector);
			vals[foundAt].push_back(v);
			source[foundAt].push_back("a");
		} else {

			/***************************************************
			 *  a short argument has 1 dash and one character (i.e --h).
			 *  they can be aggregated (i.e -h -v -> -hv)
			 ***************************************************/
			if (commandLineArg[i].size() > 1 && commandLineArg[i].substr(0,1) == "-" && commandLineArg[i].substr(1,1) != "-") {
				string aggregate = commandLineArg[i].substr(1,commandLineArg[i].size()-1);
				for (int k=0; k<aggregate.size(); k++) {
					string s = aggregate.substr(k,1);
					bool found = false;
					int foundAt = -1;
					
					// if the short option was linked to a long option (i.e. -h -> --help) do
					// substitution
					string letter = "o";
					for (unsigned int j=0; j<shortOptEquivalent.size(); j++) {
						if (s == shortOptEquivalent[j][0]) {
							s = shortOptEquivalent[j][1];
							letter = "a";
							break;
						}
					}
					/***************************************************
					 *  If the same argument was aready given ignore it
					 ***************************************************/
					for (int j=0; j<opts.size(); j++) {
						if (s == opts[j]) {
							found = true;
							foundAt = j;
							break;
						}
					}

					if (!found) {
						opts.push_back(s);
						vals.push_back(vector<string>());
						source.push_back(vector<string>());
						foundAt = opts.size() - 1;
					}
					string v = "";
					/***************************************************
					 *  Only the last letter of an aggregate set can have
					 *  arguments
					 ***************************************************/
					if (k == aggregate.size() - 1) {
						vector<string> valvector;
						while (i<commandLineArg.size()-1 && commandLineArg[i+1].substr(0,1) != "-") {
							i++;
							valvector.push_back(commandLineArg[i]);
						}
						/*
						for (int j=0; j<valvector.size(); j++) {
							if (j==0) {
								v = valvector[j];
							} else {
								v += " " + valvector[j];
							}
						}
						*/
						v = MslTools::joinLines(valvector);
					}
					vals[foundAt].push_back(v);
					source[foundAt].push_back(letter);
				}

			} else {
				freeArgs.push_back(commandLineArg[i]);
			}
		}

	}

	addFreeArgumentsToDefault();



	if (getBool("mslVersion")){
		cout << "MSLVERSION: "<<MSLVERSION << endl<<"MSLDATE: "<<MSLDATE<<endl;
		exit(0);
	}

	if (getBool("printOptions")){

		cout << "program options: "<<endl;
		for (uint i = 0; i < required.size();i++){
			cout <<"R  --"<<required[i]<<"  "<<endl;
		}
		cout <<endl;
		for (uint i = 0; i < allowed.size();i++){
			cout <<"O  --"<<allowed[i]<<"  "<<endl;
		}
		cout << endl;
		exit(0);
	}


	// Default is to turn all output off
	MSLOUT.turnAllOff();	    

	// Turn all output for all objects on
	if (getBool("speakAll")){
		MSLOUT.turnAllOn();
	}

	// Turn on output for specific objects
	std::vector<std::string> objectsToSpeak = getMultiString("speak");
	for (uint i = 0; i < objectsToSpeak.size();i++){
		MSLOUT.turnOn(objectsToSpeak[i]);
	}

	// Turn off output for specific objects
	std::vector<std::string> objectsToMute = getMultiString("mute");
	for (uint i = 0; i < objectsToMute.size();i++){
		MSLOUT.turnOff(objectsToMute[i]);
	}

}


void OptionParser::setDefaultArgument(string _argument) {
	defaultArguments.clear();
	defaultArguments.push_back(_argument);
	addFreeArgumentsToDefault();
}

void OptionParser::setDefaultArguments(vector<string> _argument) {
	defaultArguments = _argument;
	addFreeArgumentsToDefault();
}

void OptionParser::addFreeArgumentsToDefault() {
	vector<string>::iterator k;
	
	for (unsigned int i=0; i<defaultArguments.size(); i++) {
		if (freeArgs.size() == 0) {
			return;
		}

		// put on tmp the arguments to be put onto the default argument
		vector<string> tmp;
		if (i == defaultArguments.size() - 1) {
			// last argument put them all
			tmp = freeArgs;
			freeArgs.clear();
		} else {
			// put just one
			tmp.push_back(freeArgs[0]);
			k = freeArgs.begin();
			freeArgs.erase(k);
		}

		// has it been already given?
		bool found = false;
		int foundAt = 0;
		for (unsigned int j=0; j<opts.size(); j++) {
			if (defaultArguments[i] == opts[j]) {
				found = true;
				foundAt = j;
				break;
			}
		}

		if (!found) {
			opts.push_back(defaultArguments[i]);
			vals.push_back(vector<string>());
			source.push_back(vector<string>());
			foundAt = opts.size() - 1;
		}

		string v = MslTools::joinLines(tmp);
		vals[foundAt].insert(vals[foundAt].begin(), v);
		source[foundAt].insert(source[foundAt].begin(), "a");
	}

}

bool OptionParser::readFile(string filename) {
	
	errorFlag = false;
	vector<string> fileLines;
	ifstream fs;
	fs.open(filename.c_str());
	if (fs.fail()) {
		//cerr << "WARNING 1254: Cannot open " << filename << endl;
		errorFlag = true;
		return false;
	}

	fileLines.clear();	
	while (true) {
		string line;
		getline(fs,line);
		if (fs.fail()) {
			// no more lines to read from file
			break;
		}
		fileLines.push_back(line);

	}
	fs.close();

	// join lines that end with a backslash
	//fileLines = MslTools::joinBackslashedLines(fileLines);
	fileLines = MslTools::joinConnectedLines(fileLines, "\\");
	fileLines = MslTools::uncomment(fileLines);
	fileLines = MslTools::trim(fileLines);
	fileLines = MslTools::removeEmptyLines(fileLines);

	for (int i=0; i<fileLines.size(); i++) {
		vector<string> s = MslTools::tokenize(fileLines[i]);
		bool found = false;
		int foundAt = -1;

		for (int j=0; j<opts.size(); j++) {
			if (s[0] == opts[j]) {
				found = true;
				foundAt = j;
				break;
			}
		}

		if (!found) {
			opts.push_back(s[0]);
			vals.push_back(vector<string>());
			source.push_back(vector<string>());
			foundAt = opts.size() - 1;
		}
		/*
		string v = "";
		for (int j=1; j<s.size(); j++) {
			if (j==1) {
				v = s[j];
			} else {
				v += " " + s[j];
			}
		}
		*/
		vector<string>::iterator k = s.begin();
		s.erase(k);
		string v = MslTools::joinLines(s);
		vals[foundAt].push_back(v);
		source[foundAt].push_back("f");
	}
	return true;

}


void OptionParser::printConfFile() const {
	cout << getConfFile();
}


bool OptionParser::writeConfFile(string filename) {
	errorFlag = false;
	vector<string> fileLines;
	ofstream fs;
	fs.open(filename.c_str());
	if (fs.fail()) {
		errorFlag = true;
		return false;
	}
	//fs << makeRunConfFile();
	fs << getConfFile();
	return true;
}

string OptionParser::getConfFile() const {
	string out;

	out += "########################################################\n";
	out += "#  Options from command line arguments:\n";
	out += "#\n";
	out += "#  NOTE: command line options take priority over the\n";
	out += "#        ones entered with the configuration file.\n";
	out += "########################################################\n";
	for (int i=0; i<opts.size(); i++) {
		//cout << i << " " << opts[i] << " " << source[i][0].size() << endl;
		// MOD HERE
		//if (source[i] == "a") {
		for (int j=0; j<source[i].size(); j++) {
			if (source[i][j] == "a" || source[i][j] == "o") {
				//out += opts[i];
				string line = opts[i];
				while (line.size() < 30) {
					line += " ";
				}
				//out += "    ";
				// MOD HERE
				//out += vals[i];
				line += vals[i][j];
				line += "\n";
				out += line;
			}
		}
	}
	out += "\n";
	out += "########################################################\n";
	out += "#  Options from configuration file ";
	out += confFile;
	out += ":\n";
	out += "########################################################\n";
	for (int i=0; i<opts.size(); i++) {
		// MOD HERE
		//if (source[i] == "f") {
		for (int j=0; j<source[i].size(); j++) {
			if (source[i][j] == "f") {
				//out += opts[i];
				string line = opts[i];
				while (line.size() < 30) {
					line += " ";
				}
				//out += "    ";
				// MOD HERE
				//out += vals[i];
				line += vals[i][j];
				line += "\n";
				out += line;
			}
		}
	}
	return out;
}



int OptionParser::countOptions() {
	return opts.size();

}

/*********** PUBLIC FUNCTIOS TO GET ARGUMENTS ****************/

string OptionParser::getCommandName() const {
	return commandName;
}

string OptionParser::getEnv(string _env) {
	errorFlag = false;
	if (getenv(_env.c_str()) != NULL) {
		return getenv(_env.c_str());
	} else {
		errorFlag = true;
		return "";
	}
}

string OptionParser::getString(string name) {
	return getString(name, 0);
}

string OptionParser::getString(string name, int pos) {
	errorFlag = false;
	vector<string> tmp = getMultiString(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return "";
	} else {
		return tmp[pos];
	}
}

vector<string> OptionParser::getMultiString(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<string> out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(getSingleString(name, i));
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

int OptionParser::getInt(string name) {
	return getInt(name, 0);
}

int OptionParser::getInt(string name, int pos) {

	errorFlag = false;
	vector<int> tmp = getMultiInt(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return 0;
	} else {
		return tmp[pos];
	}
}

vector<int> OptionParser::getMultiInt(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<int> out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(atoi(getSingleString(name, i).c_str()));
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

unsigned int OptionParser::getUnsignedInt(string name) {
	return getInt(name, 0);
}

unsigned int OptionParser::getUnsignedInt(string name, int pos) {

	errorFlag = false;
	vector<unsigned int> tmp = getMultiUnsignedInt(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return 0;
	} else {
		return tmp[pos];
	}
}
vector<unsigned int> OptionParser::getMultiUnsignedInt(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<unsigned int> out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back((unsigned)(atoi(getSingleString(name, i).c_str())));
		if (errorFlag) {
			return out;
		}
	}
	return out;
}



double OptionParser::getDouble(string name) {
	return getDouble(name, 0);
}

double OptionParser::getDouble(string name, int pos) {
	errorFlag = false;
	vector<double> tmp = getMultiDouble(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return 0.0;
	} else {
		return tmp[pos];
	}
}

vector<double> OptionParser::getMultiDouble(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<double> out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(atof(getSingleString(name, i).c_str()));
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

bool OptionParser::getBool(string name) {
	return getBool(name, 0);
}

bool OptionParser::getBool(string name, int pos) {
	errorFlag = false;
	vector<bool> tmp = getMultiBool(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return false;
	} else {
		return tmp[pos];
	}
}
	
vector<bool> OptionParser::getMultiBool(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<bool> out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		string line = getSingleString(name, i);
		if (errorFlag) {
			return out;
		}

		if (line.length() == 0) {
			out.push_back(true);
			continue;
		}
		for (int i=0; i<line.length(); i++) {
			line[i] = toupper(line[i]);
		}
		if (line == "TRUE" || line == "T" || line == "1" || line == "Y" || line == "YES") {
			out.push_back(true);
			continue;
		}
		if (line == "FALSE" || line == "F" || line == "0" || line == "N" || line == "NO") {
			out.push_back(false);
			continue;
		}
		errorFlag = true;
		out.clear();
		return out;
	}
	return out;
}

vector<string> OptionParser::getStringVector(string name) {
	return getStringVector(name, 0);
}
vector<string> OptionParser::getStringVector(string name, int pos) {
	errorFlag = false;
	vector<vector<string> > tmp = getMultiStringVector(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return vector<string>();
	} else {
		return tmp[pos];
	}
}

vector<string> OptionParser::getStringVectorJoinAll(string name) {
	errorFlag = false;
	vector<vector<string> > tmp = getMultiStringVector(name);
	vector<string> out;
	if (tmp.size() == 0) {
		errorFlag = true;
	} else {
		for (unsigned int i=0; i<tmp.size(); i++) {
			out.insert(out.end(), tmp[i].begin(), tmp[i].end());
		}
	}
	return out;
}

vector<vector<string> > OptionParser::getMultiStringVector(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<vector<string> > out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(getArrayOfStrings(name, i));
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

vector<int> OptionParser::getIntVector(string name) {
	return getIntVector(name, 0);
}

vector<int> OptionParser::getIntVector(string name, int pos) {
	errorFlag = false;
	vector<vector<int> > tmp = getMultiIntVector(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return vector<int>();
	} else {
		return tmp[pos];
	}
}

vector<unsigned int> OptionParser::getUnsignedIntVector(string name) {
	return getUnsignedIntVector(name, 0);
}

vector<unsigned int> OptionParser::getUnsignedIntVector(string name, int pos) {
	errorFlag = false;
	vector<vector<unsigned int> > tmp = getMultiUnsignedIntVector(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return vector<unsigned int>();
	} else {
		return tmp[pos];
	}
}

vector<int> OptionParser::getIntVectorJoinAll(string name) {
	errorFlag = false;
	vector<vector<int> > tmp = getMultiIntVector(name);
	vector<int> out;
	if (tmp.size() == 0) {
		errorFlag = true;
	} else {
		for (unsigned int i=0; i<tmp.size(); i++) {
			out.insert(out.end(), tmp[i].begin(), tmp[i].end());
		}
	}
	return out;
}

vector<unsigned int> OptionParser::getUnsignedIntVectorJoinAll(string name) {
	errorFlag = false;
	vector<vector<unsigned int> > tmp = getMultiUnsignedIntVector(name);
	vector<unsigned int> out;
	if (tmp.size() == 0) {
		errorFlag = true;
	} else {
		for (unsigned int i=0; i<tmp.size(); i++) {
			out.insert(out.end(), tmp[i].begin(), tmp[i].end());
		}
	}
	return out;
}

vector<vector<int> > OptionParser::getMultiIntVector(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<vector<int> > out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(vector<int>());
		vector<string> vecs = getArrayOfStrings(name, i);
		for (int j=0; j<vecs.size(); j++) {
			out[i].push_back(atoi(vecs[j].c_str()));
		}
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

vector<vector<unsigned int> > OptionParser::getMultiUnsignedIntVector(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<vector<unsigned int> > out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(vector<unsigned int>());
		vector<string> vecs = getArrayOfStrings(name, i);
		for (int j=0; j<vecs.size(); j++) {
			out[i].push_back(atoi(vecs[j].c_str()));
		}
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

vector<double> OptionParser::getDoubleVector(string name) {
	return getDoubleVector(name, 0);
}

vector<double> OptionParser::getDoubleVector(string name, int pos) {
	errorFlag = false;
	vector<vector<double> > tmp = getMultiDoubleVector(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return vector<double>();
	} else {
		return tmp[pos];
	}
}

vector<double> OptionParser::getDoubleVectorJoinAll(string name) {
	errorFlag = false;
	vector<vector<double> > tmp = getMultiDoubleVector(name);
	vector<double> out;
	if (tmp.size() == 0) {
		errorFlag = true;
	} else {
		for (unsigned int i=0; i<tmp.size(); i++) {
			out.insert(out.end(), tmp[i].begin(), tmp[i].end());
		}
	}
	return out;
}

vector<vector<double> > OptionParser::getMultiDoubleVector(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<vector<double> > out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		out.push_back(vector<double>());
		vector<string> vecs = getArrayOfStrings(name, i);
		for (int j=0; j<vecs.size(); j++) {
			out[i].push_back(atof(vecs[j].c_str()));
		}
		if (errorFlag) {
			return out;
		}
	}
	return out;
}

vector<bool> OptionParser::getBoolVector(string name) {
	return getBoolVector(name, 0);
}

vector<bool> OptionParser::getBoolVector(string name, int pos) {
	errorFlag = false;
	vector<vector<bool> > tmp = getMultiBoolVector(name);
	if (errorFlag || tmp.size() <= pos) {
		errorFlag = true;
		return vector<bool>();
	} else {
		return tmp[pos];
	}
}

vector<bool> OptionParser::getBoolVectorJoinAll(string name) {
	errorFlag = false;
	vector<vector<bool> > tmp = getMultiBoolVector(name);
	vector<bool> out;
	if (tmp.size() == 0) {
		errorFlag = true;
	} else {
		for (unsigned int i=0; i<tmp.size(); i++) {
			out.insert(out.end(), tmp[i].begin(), tmp[i].end());
		}
	}
	return out;
}

vector<vector<bool> > OptionParser::getMultiBoolVector(string name) {
	errorFlag = false;
	int number = getOptionNumberOfMultiples(name);
	vector<vector<bool> > out;
	if (number == 0) {
		errorFlag = true;
		return out;
	}
	for (int i=0; i<number; i++) {
		vector<bool> vecb;
		vector<string> vecs = getArrayOfStrings(name, i);
		for (int k=0; k<vecs.size(); k++) {
			for (int j=0; j<(vecs[k]).length(); j++) {
				(vecs[k])[j] = toupper((vecs[k])[j]);
			}
			if (vecs[k] == "TRUE" || vecs[k] == "T" || vecs[k] == "1" || vecs[k] == "Y" || vecs[k] == "YES" ) {
				vecb.push_back(true);
			} else {
				if (vecs[k] == "FALSE" || vecs[k] == "F" || vecs[k] == "0" || vecs[k] == "N" || vecs[k] == "NO" ) {
					vecb.push_back(false);
				} else {
					vecb.clear();
					errorFlag = true;
					out.push_back(vecb);
					return out;
				}
			}
		}
		out.push_back(vecb);
		if (errorFlag) {
			return out;
		}
	}
	return out;
}


bool OptionParser::fail() {
	return errorFlag;
}



/************ REQUIRED AND ALLOWED OPTIONS ************************/

void OptionParser::setRequired(vector<string> _requiredOptions) {
	disallowed.clear();
	missing.clear();
	ambiguous.clear();
	required = _requiredOptions;
}

void OptionParser::setAllowed(vector<string> _allowedOptions) {
	disallowed.clear();
	missing.clear();
	ambiguous.clear();
	allowed = _allowedOptions;
}

void OptionParser::setMutualExclusive(vector<vector<string> > _mutualExclusiveOptions) {
	disallowedTogether.clear();
	mutualExclusive = _mutualExclusiveOptions;

}

void OptionParser::setDependentOn(vector<vector<string> > _dependentOptions) {
	missingDependency.clear();
	dependentOn = _dependentOptions;
}

void OptionParser::setInterDependent(vector<vector<string> > _interDependentOptions) {
	missingInterdipendenty.clear();
	interDependent = _interDependentOptions;
}

void OptionParser::setLinked(vector<vector<string> > _linkedOptions) {
	missingLinked.clear();
	linked = _linkedOptions;
}

void OptionParser::setOneRequired(vector<vector<string> > _oneRequiredOptions) {
	missingOneRequired.clear();
	oneRequired = _oneRequiredOptions;
}

void OptionParser::autoExtendOptions() {

	/***************************************************
	 *  If the required and/or allowed options are set,
	 *  if an option is given that matches part of one
	 *  of the allowed/required options, that options is
	 *  auto-extended
	 *
	 *  For example, the option "die" may be used to 
	 *  specify "dielectric"
	 * 
	 *  If the short options is ambiguous then it is not
	 *  extended
	 *
	 *  For example, if "use" was given and two valid
	 *  options were "useVdw" and "useElec" the options
	 *  is not extended 
	 *  
	 ***************************************************/
	
	for (int i=0; i<opts.size(); i++) {
		bool found = false;
		int index = -1;
		string inVector = "";
		bool foundWholeOption = false;
		bool foundAmbiguous = false;
		for (int j=0; j<required.size(); j++) {
			if (opts[i] == required[j]) {
				foundWholeOption = true;
				break;
			}
			if (opts[i].size() < required[j].size()) {
				if (required[j].substr(0,opts[i].size()) == opts[i]) {
					if (found == true) {
						// the short option fits two possibilities  
						ambiguous.push_back(opts[i]);
						found = false;
						foundAmbiguous = true;
						break;
					} else {
						found = true;
						index = j;
						inVector = "required";
					}
				}
			}
		}
		if (foundWholeOption || foundAmbiguous) {
			continue;
		}
		for (int j=0; j<allowed.size(); j++) {
			if (opts[i] == allowed[j]) {
				foundWholeOption = true;
				break;
			}
			if (opts[i].size() < allowed[j].size()) {
				if (allowed[j].substr(0,opts[i].size()) == opts[i]) {
					if (found == true) {
						// the short option fits two possibilities  
						ambiguous.push_back(opts[i]);
						found = false;
						foundAmbiguous = true;
						break;
					} else {
						found = true;
						index = j;
						inVector = "allowed";
					}
				}
			}
		}
		if (found && !foundWholeOption && !foundAmbiguous) {
			if (inVector == "required") {
				opts[i] = required[index];
			}
			if (inVector == "allowed") {
				opts[i] = allowed[index];
			}
		}
	}
}

bool OptionParser::checkOptions() {
	disallowed.clear();
	missing.clear();
	disallowedTogether.clear();
	missingDependency.clear();
	missingInterdipendenty.clear();
	missingOneRequired.clear();

	bool out = true;

	vector<string> amb = getAmbiguousOptions();
	if (amb.size() > 0) {
		out = false;
	}

	// all the required options need to be present in
	// the opts vector
	for (int i=0; i<required.size(); i++) {
		bool status = false;
		for (int j=0; j<opts.size(); j++) {
			if (required[i] == opts[j]) {
				status = true;
				break;
			}
		}
		if (!status) {
			missing.push_back(required[i]);
			out = false;
		}
	}

	// all the opts need to be present in
	// the allowed or required vector
	if (allowed.size() > 0) {
		for (int j=0; j<opts.size(); j++) {
			bool status = false;
			for (int i=0; i<allowed.size(); i++) {
				if (allowed[i] == opts[j]) {
					status = true;
					break;
				}
			}
			if (!status) {
				// if it is not in the allowed list it might
				// be in the required list
				for (int i=0; i<required.size(); i++) {
					if (required[i] == opts[j]) {
						status = true;
						break;
					}
				}
			}
			if (!status) {
				disallowed.push_back(opts[j]);
				out = false;
			}
		}
	}

	// options that cannot be given together
	for (unsigned int i=0; i<mutualExclusive.size(); i++) {
		unsigned int foundCounter = 0;
		vector<string> offending;
		for (unsigned int j=0; j<mutualExclusive[i].size(); j++) {
			for (int k=0; k<opts.size(); k++) {
				string line = opts[k];
				line = MslTools::toUpper(line);
				if (opts[k] == mutualExclusive[i][j]) {
					if (!getBool(opts[k])) {
						// mutually exclusive options can be given only if they are false
						// NOT THE RIGHT WAY: CREATE NEW CATEGORY: boolMutuallyExclusive
						continue;
					}
					foundCounter++;
					offending.push_back(mutualExclusive[i][j]);
					break;
				}
			}
		}
		if (foundCounter > 1) {
			disallowedTogether.push_back(offending);
			out = false;
		}
	}

	// one of these is required
	for (unsigned int i=0; i<oneRequired.size(); i++) {
		unsigned int foundCounter = 0;
		vector<string> offending;
		for (unsigned int j=0; j<oneRequired[i].size(); j++) {
			offending.push_back(oneRequired[i][j]);
			for (int k=0; k<opts.size(); k++) {
				if (opts[k] == oneRequired[i][j]) {
					foundCounter++;
					break;
				}
			}
		}
		if (foundCounter == 0) {
			missingOneRequired.push_back(offending);
			out = false;
		}
	}

	// options that depend on the presence of a specific option
	for (unsigned int i=0; i<dependentOn.size(); i++) {
		bool dependencyFound = false;
		bool dependentFound = false;
		vector<string> offending;
		offending.push_back(dependentOn[i][0]);
		for (int k=0; k<opts.size(); k++) {
			if (dependentOn[i][0] == opts[k]) {
				dependencyFound = true;
				break;
			}
		}
		if (!dependencyFound) {
			for (unsigned int j=1; j<dependentOn[i].size(); j++) {
				for (int k=0; k<opts.size(); k++) {
					if (dependentOn[i][j] == opts[k]) {
						dependentFound = true;
						offending.push_back(dependentOn[i][j]);
						break;
					}
				}
			}
		}
		if (dependentFound && !dependencyFound) {
			out = false;
			missingDependency.push_back(offending);
		}
	}

	for (unsigned int i=0; i<interDependent.size(); i++) {
		unsigned int foundCounter = 0;
		vector<string> offending = interDependent[i];
		for (unsigned int j=0; j<interDependent[i].size(); j++) {
			for (int k=0; k<opts.size(); k++) {
				if (opts[k] == interDependent[i][j]) {
					foundCounter++;
					break;
				}
			}
		}
		if (foundCounter != 0 && foundCounter != interDependent[i].size()) {
			missingInterdipendenty.push_back(offending);
			out = false;
		}
	}

	for (unsigned int i=0; i<linked.size(); i++) {
		unsigned int foundCounter = 0;
		bool differentCount = false;
		unsigned int valCounter = 0;
		vector<string> offending = linked[i];
		for (unsigned int j=0; j<linked[i].size(); j++) {
			for (int k=0; k<opts.size(); k++) {
				if (opts[k] == linked[i][j]) {
					foundCounter++;
					if (j==0) {
						valCounter = vals[k].size();
					} else {
						if (valCounter != vals[k].size()) {
							differentCount = true;
						}
					}
					break;
				}
			}
		}
		if ((foundCounter != 0 && foundCounter != linked[i].size()) || differentCount) {
			missingLinked.push_back(offending);
			out = false;
		}
	}

	return out;
}

/*
vector<string> OptionParser::getMissingOptions() const {
	return missing;
}

vector<string> OptionParser::getDisallowedOptions() const {
	return disallowed;
}

vector<string> OptionParser::getAmbiguousOptions() const {
	return ambiguous;
}
*/


/************ PRIVATE FUNCTIONS ***********************************/

vector<string> OptionParser::getArrayOfStrings(string name, int index) {
	errorFlag = false;
	vector<string> out;
	for (unsigned int i=0; i<opts.size(); i++) {
		if (opts[i] == name) {
			out = MslTools::tokenize(vals[i][index]);
			if (translateEnv_flag) {
				for (unsigned int j=0; j<out.size(); j++) {
					out[j] = translateEnv(out[j]);
				}
			}
			return out;
		}
	}
	errorFlag = true;
	return out;
}

string OptionParser::getSingleString(string name, int index) {
	errorFlag = false;
	string out = "";
	for (int i=0; i<opts.size(); i++) {
		if (opts[i] == name) {
			if (translateEnv_flag) {
				return (translateEnv(vals[i][index]));
			} else {
				return vals[i][index];
			}
		}
	}
	errorFlag = true;
	return out;
}

int OptionParser::getOptionNumberOfMultiples(string name) {
	errorFlag = false;
	for (int i=0; i<opts.size(); i++) {
		if (opts[i] == name) {
			return vals[i].size();
		}
	}
	errorFlag = true;
	return 0;
}

string OptionParser::translateEnv(string _input) {
	// translate $VARIABLES with their corresponding environmental variable
	//
	// TO DO: make it so that it doesn't translate if in single quotes
	
	string out = "";
	for (unsigned int i=0; i<_input.size(); i++) {
		if (int(_input[i]) == 36) {
			string envVar = "";
			while (i<_input.size()-1) {
				/************************************
				 *  Allowed: letters, number and underscore "_"
				 * 
				 *  ascii ranges
				 *  A-Z = 65 to 90
				 *  a-z = 97 to 122
				 *  0-9 = 48 to 57
				 *  _ = 97 
				 *************************************/
				if ((int(_input[i+1]) >= 65 && int(_input[i+1]) <= 90) || (int(_input[i+1]) >= 97 && int(_input[i+1]) <= 122) || (int(_input[i+1]) >= 48 && int(_input[i+1]) <= 57) || int(_input[i+1]) == 95) {
					i++;
					envVar += _input[i];
				} else {
					break;
				}
			}
			out += getEnv(envVar);
		} else {
			out += _input[i];
		}
	}
	return out;
}

void OptionParser::setShortOptionEquivalent(vector<vector<string> > _equivalent) {
	for (unsigned int i=0; i<_equivalent.size(); i++) {
		if (_equivalent[i].size() != 2) {
			cerr << "ERROR 5808: not 2 string in entry "<< i << " in void setShortOptionEquivalent(vector<vector<string> > _equivalent)" << endl;
			exit(5808);
		}
		setShortOptionEquivalent(_equivalent[i][0], _equivalent[i][1]);
	}
}

void OptionParser::setShortOptionEquivalent(string _oneLetter, string _longOption) {
	if (_oneLetter.size() != 1) {
		cerr << "ERROR 5812: short option " << _oneLetter << " is not a one letter option in void OptionParser::setShortOptionEquivalent(string _oneLetter, string _longOption)" << endl;
		exit(5812);
	}
	if (_longOption.size() == 0) {
		cerr << "ERROR 5817: empty long option in void OptionParser::setShortOptionEquivalent(string _oneLetter, string _longOption)" << endl;
		exit(5817);
	}
	vector<string> tmp;
	tmp.push_back(_oneLetter);
	tmp.push_back(_longOption);
	shortOptEquivalent.push_back(tmp);
	linkShortOptions();
}

void OptionParser::linkShortOptions() {

	/*********************************************************
	 *  Take all the options that are given on command line as 
	 *  short options (like -h) and substitute them to their linked
	 *  long option (like help) if it exists
	 *********************************************************/
	for (int i=0; i<opts.size(); i++) {
		if (opts[i].size() != 1) {
			// not a 1 letter option
			continue;
		}
		bool found = false;
		string longOpt = "";
		for (unsigned int j=0; j<shortOptEquivalent.size(); j++) {
			if (opts[i] == shortOptEquivalent[j][0]) {
				// found that the short option is linked
				found = true;
				longOpt = shortOptEquivalent[j][1];
			}
		}
		if (!found) {
			continue;
		}

		// look if the linked option was given
		found = false;
		for (int j=0; j<opts.size(); j++) {
			if (i == j) {
				continue;
			}
			if (longOpt == opts[j]) {
				found = true;
				// select only the values given as short option (-h -> source = "o") 
				// not long (--h -> "a") or file ("f")
				for (int k=0; k<source[i].size(); k++) {
					if (source[i][k] == "o") {
						// add them to opts[j]
						vals[j].push_back(vals[i][k]);
						source[j].push_back("a");

						// delete the value from opts[i]
						vals[i].erase(vals[i].begin()+k);
						source[i].erase(source[i].begin()+k);

						// if that was the only value for the option erase the option completely
						if (vals[i].size() == 0) {
							opts.erase(opts.begin()+i);
							vals.erase(vals.begin()+i);
							source.erase(source.begin()+i);
						}
						
					}
				}
				break;
			}
		}
		if (!found) {
			opts.push_back(longOpt);
			vals.push_back(vector<string>());
			source.push_back(vector<string>());
			for (int k=0; k<source[i].size(); k++) {
				if (source[i][k] == "o") {
					vals.back().push_back(vals[i][k]);
					source.back().push_back("a");

					// delete the value from opts[i]
					vals[i].erase(vals[i].begin()+k);
					source[i].erase(source[i].begin()+k);

					// if that was the only value for the option erase the option completely
					if (vals[i].size() == 0) {
						opts.erase(opts.begin()+i);
						vals.erase(vals.begin()+i);
						source.erase(source.begin()+i);
					}
				}
			}
		}

	}
}

string OptionParser::getErrors() const {
	string out;
	for (unsigned int i=0; i<missing.size(); i++) {
		out += "ERROR: missing required option: \"" + missing[i] + "\"\n";
	}
	for (unsigned int i=0; i<disallowed.size(); i++) {
		out += "ERROR: unknown option: \"" + disallowed[i] + "\"\n";
	}
	if (ambiguous.size() > 0) {
		for (unsigned int i=0; i<ambiguous.size(); i++) {
			out += "ERROR: ambiguous option: \"" + ambiguous[i] + "\"\n";
		}
	}
	for (unsigned int i=0; i<missingOneRequired.size(); i++) {
		out += "ERROR: missing one required options among the following:";
		for (unsigned int j=0; j<missingOneRequired[i].size(); j++) {
			out += (string)" \"" + missingOneRequired[i][j] + (string)"\"";
		}
		out += "\n";
	}
	for (unsigned int i=0; i<disallowedTogether.size(); i++) {
		out += "ERROR: incompatible options:";
		for (unsigned int j=0; j<disallowedTogether[i].size(); j++) {
			out += (string)" \"" + disallowedTogether[i][j] + (string)"\"";
		}
		out += "\n";
	}
	for (unsigned int i=0; i<missingDependency.size(); i++) {
		out += "ERROR: options: ";
		for (unsigned int j=1; j<missingDependency[i].size(); j++) {
			out += (string)" \"" + missingDependency[i][j] + (string)"\"";
		}
		out +=  "; depend on the missing option: \"" + missingDependency[i][0] + (string)"\"\n";
	}
	for (unsigned int i=0; i<missingInterdipendenty.size(); i++) {
		out += "ERROR: missing one or more of these interdependent options:";
		for (unsigned int j=0; j<missingInterdipendenty[i].size(); j++) {
			out += (string)" \"" + missingInterdipendenty[i][j] + (string)"\"";
		}
		out += "\n";
	}
	for (unsigned int i=0; i<missingLinked.size(); i++) {
		out += "ERROR: linked options should be given the same number of times:";
		for (unsigned int j=0; j<missingLinked[i].size(); j++) {
			out += (string)" \"" + missingLinked[i][j] + (string)"\"";
		}
		out += "\n";
	}
	for (unsigned int i=0; i<errorMessages.size(); i++) {
		out += errorMessages[i]  + (string)"\n";
	}
	return out;
}

string OptionParser::getWarnings() const {
	string out;
	for (unsigned int i=0; i<warningMessages.size(); i++) {
		out += warningMessages[i]  + (string)"\n";
	}
	return out;
}

