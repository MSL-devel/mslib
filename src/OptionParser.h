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

#ifndef ARGUMENT_H
#define ARGUMENT_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include "MslTools.h"

using namespace std;

/*! \brief Object that parses command line and/or configuration file options
 *
 *  command line:
 *    % command --opt1 bar --opt2 2.343 --opt3 --opt4 true
 *
 *  file:
 *    opt1    bar
 *    opt2    2.343
 *    opt3           # boolean, true with no val
 *    opt4    true   # takes value true, false, 0, 1
 *
 *  THE COMMAND LINE OPTIONS OVERRIDE THE CONF FILE VALUES!
 *
 *  To do: stuff in quote should be parse as unit
 *  
 */

class OptionParser {

	public:
		OptionParser();
		OptionParser(const OptionParser & _OP);
		~OptionParser();

		void operator=(const OptionParser & _OP);

		/****************************************************
		 *  INITIALIZERS WITH DATA FROM A FILE OR THE COMMAND LINE ARGUMENTS 
		 *
		 *  They can both be given and the command line superseeds the
		 *  configuration file
		 ****************************************************/
		bool readFile(string filename);
		void readArgv(int theArgc, char *theArgv[]);

		/****************************************************
		 *  FUNCTION TO GET THE VALUE ARGUMENT GIVEN ITS NAME
		 *
		 *  The data can be parsed to give
		 *   - single int, or double, or string, or bool
		 *   - a vector of int, double, string or bool
		 *
		 *  If the argument was not found an empty value is returned
		 *  (int=0, double=0.0, string="", bool=false, vectors all empty)
		 *  and the fail() function will return true.
		 *
		 *  Again, the command line superseeds the configuration file
		 ****************************************************/
		int getInt(string name);
		double getDouble(string name);
		string getString(string name);
		bool getBool(string name);

		int getInt(string name, int pos);
		double getDouble(string name, int pos);
		string getString(string name, int pos);
		bool getBool(string name, int pos);

		vector<int> getIntVector(string name, int pos);
		vector<double> getDoubleVector(string name, int pos);
		vector<string> getStringVector(string name, int pos);
		vector<bool> getBoolVector(string name, int pos);

		vector<int> getIntVector(string name);
		vector<double> getDoubleVector(string name);
		vector<string> getStringVector(string name);
		vector<bool> getBoolVector(string name);

		/* FOR REPEATED OPTIONS: "-opt ARG1 -opt ARG2" RETURN AS vector<type> (ARG1, ARG2) */
		vector<int> getMultiInt(string name);
		vector<double> getMultiDouble(string name);
		vector<string> getMultiString(string name);
		vector<bool> getMultiBool(string name);

		/* FOR MULTIPLE REPEATED OPTIONS: "-opt ARG1 ARG2 -opt ARG3 ARG4" RETURN AS vector<vector<type> > ((ARG1, ARG2), (ARG3, ARG3)) */
		vector<vector<int> > getMultiIntVector(string name);
		vector<vector<double> > getMultiDoubleVector(string name);
		vector<vector<string> > getMultiStringVector(string name);
		vector<vector<bool> > getMultiBoolVector(string name);
		
		/* JOIN REPEATED OPTIONS: "--opt ARG1 ARG3 --opt ARG3 ARG4" into a single vetcor as (ARG1, ARG2, ARG3, ARG4) */
		vector<int> getIntVectorJoinAll(string name);
		vector<double> getDoubleVectorJoinAll(string name);
		vector<string> getStringVectorJoinAll(string name);
		vector<bool> getBoolVectorJoinAll(string name);

		bool fail();
		

		// ADD A "MUTUALLY EXCLUSIVE BUT ONE REQUIRED" OPTION?


		void setRequired(vector<string> _requiredOptions); // those that can and have to be there
		void setAllowed(vector<string> _allowedOptions); // those that can be there but don't have to
		void setMutualExclusive(vector<vector<string> > _mutualExclusiveOptions); // those that cannot be given together
		void setDependentOn(vector<vector<string> > _dependentOptions); // if the first is not given, the other should not too
		void setInterDependent(vector<vector<string> > _interDependentOptions); // all given together or none
		void setLinked(vector<vector<string> > _linkedOptions); // these are interdependent but also require to be given the same number of times
		void setOneRequired(vector<vector<string> > _oneRequiredOptions); // one of these options must be given

		void setDefaultArgument(string _argument);
		void setDefaultArguments(vector<string> _arguments);
		bool checkOptions();
		vector<string> getMissingOptions() const {return missing;};
		vector<string> getDisallowedOptions() const {return disallowed;};
		vector<vector<string> > getDisallowedTogetherOptions() const {return disallowedTogether;};
		vector<vector<string> > getMissingDependencyOptions() const {return missingDependency;};
		vector<vector<string> > getMissingInterdependencyOptions() const {return missingInterdipendenty;};
		vector<vector<string> > getMissingLinkedOptions() const {return missingLinked;};
		vector<vector<string> > getMissingOneRequiredOptions() const {return missingOneRequired;};
		vector<string> getAmbiguousOptions() const {return ambiguous;};

		string getCommandName() const;
		string getEnv(string _env);
		int countOptions();
		void printConfFile() const; // print a configuration file 
		bool writeConfFile(string filename); // write a conf file
		string getConfFile() const; // get a string with the conf file
//		void printOptions() const; // write a configuration file 
//		bool writeOptions(string filename);

		void setAutoExtend(bool flag);
		bool getAutoExtend() const;		

		void autoExtendOptions();

		void setShortOptionEquivalent(vector<vector<string> >);
		void setShortOptionEquivalent(string _oneLetter, string _longOption);
		void linkShortOptions();
		string getErrors() const;
		string getWarnings() const;
		void addError(string _error);
		void addWarning(string _warning);

		friend ofstream & operator<<(ofstream &_of, OptionParser & _opts) { _of << _opts.getConfFile(); return _of; }
		friend ostream  & operator<<(ostream  &_os, OptionParser & _opts) { _os << _opts.getConfFile(); return _os; }

	private:
		void copy(const OptionParser & _OP);
		string getSingleString(string name, int index);
		vector<string> getArrayOfStrings(string name, int index);
		int getOptionNumberOfMultiples(string name);
		string translateEnv(string _input);
		void addFreeArgumentsToDefault();

		vector<string> errorMessages;
		vector<string> warningMessages;


		//string makeRunConfFile() const;

		/***********************************************
		 *  Internal storage of the options, 3 string vectors:
		 *
		 *  opts      the options
		 *  vals      the values of the options
		 *  source    source of the value is "a" for argument or "f" for file
		 ***********************************************/
		vector<string> opts;
		vector<vector<string> > vals;
		vector<vector<string> > source;
		vector<string> freeArgs;

		vector<vector<string> > shortOptEquivalent;

		string commandName;

		string confFile;

		vector<string> required;
		vector<string> allowed;
		vector<vector<string> > mutualExclusive;
		vector<vector<string> > dependentOn;
		vector<vector<string> > interDependent;
		vector<vector<string> > linked;
		vector<vector<string> > oneRequired;

		vector<string> missing;
		vector<string> disallowed;
		vector<string> ambiguous;
		vector<vector<string> > disallowedTogether;
		vector<vector<string> > missingDependency;
		vector<vector<string> > missingInterdipendenty;
		vector<vector<string> > missingLinked;
		vector<vector<string> > missingOneRequired;

		vector<string> defaultArguments;
		
		bool errorFlag;
		bool translateEnv_flag;	
};
inline void OptionParser::addError(string _error) {
	errorMessages.push_back(_error);
}
inline void OptionParser::addWarning(string _warning) {
	warningMessages.push_back(_warning);
}
#endif
