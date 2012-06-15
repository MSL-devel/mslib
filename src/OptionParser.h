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

#ifndef ARGUMENT_H
#define ARGUMENT_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include "MslTools.h"


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

namespace MSL { 
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
		bool readFile(std::string filename);
		void readArgv(int theArgc, char *theArgv[]);

		/****************************************************
		 *  FUNCTION TO GET THE VALUE ARGUMENT GIVEN ITS NAME
		 *
		 *  The data can be parsed to give
		 *   - single int, or double, or string, or bool
		 *   - a std::vector of int, double, std::string or bool
		 *
		 *  If the argument was not found an empty value is returned
		 *  (int=0, double=0.0, std::string="", bool=false, vectors all empty)
		 *  and the fail() function will return true.
		 *
		 *  Again, the command line superseeds the configuration file
		 ****************************************************/
		int getInt(std::string name);
		unsigned int getUnsignedInt(std::string name);
		double getDouble(std::string name);
		std::string getString(std::string name);
		bool getBool(std::string name);

		int getInt(std::string name, int pos);
		unsigned int getUnsignedInt(std::string name, int pos);
		double getDouble(std::string name, int pos);
		std::string getString(std::string name, int pos);
		bool getBool(std::string name, int pos);

		std::vector<int> getIntVector(std::string name, int pos);
		std::vector<unsigned int> getUnsignedIntVector(std::string name, int pos);
		std::vector<double> getDoubleVector(std::string name, int pos);
		std::vector<std::string> getStringVector(std::string name, int pos);
		std::vector<bool> getBoolVector(std::string name, int pos);

		std::vector<int> getIntVector(std::string name);
		std::vector<unsigned int> getUnsignedIntVector(std::string name);
		std::vector<double> getDoubleVector(std::string name);
		std::vector<std::string> getStringVector(std::string name);
		std::vector<bool> getBoolVector(std::string name);

		/* FOR REPEATED OPTIONS: "-opt ARG1 -opt ARG2" RETURN AS std::vector<type> (ARG1, ARG2) */
		std::vector<int> getMultiInt(std::string name);
		std::vector<unsigned int> getMultiUnsignedInt(std::string name);
		std::vector<double> getMultiDouble(std::string name);
		std::vector<std::string> getMultiString(std::string name);
		std::vector<bool> getMultiBool(std::string name);

		/* FOR MULTIPLE REPEATED OPTIONS: "-opt ARG1 ARG2 -opt ARG3 ARG4" RETURN AS std::vector<std::vector<type> > ((ARG1, ARG2), (ARG3, ARG3)) */
		std::vector<std::vector<int> > getMultiIntVector(std::string name);
		std::vector<std::vector<unsigned int> > getMultiUnsignedIntVector(std::string name);
		std::vector<std::vector<double> > getMultiDoubleVector(std::string name);
		std::vector<std::vector<std::string> > getMultiStringVector(std::string name);
		std::vector<std::vector<bool> > getMultiBoolVector(std::string name);
		
		/* JOIN REPEATED OPTIONS: "--opt ARG1 ARG3 --opt ARG3 ARG4" into a single vetcor as (ARG1, ARG2, ARG3, ARG4) */
		std::vector<int> getIntVectorJoinAll(std::string name);
		std::vector<unsigned int> getUnsignedIntVectorJoinAll(std::string name);
		std::vector<double> getDoubleVectorJoinAll(std::string name);
		std::vector<std::string> getStringVectorJoinAll(std::string name);
		std::vector<bool> getBoolVectorJoinAll(std::string name);

		bool fail();
		

		// ADD A "MUTUALLY EXCLUSIVE BUT ONE REQUIRED" OPTION?


		void setRequired(std::vector<std::string> _requiredOptions); // those that can and have to be there
		void setAllowed(std::vector<std::string> _allowedOptions); // those that can be there but don't have to
		void setMutualExclusive(std::vector<std::vector<std::string> > _mutualExclusiveOptions); // those that cannot be given together
		void setDependentOn(std::vector<std::vector<std::string> > _dependentOptions); // if the first is not given, the other should not too
		void setInterDependent(std::vector<std::vector<std::string> > _interDependentOptions); // all given together or none
		void setLinked(std::vector<std::vector<std::string> > _linkedOptions); // these are interdependent but also require to be given the same number of times
		void setOneRequired(std::vector<std::vector<std::string> > _oneRequiredOptions); // one of these options must be given

		void setDefaultArgument(std::string _argument);
		void setDefaultArguments(std::vector<std::string> _arguments);
		bool checkOptions();
		std::vector<std::string> getMissingOptions() const {return missing;};
		std::vector<std::string> getDisallowedOptions() const {return disallowed;};
		std::vector<std::vector<std::string> > getDisallowedTogetherOptions() const {return disallowedTogether;};
		std::vector<std::vector<std::string> > getMissingDependencyOptions() const {return missingDependency;};
		std::vector<std::vector<std::string> > getMissingInterdependencyOptions() const {return missingInterdipendenty;};
		std::vector<std::vector<std::string> > getMissingLinkedOptions() const {return missingLinked;};
		std::vector<std::vector<std::string> > getMissingOneRequiredOptions() const {return missingOneRequired;};
		std::vector<std::string> getAmbiguousOptions() const {return ambiguous;};

		std::string getCommandName() const;
		std::string getEnv(std::string _env);
		int countOptions();
		void printConfFile() const; // print a configuration file 
		bool writeConfFile(std::string filename); // write a conf file
		std::string getConfFile() const; // get a std::string with the conf file
//		void printOptions() const; // write a configuration file 
//		bool writeOptions(std::string filename);

		void setAutoExtend(bool flag);
		bool getAutoExtend() const;		

		void autoExtendOptions();

		void setShortOptionEquivalent(std::vector<std::vector<std::string> >);
		void setShortOptionEquivalent(std::string _oneLetter, std::string _longOption);
		void linkShortOptions();
		std::string getErrors() const;
		std::string getWarnings() const;
		void addError(std::string _error);
		void addWarning(std::string _warning);

		friend std::ofstream & operator<<(std::ofstream &_of, OptionParser & _opts) { _of << _opts.getConfFile(); return _of; }
		friend std::ostream  & operator<<(std::ostream  &_os, OptionParser & _opts) { _os << _opts.getConfFile(); return _os; }

	private:
		void copy(const OptionParser & _OP);
		std::string getSingleString(std::string name, int index);
		std::vector<std::string> getArrayOfStrings(std::string name, int index);
		int getOptionNumberOfMultiples(std::string name);
		std::string translateEnv(std::string _input);
		void addFreeArgumentsToDefault();

		std::vector<std::string> errorMessages;
		std::vector<std::string> warningMessages;


		//std::string makeRunConfFile() const;

		/***********************************************
		 *  Internal storage of the options, 3 std::string vectors:
		 *
		 *  opts      the options
		 *  vals      the values of the options
		 *  source    source of the value is "a" for argument or "f" for file
		 ***********************************************/
		std::vector<std::string> opts;
		std::vector<std::vector<std::string> > vals;
		std::vector<std::vector<std::string> > source;
		std::vector<std::string> freeArgs;

		std::vector<std::vector<std::string> > shortOptEquivalent;

		std::string commandName;

		std::string confFile;

		std::vector<std::string> required;
		std::vector<std::string> allowed;
		std::vector<std::vector<std::string> > mutualExclusive;
		std::vector<std::vector<std::string> > dependentOn;
		std::vector<std::vector<std::string> > interDependent;
		std::vector<std::vector<std::string> > linked;
		std::vector<std::vector<std::string> > oneRequired;

		std::vector<std::string> missing;
		std::vector<std::string> disallowed;
		std::vector<std::string> ambiguous;
		std::vector<std::vector<std::string> > disallowedTogether;
		std::vector<std::vector<std::string> > missingDependency;
		std::vector<std::vector<std::string> > missingInterdipendenty;
		std::vector<std::vector<std::string> > missingLinked;
		std::vector<std::vector<std::string> > missingOneRequired;

		std::vector<std::string> defaultArguments;
		
		bool errorFlag;
		bool translateEnv_flag;	
};
inline void OptionParser::addError(std::string _error) {
	errorMessages.push_back(_error);
}
inline void OptionParser::addWarning(std::string _warning) {
	warningMessages.push_back(_warning);
}
}

#endif
