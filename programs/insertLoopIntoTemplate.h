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
#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("template");
		optional.push_back("templateChain");
		required.push_back("fragment");
		required.push_back("fragmentChain");

		optional.push_back("templateStem1");
		optional.push_back("templateStem2");
		optional.push_back("clashCheck");
		optional.push_back("numClashes");
		optional.push_back("includeTemplateStems");
		optional.push_back("checkCaCaDistances");
		optional.push_back("outputRosettaFiles");
	}




	// Storage for the vales of each option
	string templatePDB;
	string templateChain;
        string templateStem1;
        string templateStem2;
        string fragmentPDB;
        bool clashCheck;
        bool includeTemplateStems;
        bool checkCaCaDistances;
	int numClashes;
        string fragmentChain;
        string outputRosettaFiles;
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
