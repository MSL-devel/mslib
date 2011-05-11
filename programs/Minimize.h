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
#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Required
		required.push_back("pdb");

		// Optional
		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("selection");
		optional.push_back("springconstant");

		optional.push_back("method");
		optional.push_back("steps");
		optional.push_back("stepsize");
		optional.push_back("dielectric");
		optional.push_back("distanceDielectric");
		optional.push_back("cuton");
		optional.push_back("cutoff");
		optional.push_back("cutnb");

		// Configuration file..
		defaultArgs.push_back("configfile");

	}




	// Storage for the vales of each option
	string configFile;
	string pdb;
	string topfile;
	string parfile;
	string method;
	string outpdb;
	string selection;
	double springConstant;
	int steps;
	double stepSize;
	double dielectric;
	double ctonnb;
	double ctofnb;
	double cutnb;
	bool distanceDependentElectrostatics;


	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);


