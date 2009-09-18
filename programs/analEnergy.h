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
struct Options {

	// Set up options here...
	Options(){

		// Energy Table
		required.push_back("pdb");
		required.push_back("positions");

		/************************
		     MC options
		*************************/

		// Debug,help options
		optional.push_back("parfile");
		optional.push_back("topfile");

		optional.push_back("debug");
		optional.push_back("help");

		// Configuration file..
		defaultArgs.push_back("configfile");

	}




	// Storage for the vales of each option
	string pdb;
	string topfile;
	string parfile;
	string configfile;
	vector<string> positions;

	bool debug;
	bool help;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

