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
struct Options {

	// Set up options here...
	Options(){


		required.push_back("pdblist");
		required.push_back("regex");

		/************************
		     Optionals
		*************************/
		optional.push_back("outdir");
		optional.push_back("config");

		// Debug,help options
		optional.push_back("debug");
		optional.push_back("help");


		// Configuration file..
		defaultArgs.push_back("config");

	}




	// Storage for the vales of each option
	std::string pdbs;
	std::string regex;
	std::string outdir;
	std::string configFile;
	bool debug;
	bool help;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
