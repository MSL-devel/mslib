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

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("pdb");
		required.push_back("doffile");
		required.push_back("selection");
		optional.push_back("phiPsiTable");
		optional.push_back("computeDeltaNormalizedSASA");
		optional.push_back("charmmPar");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	std::string pdb;
	std::string dofFile;
	std::string selection;
	std::string phiPsiTable;
        bool computeDeltaNormalizedSASA;
        std::string charmmParameterFile;
	bool debug;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

