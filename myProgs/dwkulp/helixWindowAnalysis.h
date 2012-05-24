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
		required.push_back("rscript");

		optional.push_back("sel");	
		optional.push_back("window");
		optional.push_back("pymol");
		optional.push_back("debug");
		optional.push_back("windowAtA");
		optional.push_back("windowAtB");
	}

	// Storage for the vales of each optional
	string pdb;
	vector<string> sel;
	int window;
	bool pymol;
	bool debug;
	string rscript;
	int windowAtA;
	int windowAtB;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
vector<double> getGeometricParameters(Line &_axis1, Line &_axis2, Line &_axisCenter, CartesianPoint &_phasePoint1, CartesianPoint &_phasePoint2);
