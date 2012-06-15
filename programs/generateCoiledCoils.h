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
#include <string>
#include <vector>

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Type of symmetry to apply
		required.push_back("symmetry");

		// super-helical Radius parameters (LOW, HIGH, STEP)
		required.push_back("superHelicalRadius");

		// alpha-helical phase angle parameters (LOW, HIGH, STEP)
		required.push_back("alphaHelicalPhaseAngle");

		// super-helical pitch angle (LOW, HIGH, STEP)
		optional.push_back("superHelicalPitchAngle");

		// Number of residues
		required.push_back("numberOfResidues");

		// Optional arguments
		optional.push_back("d2zTranslation");
		optional.push_back("superHelicalPhaseAngle");
		optional.push_back("name");


	}

	// Storage for the vales of each option
	std::string symmetry;
	std::vector<double> superHelicalRadius;
	std::vector<double> alphaHelicalPhaseAngle;
	int numberOfResidues;

	std::vector<double> d2zTranslation;
	std::vector<double> superHelicalPhaseAngle;
        std::vector<double> superHelicalPitchAngle;
	std::string name;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
