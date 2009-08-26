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
		required.push_back("energyTable");

		/************************
		     MC options
		*************************/

		// Annealing Schedule
		optional.push_back("annealScheduleType");
		optional.push_back("annealScheduleStartTemp");
		optional.push_back("annealScheduleEndTemp");
		optional.push_back("numberOfAnnealCycles");

		// Number of MC cycles (= number of single point mutations [ configuations ] to try)
		optional.push_back("numberOfCycles");

		// Initialization algorithm
		optional.push_back("initializationAlgorithm");
		optional.push_back("initializationConfiguration");

		// Number of configurations to store as we sample the energy landscape
		optional.push_back("numberOfStoredConfigurations");
		
		// Random seed -1 means time based
		optional.push_back("randomSeed");

		// Dead-End Elimination
		optional.push_back("DEE");

		// Debug,help options
		optional.push_back("debug");
		optional.push_back("help");

		// Configuration file..
		defaultArgs.push_back("configfile");

	}




	// Storage for the vales of each option
	string configfile;
	string energyTable;
	string annealType;
	double annealStart;
	double annealEnd;
	int    numberOfAnnealCycles;
	int    numCycles;
	string initAlgorithm;
	string initConfiguration;
	int    numStoredConfigurations;
	int    randomSeed;
	bool   DEE;

	bool debug;
	bool help;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

void cleanExit(int sig);
