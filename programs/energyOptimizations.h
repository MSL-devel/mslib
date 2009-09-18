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

#include "System.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "PairwiseEnergyCalculator.h"
#include "OptionParser.h"

struct StructureOptions {

	// Set up options here...
	StructureOptions(){


		required.push_back("pdb");
		required.push_back("rotlib");
		required.push_back("position");
		required.push_back("identities");
		required.push_back("rotNumbers");

		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("linkedPosition");
		optional.push_back("flexNeighborDistance");
		optional.push_back("flexNeighborRotamers");

		// Debug,help options
		optional.push_back("help");

		// Configuration file..
		defaultArgs.push_back("config");

	};

	// Storage for the values of each option
	string pdb;
	string rotlib;
	string topfile;
	string parfile;

	vector<string> positions;
	vector<vector<string> > identities;
	vector<vector<int> > rotNumbers;
	vector<vector<string> > linkedPositions;
	double flexNeighborDistance;
	int flexNeighborRotamers;
	string configFile;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;

};



struct EnergyTableOptions {


	EnergyTableOptions(){
		/************************
		     Optionals
		*************************/
		optional.push_back("energyTableName");
		optional.push_back("partialEnergyTable");
		optional.push_back("dielectric");
		optional.push_back("distanceDielectric");
		optional.push_back("vdwScale");


		// Configuration files
		optional.push_back("structureConfig");
		defaultArgs.push_back("config");

	}




	// Storage for the values of each option
	StructureOptions structOpt;
	string structureConfig;
	string configFile;
	string energyTableName;
	string partialEnergyTable;
	double dielectric;
	bool distanceDependentElectrostatics;
	double vdwScale;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


struct MonteCarloOptions {

	// Set up options here...
	MonteCarloOptions(){

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

		// Configuration files
		optional.push_back("structureConfig");
		defaultArgs.push_back("configfile");

	}




	// Storage for the values of each option
	StructureOptions structOpt;
	string structureConfig;

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


struct LinearProgrammingOptions {

	// Set up options here...
	LinearProgrammingOptions(){

		// Energy Table
		required.push_back("energyTable");

		/************************
		     LP options
		*************************/

		// Configuration file..
		optional.push_back("structureConfig");
		defaultArgs.push_back("configfile");


	}

	// Storage for the values of each option
	StructureOptions structOpt;
	string structureConfig;

	string configfile;
	string energyTable;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};



EnergyTableOptions setupEnergyTableOptions(int theArgc, char * theArgv[]);
MonteCarloOptions setupMonteCarloOptions(int theArgc, char * theArgv[]);
LinearProgrammingOptions setupLinearProgrammingOptions(int theArgc, char * theArgv[]);

StructureOptions setupStructureOptions(string _confFile);


void createSystem(StructureOptions &_structOpt, System &_sys);

void cleanExit(int sig);


StructureOptions setupStructureOptions(string _confFile){

	// Create the options
	StructureOptions opt;
	

	// Parse the options
	OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	opt.configFile = _confFile;
	
	OP.readFile(opt.configFile);
	if (OP.fail()) {
		string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
		cerr << "ERROR 1111 "<<errorMessages<<endl;
	}


	if (OP.countOptions() == 0){
		cout << "Options:" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << "rotlib ROTLIB\n";
		cout << "topfile TOPFILE\n";
		cout << "parfile PARFILE\n";
		cout << "\n#For each variable position..\n";
		cout << "position     CHAIN_RESIDUE\n";
		cout << "identities   RES\n";
		cout << "rotNumbers   NUM\n";
		cout << "# linkedPosition CHAIN_RESIDUE_ORIG CHAIN_RESIDUE_LINK1 CHAIN_RESIDUE_LINK2 \n";
		cout << "# flexNeighborDistance DISTANCE\n";
		cout << "# flexNeighborRotamers NUMBER\n";
		cout << endl;
		exit(0);
	}




	// Input structure file
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}


	// Input rotamer library
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "ERROR 1111 no rotlib specified."<<endl;	
		exit(1111);
	}



	// Positions
	int index = 0;
	while (true) {
		string tmp = OP.getString("position",index);
		if (OP.fail()) {
			break;
		}

		opt.positions.push_back(tmp);
		index++;
	}
	if (opt.positions.size() == 0) {
		cerr << "The number of entries for option 'position' needs to be at least \n";
	}
	
	

	// Linked Positions
	index = 0;
	while (true) {
		vector<string> tmp = OP.getStringVector("linkedPosition", index);
		if (OP.fail()) {
			break;
		}

		opt.linkedPositions.push_back(tmp);
		index++;
	}

	// Identities
	index = 0;
	while (true) {
		vector<string> tmp = OP.getStringVector("identities", index);
		if (OP.fail()) {
			break;
		}
		opt.identities.push_back(tmp);
		index++;
	}
	if (opt.positions.size() != opt.identities.size()) {
		cerr << "The number of entries for options positions and identities need to be identical. Since its not, the program will take WT Identities\n";
	}

	// RotNumbers
	index = 0;
	while (true) {
		vector<int> tmp = OP.getIntVector("rotNumbers", index);
		if (OP.fail()) {
			break;
		}
		opt.rotNumbers.push_back(tmp);
		index++;
	}
	if (opt.positions.size() != opt.rotNumbers.size()) {
		cerr << "The number of entries for options positions and rotNumbers need to be identical\n";
	} else {
		for (unsigned i=0; i<opt.identities.size(); i++) {
			if (opt.rotNumbers[i].size() != opt.identities[i].size()) {
				cerr<< "The series of number of rotamers \"rotNumbers\" needs to be identical to the series of identities at position " + opt.positions[i] + "\n";
			}
		}
	}


	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		cerr << "WARNING no topfile specified, using default /library/charmmTopPar/top_all22_prot.inp\n";
		opt.topfile = "/library/charmmTopPar/top_all22_prot.inp";
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		cerr << "WARNING no parfile specified, using default /library/charmmTopPar/par_all22_prot.inp\n";
		opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
	}

	opt.flexNeighborDistance = OP.getDouble("flexNeighborDistance");
	if (OP.fail()){
		opt.flexNeighborDistance = -1;
	} else {
		cout << "Using flexible neighbor flag, defining neighbors by distance: "<<opt.flexNeighborDistance<<endl;
	}

	opt.flexNeighborRotamers = OP.getInt("flexNeighborRotamers");
	if (OP.fail()){
		if (opt.flexNeighborRotamers != -1){
			opt.flexNeighborRotamers = 10;
		} else {
			opt.flexNeighborRotamers = -1;
		}



	}else {
		cout << "Using flexible neighbor flag, number of neighbor rotamers "<<opt.flexNeighborRotamers<<endl;
	}

	return opt;
}

EnergyTableOptions setupEnergyTableOptions(int theArgc, char * theArgv[]){

	// Create the options
	EnergyTableOptions opt;
	

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: energyTable conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "structureConfig STRUCTURE_CONF_FILE\n";
		cout << "#energyTableName FILE\n";
		cout << "#dielectric 80\n";
		cout << "#distanceDielectric\n";
		cout << "#vdwScale 1.0\n";
		cout << endl;
		exit(0);
	}

	opt.configFile = OP.getString("config");
	
	if (opt.configFile != "") {
		OP.readFile(opt.configFile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (OP.fail()){
		cerr << "ERROR 111 structureConfig file required\n";
		exit(1111);
	} 
	opt.structOpt = setupStructureOptions(opt.structureConfig);
	

	opt.energyTableName = OP.getString("energyTableName");
	if (OP.fail()){
		opt.energyTableName = "energyTable.txt";
	}


	opt.dielectric = OP.getDouble("dielectric");
	if (OP.fail()){
		cerr << "WARNING no dielectric specified, using 80 as default.\n";
		opt.dielectric = 80;

	}
	opt.distanceDependentElectrostatics = OP.getBool("distanceDielectric");
	if (OP.fail()){
		cerr << "WARNING no distanceDielectric specified, using 'false' as default.\n";
		opt.distanceDependentElectrostatics = false;
	}

	opt.vdwScale = OP.getDouble("vdwScale");
	if (OP.fail()){
		opt.vdwScale = 1.0;
	}
	return opt;
}

MonteCarloOptions setupMonteCarloOptions(int theArgc, char * theArgv[]){
	// Create the options
	MonteCarloOptions opt;
    

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);    
	OP.setAllowed(opt.optional);
	//    OP.setShortOptionEquivalent(opt.equivalent);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"

	opt.configfile = OP.getString("configfile");
    
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}


	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "optimizeMC CONF\n";
		exit(0);
	}



	if (OP.getBool("help")){
	    
		cout << "# Options for optimizeMC\n\n";
		cout << "# Energy Table\n";
		cout << "energyTable energy.txt\n\n";
		cout << "# Annealing Schedule\n";
		cout << "#     can be LINEAR, EXPONENTIAL, LINEAR_CYCLES, EXPONENTIAL_CYCLES\n";
		cout << "annealScheduleType EXPONENTIAL\n";
		cout << "annealScheduleStartTemp 1000\n";
		cout << "annealScheduleEndTemp   1 \n\n";
		cout << "numberOfAnnealCycles    1 \n\n";
		cout << "# Number of MC Cycles \n";
		cout << "numberOfCycles 10000\n\n";
		cout << "# Initial Rotamer Configuration set via initializationAlgorithm\n";
		cout << "#     can be RANDOM LOWESTSELF QUICKSCAN USERINPUT\n";
		cout << "initializationAlgorithm LOWESTSELF\n\n";
		cout << "# Initial Configuration set by user, RotamerNumberForPosition1:RotamerNumberForPosition2:\n";
		cout << "#   Notice the last character is a ':'\n";
		cout << "# initializationConfiguration 0:0:0:0:\n\n";
		cout << "# Number of Rotamer Configurations to report , top 100 lowest energy scoring found\n";
		cout << "numberOfStoredConfigurations 100\n\n";
		cout << "# Random Seed to use, -1 means create a time-based seed\n";
		cout << "randomSeed 838201\n\n";
		cout << "# Structure Configuration File to Output PDBs\n";
		cout << "structureConfig FILENAME\n";
		exit(0);
        
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (!OP.fail()){
		opt.structOpt = setupStructureOptions(opt.structureConfig);   
	}

	opt.debug = OP.getBool("debug");

	opt.energyTable = OP.getString("energyTable");
	if (OP.fail()){
		cerr << "ERROR 1111 energyTable not specified."<<endl;	
		exit(1111);
	}


	opt.annealType = OP.getString("annealScheduleType");
	if (OP.fail()){
		opt.annealType = "LINEAR";
	}

	opt.annealStart = OP.getDouble("annealScheduleStartTemp");
	if (OP.fail()){
		opt.annealStart = 1000;
	}

	opt.annealEnd = OP.getDouble("annealScheduleEndTemp");
	if (OP.fail()){
		opt.annealEnd = 1;
	}

	opt.numberOfAnnealCycles = OP.getInt("numberOfAnnealCycles");
	if (OP.fail()){
		opt.numberOfAnnealCycles = 1;
	}

	opt.numCycles = OP.getInt("numberOfCycles");
	if (OP.fail()){
		opt.numCycles = 10000;
	}

	opt.initAlgorithm = OP.getString("initializationAlgorithm");
	if (OP.fail()){
		opt.initAlgorithm = "LOWESTSELF";
	}
	opt.initConfiguration = OP.getString("initializationConfiguration");
	if (OP.fail()){
		opt.initConfiguration = "";
	}

	opt.numStoredConfigurations = OP.getInt("numberOfStoredConfigurations");
	if (OP.fail()){
		opt.numStoredConfigurations = 50;
	}
    
	opt.randomSeed = OP.getInt("randomSeed");
	if (OP.fail()){
		opt.randomSeed = -1;
	}

    
	opt.DEE  = OP.getBool("DEE");
	if (OP.fail()){
		opt.DEE = false;
	}


	cout << OP<<endl;
	return opt;
}

LinearProgrammingOptions setupLinearProgrammingOptions(int theArgc, char * theArgv[]) {
	// Create the options
	LinearProgrammingOptions opt;
    

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);    
	OP.setAllowed(opt.optional);
	//    OP.setShortOptionEquivalent(opt.equivalent);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"


	opt.configfile = OP.getString("configfile");
    
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}


	if (OP.getBool("help")){
	    
		cout << "# Options for optimizeLP\n\n";
		cout << "# Energy Table\n";
		cout << "energyTable energy.txt\n\n";
		cout << "#structureConfig STRUCT_CONFIG_FILE\n";

		exit(0);
        
	}

	opt.energyTable = OP.getString("energyTable");
	if (OP.fail()){
		cerr << "ERROR 1111 energyTable not specified."<<endl;	
		exit(1111);
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (!OP.fail()){
		opt.structOpt = setupStructureOptions(opt.structureConfig);
	}


	cout << OP<<endl;
	return opt;
}

void createSystem(StructureOptions &_opt, System &_sys) {


	// Read-in PDB as initialSystem
	cout << "Read in PDB"<<endl;
	System initialSystem;
	initialSystem.readPdb(_opt.pdb);

	// Build Map of Variable Positions
	cout << "Build Designable PolymerSequence"<<endl;
	map<string,map<int,int> > variablePositionMap;
	uint psize = _opt.positions.size();
	for (uint v = 0; v < psize;v++){
		cout << "Position "<<v<<" '"<<_opt.positions[v]<<"'"<<endl;
		vector<string> tmp = MslTools::tokenizeAndTrim(_opt.positions[v],"_");
		if (tmp.size() < 2){
			cerr << "ERROR 2222: Position "<<v<<" '"<<_opt.positions[v]<<"' does not have CHAIN_RESIDUE format\n";
			exit(2222);
		}

		if (initialSystem.exists(tmp[0],tmp[1])){
			cout << "\t"<<tmp[0]<<"."<<tmp[1]<<"."<<endl;
			variablePositionMap[tmp[0]][MslTools::toInt(tmp[1])] = v;					


			// If add neighbors automatically
			if (_opt.flexNeighborDistance != -1){
				Residue &res = initialSystem.getResidue(initialSystem.getPositionIndex(tmp[0],tmp[1]));

				cout << "Finding Neighbors"<<endl;
				vector<int> neighbors = res.findNeighbors(_opt.flexNeighborDistance);
				cout << "Neighbors are: "<<neighbors.size()<<endl;

				
				for (uint n = 0; n < neighbors.size();n++){
					Residue &neighbor = initialSystem.getResidue(neighbors[n]);
					cout << "\t"<<neighbors[n]<<" "<<neighbor.getResidueNumber()<<endl;
				}
				for (uint n = 0; n < neighbors.size();n++){

					Residue &neighbor = initialSystem.getResidue(neighbors[n]);


					if (neighbor.getResidueName() == "ALA" || neighbor.getResidueName() == "PRO" || neighbor.getResidueName() ==  "GLY"){
						continue;
					}

					// Skip if in map already...
					if (variablePositionMap[neighbor.getChainId()][neighbor.getResidueNumber()]){
						continue;
					}

					// Add to _opt.positions, _opt.rotNumbers, _opt.identities and variablePositionMap
					char tmpStr[80];
					sprintf(tmpStr, "%1s_%d",neighbor.getChainId().c_str(),neighbor.getResidueNumber());
					_opt.positions.push_back(tmpStr);
					_opt.identities.push_back(vector<string>(1,neighbor.getResidueName()));

					if (neighbor.getResidueName() == "SER" || neighbor.getResidueName() == "THR" || neighbor.getResidueName() == "CYS" || neighbor.getResidueName() == "VAL"){
						_opt.rotNumbers.push_back(vector<int>(1,10));
					} else {
						_opt.rotNumbers.push_back(vector<int>(1,_opt.flexNeighborRotamers));
					}

					
					variablePositionMap[neighbor.getChainId()][neighbor.getResidueNumber()] = _opt.positions.size() - 1;
					
				}
			}

		}
	}


	// Add linked positions, first position in each linkedPositions[INDEX] is the "master" position
	initialSystem.setLinkedPositions(_opt.linkedPositions);


	// Build PolymerSequence.
	cout << "build it"<<endl;
	PolymerSequence pseq(initialSystem, variablePositionMap, _opt.identities);
	cout << "PolymerSequence: "<<endl<<pseq.toString()<<endl;

	// Build a new system from polymer sequence and create energySet from energy terms
	cout << "Build Charmm System"<<endl;
	string topFile = _opt.topfile;
	string parFile = _opt.parfile;
	cout << "Use toppar " << topFile << ", " << parFile << endl;
	CharmmSystemBuilder CSB(topFile,parFile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(_sys,pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	int numAssignedAtoms = _sys.assignCoordinates(initialSystem.getAtoms(),false);
	fprintf(stdout, "Assigned %8d atoms\n",numAssignedAtoms);

	// Build the all atoms without coordinates (not in initial PDB)
	cout <<"Build Atoms"<<endl;
	_sys.buildAllAtoms(); 

	// Set linked positions in our new sequence-built system
	_sys.setLinkedPositions(_opt.linkedPositions);


	// Write initially built file.
	string filename = "/tmp/initialBuild.pdb";
	cout << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(_sys.getAtoms())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();

	
	// Build rotamers
	cout << "Read rotamer library " << _opt.rotlib << " and load rotamers"<<endl;	
	SystemRotamerLoader sysRot(_sys, _opt.rotlib);
	map<string,map<int,int> >::iterator it1;
	map<int,int>::iterator it2;
	for (uint i = 0; i < _sys.positionSize();i++){
		Position &pos = _sys.getPosition(i);

		string chainId = pos.getChainId();
		int resNum     = pos.getResidueNumber();

		// If a "slave" position, then take its "master" positions chainId and residue number.
		if (pos.getLinkedPositionType() == Position::SLAVE){
			chainId = pos.getLinkedPositions()[0]->getChainId();
			resNum  = pos.getLinkedPositions()[0]->getResidueNumber();
			cout << " FOUND LINKED POS: "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" to "<<chainId<<" " <<resNum<<endl;
		}	
		
		int index = -1;
		it1 = variablePositionMap.find(chainId);
		if (it1 != variablePositionMap.end()){
			
			it2 = it1->second.find(resNum);
			if (it2 != it1->second.end()){
				index = it2->second;
			}
		}
		
		if (index != -1){
			for (uint j = 0; j < _opt.identities[index].size();j++){
				cout << "Loading "<<_opt.rotNumbers[index][j]<<" rotamers of type "<<_opt.identities[index][j]<<" at postion "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName();
				if (pos.getLinkedPositionType() == Position::SLAVE){
					cout << " LINKED TO "<< chainId<<" "<<resNum;
				}
				cout <<endl;

				// Test for position specific library...
				char posSpecificLib[80];
				sprintf(posSpecificLib,"POS_%1s_%04d",pos.getChainId().c_str(),pos.getResidueNumber());


				// Load rotamers
				if (sysRot.getRotamerLibrary()->libraryExists((string)posSpecificLib)) {
					int lastRotamerIndex = _opt.rotNumbers[index][j]-1;
					if (sysRot.getRotamerLibrary()->size((string)posSpecificLib,_opt.identities[index][j]) < _opt.rotNumbers[index][j]){
						cerr << "WARNING "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName()<<" has a specific rotamer library but not enough rotamers it has: "<<sysRot.getRotamerLibrary()->size((string)posSpecificLib,_opt.identities[index][j])<<endl;
						lastRotamerIndex =sysRot.getRotamerLibrary()->size((string)posSpecificLib,_opt.identities[index][j])-1;
					}
					
					sysRot.loadRotamers(&pos, posSpecificLib, _opt.identities[index][j], 0, lastRotamerIndex); 
				} else {
					sysRot.loadRotamers(&pos, "BALANCED-200", _opt.identities[index][j], 0, _opt.rotNumbers[index][j]-1); 
				}
				
				
			}
		}
	}
}


void changeRotamerState(StructureOptions &_opt, System &_sys, vector<int> &_rotamerState){


		for (uint i = 0; i < _opt.positions.size();i++){
			//cout << "Working on absres: "<<_opt.positions[i]<<endl;

			vector<string> tmp = MslTools::tokenizeAndTrim(_opt.positions[i],"_");
			if (tmp.size() < 2){
				cerr << "ERROR 2222: Position "<<i<<" '"<<_opt.positions[i]<<"' does not have CHAIN_RESIDUE format\n";
				exit(2222);
			}

			if (!_sys.exists(tmp[0],tmp[1])){
				cerr << "ERROR 2222: Chain '"<<tmp[0]<<"' and residue '"<<tmp[1]<<"' doesn't exist in pdb file: '"<<_opt.pdb<<"'"<<endl;
				exit(2222);
			}
			Position &pos = _sys.getPosition(tmp[0],tmp[1]);

			pos.setActiveRotamer(_rotamerState[i]);
		}
}
