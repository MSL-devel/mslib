/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

using namespace std;
using namespace MSL;

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
	std::string pdb;
	std::string rotlib;
	std::string topfile;
	std::string parfile;

	std::vector<std::string> positions;
	std::vector<std::vector<std::string> > identities;
	std::vector<std::vector<int> > rotNumbers;
	std::vector<std::vector<std::string> > linkedPositions;
	double flexNeighborDistance;
	int flexNeighborRotamers;
	std::string configFile;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;

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
		required.push_back("structureConfig");
		defaultArgs.push_back("config");

	}




	// Storage for the values of each option
	StructureOptions structOpt;
	std::string structureConfig;
	std::string configFile;
	std::string energyTableName;
	std::string partialEnergyTable;
	double dielectric;
	bool distanceDependentElectrostatics;
	double vdwScale;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;
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
		defaultArgs.push_back("config");

	}




	// Storage for the values of each option
	StructureOptions structOpt;
	std::string structureConfig;

	std::string configfile;
	std::string energyTable;
	std::string annealType;
	double annealStart;
	double annealEnd;
	int    numberOfAnnealCycles;
	int    numCycles;
	std::string initAlgorithm;
	std::string initConfiguration;
	int    numStoredConfigurations;
	int    randomSeed;
	bool   DEE;

	bool debug;
	bool help;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;
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
		defaultArgs.push_back("config");


	}

	// Storage for the values of each option
	StructureOptions structOpt;
	std::string structureConfig;

	std::string configfile;
	std::string energyTable;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;
};



EnergyTableOptions setupEnergyTableOptions(int theArgc, char * theArgv[]);
MonteCarloOptions setupMonteCarloOptions(int theArgc, char * theArgv[]);
LinearProgrammingOptions setupLinearProgrammingOptions(int theArgc, char * theArgv[]);

StructureOptions setupStructureOptions(std::string _confFile);


void createSystem(StructureOptions &_structOpt, MSL::System &_sys);

void cleanExit(int sig);


StructureOptions setupStructureOptions(std::string _confFile){

	// Create the options
	StructureOptions opt;
	

	// Parse the options
	MSL::OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	opt.configFile = _confFile;
	
	OP.readFile(opt.configFile);
	if (OP.fail()) {
		std::string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
		std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
	}


	if (OP.countOptions() == 0){
		std::cout << "Options:" << std::endl;
		std::cout << std::endl;
		std::cout << "\n";
		std::cout << "pdb PDB\n";
		std::cout << "rotlib ROTLIB\n";
		std::cout << "topfile TOPFILE\n";
		std::cout << "parfile PARFILE\n";
		std::cout << "\n#For each variable position..\n";
		std::cout << "position     CHAIN_RESIDUE\n";
		std::cout << "identities   RES\n";
		std::cout << "rotNumbers   NUM\n";
		std::cout << "# linkedPosition CHAIN_RESIDUE_ORIG CHAIN_RESIDUE_LINK1 CHAIN_RESIDUE_LINK2 \n";
		std::cout << "# flexNeighborDistance DISTANCE\n";
		std::cout << "# flexNeighborRotamers NUMBER\n";
		std::cout << std::endl;
		exit(0);
	}




	// Input structure file
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		std::cerr << "ERROR 1111 no pdb specified."<<std::endl;
		exit(1111);
	}


	// Input rotamer library
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		std::cerr << "ERROR 1111 no rotlib specified."<<std::endl;	
		exit(1111);
	}



	// Positions
	int index = 0;
	while (true) {
		std::string tmp = OP.getString("position",index);
		if (OP.fail()) {
			break;
		}

		opt.positions.push_back(tmp);
		index++;
	}
	if (opt.positions.size() == 0) {
		std::cerr << "The number of entries for option 'position' needs to be at least \n";
	}
	
	

	// Linked Positions
	index = 0;
	while (true) {
		std::vector<std::string> tmp = OP.getStringVector("linkedPosition", index);
		if (OP.fail()) {
			break;
		}

		opt.linkedPositions.push_back(tmp);
		index++;
	}

	// Identities
	index = 0;
	while (true) {
		std::vector<std::string> tmp = OP.getStringVector("identities", index);
		if (OP.fail()) {
			break;
		}
		opt.identities.push_back(tmp);
		index++;
	}
	if (opt.positions.size() != opt.identities.size()) {
		std::cerr << "The number of entries for options positions and identities need to be identical. Since its not, the program will take WT Identities\n";
	}

	// RotNumbers
	index = 0;
	while (true) {
		std::vector<int> tmp = OP.getIntVector("rotNumbers", index);
		if (OP.fail()) {
			break;
		}
		opt.rotNumbers.push_back(tmp);
		index++;
	}
	if (opt.positions.size() != opt.rotNumbers.size()) {
		std::cerr << "The number of entries for options positions and rotNumbers need to be identical\n";
	} else {
		for (unsigned i=0; i<opt.identities.size(); i++) {
			if (opt.rotNumbers[i].size() != opt.identities[i].size()) {
				std::cerr<< "The series of number of rotamers \"rotNumbers\" needs to be identical to the series of identities at position " + opt.positions[i] + "\n";
			}
		}
	}


	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		std::cerr << "WARNING no topfile specified, using default /library/charmmTopPar/top_all22_prot.inp\n";
		opt.topfile = "/library/charmmTopPar/top_all22_prot.inp";
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		std::cerr << "WARNING no parfile specified, using default /library/charmmTopPar/par_all22_prot.inp\n";
		opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
	}

	opt.flexNeighborDistance = OP.getDouble("flexNeighborDistance");
	if (OP.fail()){
		opt.flexNeighborDistance = -1;
	} else {
		std::cout << "Using flexible neighbor flag, defining neighbors by distance: "<<opt.flexNeighborDistance<<std::endl;
	}

	opt.flexNeighborRotamers = OP.getInt("flexNeighborRotamers");
	if (OP.fail()){
		if (opt.flexNeighborRotamers != -1){
			opt.flexNeighborRotamers = 10;
		} else {
			opt.flexNeighborRotamers = -1;
		}



	}else {
		std::cout << "Using flexible neighbor flag, number of neighbor rotamers "<<opt.flexNeighborRotamers<<std::endl;
	}

	return opt;
}

EnergyTableOptions setupEnergyTableOptions(int theArgc, char * theArgv[]){

	// Create the options
	EnergyTableOptions opt;
	

	// Parse the options
	MSL::OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	OP.readArgv(theArgc, theArgv);


	if (OP.countOptions() == 0){
		std::cout << "Usage: energyTable conf" << std::endl;
		std::cout << std::endl;
		std::cout << "\n";
		std::cout << "structureConfig STRUCTURE_CONF_FILE\n";
		std::cout << "#energyTableName FILE\n";
		std::cout << "#dielectric 80\n";
		std::cout << "#distanceDielectric\n";
		std::cout << "#vdwScale 1.0\n";
		std::cout << std::endl;
		exit(0);
	}

	opt.configFile = OP.getString("config");
	
	if (opt.configFile != "") {
		OP.readFile(opt.configFile);
		if (OP.fail()) {
			std::string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
			std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
		}
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (OP.fail()){
		std::cerr << "ERROR 111 structureConfig file required\n";
		exit(1111);
	} 
	opt.structOpt = setupStructureOptions(opt.structureConfig);
	

	opt.energyTableName = OP.getString("energyTableName");
	if (OP.fail()){
		opt.energyTableName = "energyTable.txt";
	}


	opt.dielectric = OP.getDouble("dielectric");
	if (OP.fail()){
		std::cerr << "WARNING no dielectric specified, using 80 as default.\n";
		opt.dielectric = 80;

	}
	opt.distanceDependentElectrostatics = OP.getBool("distanceDielectric");
	if (OP.fail()){
		std::cerr << "WARNING no distanceDielectric specified, using 'false' as default.\n";
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
	MSL::OptionParser OP;

	OP.setRequired(opt.required);    
	OP.setAllowed(opt.optional);
	//    OP.setShortOptionEquivalent(opt.equivalent);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"

	OP.readArgv(theArgc, theArgv);
	opt.configfile = OP.getString("configfile");
    
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			std::string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
			std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
		}
	}


	if (OP.countOptions() == 0){
		std::cout << "Usage:" << std::endl;
		std::cout << std::endl;
		std::cout << "optimizeMC CONF\n";
		exit(0);
	}



	if (OP.getBool("help")){
	    
		std::cout << "# Options for optimizeMC\n\n";
		std::cout << "# Energy Table\n";
		std::cout << "energyTable energy.txt\n\n";
		std::cout << "# Annealing Schedule\n";
		std::cout << "#     can be LINEAR, EXPONENTIAL, LINEAR_CYCLES, EXPONENTIAL_CYCLES\n";
		std::cout << "annealScheduleType EXPONENTIAL\n";
		std::cout << "annealScheduleStartTemp 1000\n";
		std::cout << "annealScheduleEndTemp   1 \n\n";
		std::cout << "numberOfAnnealCycles    1 \n\n";
		std::cout << "# Number of MC Cycles \n";
		std::cout << "numberOfCycles 10000\n\n";
		std::cout << "# Initial Rotamer Configuration set via initializationAlgorithm\n";
		std::cout << "#     can be RANDOM LOWESTSELF QUICKSCAN USERINPUT\n";
		std::cout << "initializationAlgorithm LOWESTSELF\n\n";
		std::cout << "# Initial Configuration set by user, RotamerNumberForPosition1:RotamerNumberForPosition2:\n";
		std::cout << "#   Notice the last character is a ':'\n";
		std::cout << "# initializationConfiguration 0:0:0:0:\n\n";
		std::cout << "# Number of Rotamer Configurations to report , top 100 lowest energy scoring found\n";
		std::cout << "numberOfStoredConfigurations 100\n\n";
		std::cout << "# Random Seed to use, -1 means create a time-based seed\n";
		std::cout << "randomSeed 838201\n\n";
		std::cout << "# Structure Configuration File to Output PDBs\n";
		std::cout << "structureConfig FILENAME\n";
		exit(0);
        
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (!OP.fail()){
		opt.structOpt = setupStructureOptions(opt.structureConfig);   
	}

	opt.debug = OP.getBool("debug");

	opt.energyTable = OP.getString("energyTable");
	if (OP.fail()){
		std::cerr << "ERROR 1111 energyTable not specified."<<std::endl;	
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


	std::cout << OP<<std::endl;
	return opt;
}

LinearProgrammingOptions setupLinearProgrammingOptions(int theArgc, char * theArgv[]) {
	// Create the options
	LinearProgrammingOptions opt;
    

	// Parse the options
	MSL::OptionParser OP;

	OP.setRequired(opt.required);    
	OP.setAllowed(opt.optional);
	//    OP.setShortOptionEquivalent(opt.equivalent);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"

	OP.readArgv(theArgc, theArgv);
	opt.configfile = OP.getString("configfile");
    
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			std::string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
			std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
		}
	}


	if (OP.getBool("help")){
	    
		std::cout << "# Options for optimizeLP\n\n";
		std::cout << "# Energy Table\n";
		std::cout << "energyTable energy.txt\n\n";
		std::cout << "#structureConfig STRUCT_CONFIG_FILE\n";

		exit(0);
        
	}

	opt.energyTable = OP.getString("energyTable");
	if (OP.fail()){
		std::cerr << "ERROR 1111 energyTable not specified."<<std::endl;	
		exit(1111);
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (!OP.fail()){
		opt.structOpt = setupStructureOptions(opt.structureConfig);
	}


	std::cout << OP<<std::endl;
	return opt;
}

void createSystem(StructureOptions &_opt, MSL::System &_sys) {


/*
			string chain = "";
			int resnum = 0;
			string icode = "";
			bool OK = MslTools::parsePositionId(_linkedPositions[v][t], chain, resnum, icode, 0);
			if (!OK) {
				cerr << "DEPRECATED USE OF UNDERSCORE SEPARATED IDENTIFIERS (I.E. \"A_37\"), USE COMMA SEPARATION (\"A,37\") in void System::setLinkedPositions(vector<vector<string> > &_linkedPositions)" << endl;
				vector<string> tmp = MslTools::tokenizeAndTrim(_linkedPositions[v][t],"_");
				string newId;
				if (tmp.size() > 1) {
					newId = tmp[0] + "," + tmp[1];
				}
				OK = MslTools::parsePositionId(newId, chain, resnum, icode, 0);

			}
*/

	// Read-in PDB as initialSystem
	std::cout << "Read in PDB"<<std::endl;
	MSL::System initialSystem;
	initialSystem.readPdb(_opt.pdb);

	// Build Map of Variable Positions
	std::cout << "Build Designable PolymerSequence"<<std::endl;
	std::map<std::string,std::map<int,int> > variablePositionMap;
	uint psize = _opt.positions.size();
	for (uint v = 0; v < psize;v++){
		std::cout << "Position "<<v<<" '"<<_opt.positions[v]<<"'"<< std::endl;

		std::string chain = "";
		int resnum = 0;
		std::string icode = "";
		bool OK = MslTools::parsePositionId(_opt.positions[v], chain, resnum, icode, 0);
		/*
		if (!OK) {
			std::cerr << "DEPRECATED USE OF UNDERSCORE SEPARATED IDENTIFIERS (I.E. \"A_37\"), USE COMMA SEPARATION (\"A,37\")" << std::endl;
			std::vector<std::string> tmp = MslTools::tokenizeAndTrim(_opt.positions[v],"_");
			std::string newId;
			if (tmp.size() > 1) {
				newId = tmp[0] + "," + tmp[1];
			}
			OK = MslTools::parsePositionId(newId, chain, resnum, icode, 0);

		}
		*/

	//	std::vector<std::string> tmp = MslTools::tokenizeAndTrim(_opt.positions[v],"_");
		//if (tmp.size() < 2){
		if (!OK){
			//std::cerr << "ERROR 2222: Position "<<v<<" '"<<_opt.positions[v]<<"' does not have CHAIN_RESIDUE format\n";
			std::cerr << "ERROR 2222: Position "<<v<<" '"<<_opt.positions[v]<<"' does not have CHAIN,RESIDUE format\n";
			exit(2222);
		}

		//if (initialSystem.exists(tmp[0],tmp[1])){
		if (initialSystem.positionExists(chain, resnum, icode)){
		//	std::out << "\t"<<tmp[0]<<"."<<tmp[1]<<"."<<std::endl;
			//variablePositionMap[tmp[0]][MslTools::toInt(tmp[1])] = v;					
			variablePositionMap[chain][resnum] = v;					


			// If add neighbors automatically
			if (_opt.flexNeighborDistance != -1){
			//	Residue &res = initialSystem.getResidue(initialSystem.getPositionIndex(tmp[0],tmp[1]));
				MSL::Residue &res = initialSystem.getResidue(chain, resnum, icode); // get the current identity

				std::cout << "Finding Neighbors"<<std::endl;
				std::vector<int> neighbors = res.findNeighbors(_opt.flexNeighborDistance);
				std::cout << "Neighbors are: "<<neighbors.size()<<std::endl;

				
				for (uint n = 0; n < neighbors.size();n++){
					MSL::Residue &neighbor = initialSystem.getResidue(neighbors[n]);
					std::cout << "\t"<<neighbors[n]<<" "<<neighbor.getResidueNumber()<<std::endl;
				}
				for (uint n = 0; n < neighbors.size();n++){

					MSL::Residue &neighbor = initialSystem.getResidue(neighbors[n]);


					if (neighbor.getResidueName() == "ALA" || neighbor.getResidueName() == "PRO" || neighbor.getResidueName() ==  "GLY"){
						continue;
					}

					// Skip if in std::map already...
					if (variablePositionMap[neighbor.getChainId()][neighbor.getResidueNumber()]){
						continue;
					}

					// Add to _opt.positions, _opt.rotNumbers, _opt.identities and variablePositionMap
				//	char tmpStr[80];
				//	sprintf(tmpStr, "%1s_%d",neighbor.getChainId().c_str(),neighbor.getResidueNumber());
				//	_opt.positions.push_back(tmpStr);
					_opt.positions.push_back(MslTools::getPositionId(neighbor.getChainId(), neighbor.getResidueNumber(), neighbor.getResidueIcode()));
					_opt.identities.push_back(std::vector<std::string>(1,neighbor.getResidueName()));

					if (neighbor.getResidueName() == "SER" || neighbor.getResidueName() == "THR" || neighbor.getResidueName() == "CYS" || neighbor.getResidueName() == "VAL"){
						_opt.rotNumbers.push_back(std::vector<int>(1,10));
					} else {
						_opt.rotNumbers.push_back(std::vector<int>(1,_opt.flexNeighborRotamers));
					}

					
					variablePositionMap[neighbor.getChainId()][neighbor.getResidueNumber()] = _opt.positions.size() - 1;
					
				}
			}

		}
	}


	// Add linked positions, first position in each linkedPositions[INDEX] is the "master" position
	initialSystem.setLinkedPositions(_opt.linkedPositions);


	// Build PolymerSequence.
	std::cout << "build it"<<std::endl;
	MSL::PolymerSequence pseq(initialSystem, variablePositionMap, _opt.identities);
	std::cout << "PolymerSequence: "<<std::endl<<pseq.toString()<<std::endl;

	// Build a new system from polymer sequence and create energySet from energy terms
	std::cout << "Build Charmm System"<<std::endl;
	std::string topFile = _opt.topfile;
	std::string parFile = _opt.parfile;
	std::cout << "Use toppar " << topFile << ", " << parFile << std::endl;
	MSL::CharmmSystemBuilder CSB(_sys,topFile,parFile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	int numAssignedAtoms = _sys.assignCoordinates(initialSystem.getAtomPointers(),false);
	fprintf(stdout, "Assigned %8d atoms\n",numAssignedAtoms);

	// Build the all atoms without coordinates (not in initial PDB)
	std::cout <<"Build Atoms"<<std::endl;
	_sys.buildAllAtoms(); 

	// Set linked positions in our new sequence-built system
	_sys.setLinkedPositions(_opt.linkedPositions);


	// Write initially built file.
	std::string filename = "/tmp/initialBuild.pdb";
	std::cout << "Write pdb " << filename << std::endl;
	MSL::PDBWriter writer;
	writer.open(filename);
	if (!writer.write(_sys.getAtomPointers())) {
		std::cerr << "Problem writing " << filename << std::endl;
	}
	writer.close();

	
	// Build rotamers
	std::cout << "Read rotamer library " << _opt.rotlib << " and load rotamers"<<std::endl;	
	MSL::SystemRotamerLoader sysRot(_sys, _opt.rotlib);
	std::map<std::string,std::map<int,int> >::iterator it1;
	std::map<int,int>::iterator it2;
	for (uint i = 0; i < _sys.positionSize();i++){
		MSL::Position &pos = _sys.getPosition(i);

		std::string chainId = pos.getChainId();
		int resNum     = pos.getResidueNumber();

		// If a "slave" position, then take its "master" positions chainId and residue number.
		if (pos.getLinkedPositionType() == MSL::Position::SLAVE){
			chainId = pos.getLinkedPositions()[0]->getChainId();
			resNum  = pos.getLinkedPositions()[0]->getResidueNumber();
			std::cout << " FOUND LINKED POS: "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" to "<<chainId<<" " <<resNum<<std::endl;
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
				std::cout << "Loading "<<_opt.rotNumbers[index][j]<<" rotamers of type "<<_opt.identities[index][j]<<" at postion "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName();
				if (pos.getLinkedPositionType() == MSL::Position::SLAVE){
					std::cout << " LINKED TO "<< chainId<<" "<<resNum;
				}
				std::cout <<std::endl;

				// Test for position specific library...
				char posSpecificLib[80];
				sprintf(posSpecificLib,"POS_%1s_%04d",pos.getChainId().c_str(),pos.getResidueNumber());


				// Load rotamers
				if (sysRot.getRotamerLibrary()->libraryExists((std::string)posSpecificLib)) {
					int lastRotamerIndex = _opt.rotNumbers[index][j]-1;
					if (sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_opt.identities[index][j]) < _opt.rotNumbers[index][j]){
						std::cerr << "WARNING "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName()<<" has a specific rotamer library but not enough rotamers it has: "<<sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_opt.identities[index][j])<<std::endl;
						lastRotamerIndex =sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_opt.identities[index][j])-1;
					}
					
					//sysRot.loadRotamers(&pos, posSpecificLib, _opt.identities[index][j], 0, lastRotamerIndex); 
					sysRot.loadRotamers(&pos, _opt.identities[index][j], 0, lastRotamerIndex, posSpecificLib); 
				} else {
					//sysRot.loadRotamers(&pos, "BALANCED-200", _opt.identities[index][j], 0, _opt.rotNumbers[index][j]-1); 
					sysRot.loadRotamers(&pos, _opt.identities[index][j], 0, _opt.rotNumbers[index][j]-1, ""); 
				}
				
				
			}
		}
	}



}


void changeRotamerState(StructureOptions &_opt, MSL::System &_sys, std::vector<int> &_rotamerState){


		for (uint i = 0; i < _opt.positions.size();i++){
			//std::cout << "Working on absres: "<<_opt.positions[i]<<std::endl;

		//	std::string chain = "";
		//	int resnum = 0;
		//	std::string icode = "";
		//	bool OK = MslTools::parsePositionId(_opt.positions[i], chain, resnum, icode, 0);
			/*
			std::vector<std::string> tmp = MslTools::tokenizeAndTrim(_opt.positions[i],"_");
			if (tmp.size() < 2){
				std::cerr << "ERROR 2222: Position "<<i<<" '"<<_opt.positions[i]<<"' does not have CHAIN_RESIDUE format\n";
				exit(2222);
			}
			*/
			//if (!_sys.positionExists(tmp[0],tmp[1])){
			if (!_sys.positionExists(_opt.positions[i])){
				std::cerr << "ERROR 2222: Position '"<<_opt.positions[i]<<"' doesn't exist in pdb file: '"<<_opt.pdb<<"'"<<std::endl;
				exit(2222);
			}
			//Position &pos = _sys.getPosition(tmp[0],tmp[1]);
			MSL::Position &pos = _sys.getPosition(_opt.positions[i]);

			pos.setActiveRotamer(_rotamerState[i]);
		}
}
