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

#include "System.h"
#include "SysEnv.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "BaselineEnergyBuilder.h"
#include "MonteCarloOptimization.h"
#include "MonteCarloManager.h"
#include "SystemRotamerLoader.h"
#include "OnTheFlyManager.h"
#include "MslOut.h"
#include "OptionParser.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;
static MslOut MSLOUT_OPT("energyOptimization");

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
		optional.push_back("solvfile");
		optional.push_back("solvent");
		optional.push_back("hbondfile");
		optional.push_back("baselinefile");
		optional.push_back("linkedPosition");
		optional.push_back("flexNeighborDistance");
		optional.push_back("flexNeighborRotamers");

		optional.push_back("includecrystalrotamer");


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
	std::string solvfile;
	std::string solvent;
	std::string hbondfile;
	std::string baselinefile;

	std::vector<std::string> positions;
	std::vector<std::vector<std::string> > identities;
	std::vector<std::vector<int> > rotNumbers;
	std::vector<std::vector<std::string> > linkedPositions;
	double flexNeighborDistance;
	int flexNeighborRotamers;
	std::string configFile;
        bool includeCR;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;

};



struct EnergyTableOptions {

	// TODO : add cutoff options 
	// we should be able to build systems with or without cutoffs
	// also we should be able to turn on/off specific energy terms

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
	        optional.push_back("energyTable");

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
		optional.push_back("energyConfig");
		optional.push_back("maxRejections");
		optional.push_back("deltaSteps");
		optional.push_back("minDeltaE");
		defaultArgs.push_back("config");


		annealShapeMap["CONSTANT"]    = MonteCarloManager::CONSTANT;
		annealShapeMap["LINEAR"]      = MonteCarloManager::LINEAR;
		annealShapeMap["SIGMOIDAL"]   = MonteCarloManager::SIGMOIDAL;
		annealShapeMap["EXPONENTIAL"] = MonteCarloManager::EXPONENTIAL;
		annealShapeMap["SOFT"]        = MonteCarloManager::SOFT;
		annealShape = annealShapeMap["LINEAR"];

		initAlgoMap["LOWESTSELF"]     = MonteCarloOptimization::LOWESTSELF;
		initAlgoMap["RANDOM"]         = MonteCarloOptimization::RANDOM;
		initAlgoMap["QUICKSCAN"]      = MonteCarloOptimization::QUICKSCAN;
		initAlgoMap["USERDEF"]        = MonteCarloOptimization::USERDEF;
		initAlgo = initAlgoMap["LOWESTSELF"];
 
	}




	// Storage for the values of each option
	EnergyTableOptions energyOpt;
	std::string energyConfig;

	std::string configfile;
	std::string energyTable;
	std::string annealType;
	double annealStart;
	double annealEnd;
	int    numberOfAnnealCycles;
	int    numCycles;
	int maxRejections;
	double minDeltaE;
	int deltaSteps;
	std::string initAlgorithm;
	std::string initConfiguration;
	int    numStoredConfigurations;
	int    randomSeed;
	bool   DEE;

        MonteCarloManager::ANNEALTYPES annealShape;
        map<string,MonteCarloManager::ANNEALTYPES> annealShapeMap;
        MonteCarloOptimization::INITTYPE initAlgo;
        map<string, MonteCarloOptimization::INITTYPE> initAlgoMap;

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


struct AnalysisOptions {

	// Set up options here...
        AnalysisOptions(){

		// Configuration file..
		required.push_back("structureConfig");
		optional.push_back("position");
		optional.push_back("hbondParFile");
		optional.push_back("posToCheck");
		optional.push_back("aaToCheck");
		defaultArgs.push_back("config");

	}

	// Storage for the values of each option
	StructureOptions structOpt;
	std::string structureConfig;
        std::string position;
        std::string hbondParFile;
        std::string posToCheck;
        std::string aaToCheck;

	std::string configfile;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	std::vector<std::string> defaultArgs;
};

EnergyTableOptions setupEnergyTableOptions(int theArgc, char * theArgv[]);
EnergyTableOptions setupEnergyTableOptions(std::string _confFile);
MonteCarloOptions setupMonteCarloOptions(int theArgc, char * theArgv[]);
LinearProgrammingOptions setupLinearProgrammingOptions(int theArgc, char * theArgv[]);
AnalysisOptions setupAnalysisOptions(int theArgc, char * theArgv[]);

StructureOptions setupStructureOptions(std::string _confFile);

void createSystem(StructureOptions &_structOpt, MSL::System &_sys,EnergyTableOptions *_eneOpt=NULL);

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
		std::cout << "solvfile SOLVFILE\n";
		std::cout << "solvent <solvent>\n";
		std::cout << "hbondfile HBONDFILE\n";
		std::cout << "baselinefile BASELINEFILE\n";
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
	cout << "GET ROTLIB"<<endl;
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
	   // Default it to MSL_EBL (Energy-based library)
	  opt.rotlib = SYSENV.getEnv("MSL_ROTLIB"); 
	  cout << "ROTLIB: "<<opt.rotlib<<endl;
	  if (opt.rotlib == "UNDEF"){
		std::cerr << "ERROR 1111 no rotlib specified."<<std::endl;	
		exit(1111);
	  } else {
	      cerr << "WARNING rotlib defaulted to: "<<opt.rotlib<<endl;
	  }

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
		// Short cut, if ALL is found in identities string, then add all Amino Acids and Rotamer numbers
		if (tmp.size() == 1 && opt.identities[index].size() == 1 && opt.identities[index][0] == "ALL"){
		  std::vector<std::string> aaTmp(20,"");
		  aaTmp[0] = "ALA";
		  aaTmp[1] = "ARG";
		  aaTmp[2] = "ASN";
		  aaTmp[3] = "ASP";
		  aaTmp[4] = "CYS";
		  aaTmp[5] = "GLN";
		  aaTmp[6] = "GLU";
		  aaTmp[7] = "GLY";
		  aaTmp[8] = "HSD";
		  aaTmp[9] = "HSE";
		  aaTmp[10] = "ILE";
		  aaTmp[11] = "LEU";
		  aaTmp[12] = "LYS";
		  aaTmp[13] = "MET";
		  aaTmp[14] = "PHE";
		  //aaTmp[15] = "PRO";
		  aaTmp[15] = "SER";
		  aaTmp[16] = "THR";
		  aaTmp[17] = "TRP";
		  aaTmp[18] = "TYR";
		  aaTmp[19] = "VAL";
		  
		  opt.identities[index] = aaTmp;

		  std::vector<int> rotTmp(20,tmp[0]);
		  tmp = rotTmp;

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
	    opt.topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	    if (opt.topfile == "UNDEF"){
	      cerr << "ERROR 1111 topfile not defined\n";
	      exit(1111);
	    } else {
	      cerr << "WARNING topfile defaulted to: "<<opt.topfile<<endl;
	    }
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		opt.parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
		if (opt.parfile == "UNDEF"){
		  cerr << "ERROR 1111 parfile not defined\n";
		  exit(1111);
		}else {
		  cerr << "WARNING parfile defaulted to: "<<opt.parfile<<endl;
		}
	}

	opt.solvfile = OP.getString("solvfile");
	if (OP.fail()){
		opt.solvfile = "";
	        cerr << "WARNING not using solvation " <<endl;
	}

	opt.solvent = OP.getString("solvent");
	if (OP.fail() && opt.solvfile != ""){
		opt.solvent = "WATER";
	        cerr << "WARNING solvent defaulted to WATER " <<endl;
	}

	opt.hbondfile = OP.getString("hbondfile");
	if (OP.fail()){
		opt.hbondfile = SYSENV.getEnv("MSL_HBOND_PAR");
		if (opt.hbondfile == "UNDEF"){
		  cerr << "WARNING 1111 hbondfile not defined - not building hydrogen bond interactions\n";
		  opt.hbondfile = "";
		}else {
		  cerr << "WARNING hbondfile defaulted to: "<<opt.hbondfile<<endl;
		}
	}

	opt.baselinefile = OP.getString("baselinefile");
	if (OP.fail()){
		opt.baselinefile = "";
		cerr << "WARNING no baselinefile specified " <<endl;
	}

	opt.flexNeighborDistance = OP.getDouble("flexNeighborDistance");
	if (OP.fail()){
		opt.flexNeighborDistance = -1;
	} else {
	  cout << "Using flexible neighbor flag, defining neighbors by distance: "<<opt.flexNeighborDistance<<std::endl;
	}

	opt.flexNeighborRotamers = OP.getInt("flexNeighborRotamers");
	if (OP.fail()){
		if (opt.flexNeighborRotamers != -1){
			opt.flexNeighborRotamers = 10;
		} else {
			opt.flexNeighborRotamers = -1;
		}
	} else {
	  cout << "Using flexible neighbor flag, number of neighbor rotamers "<<opt.flexNeighborRotamers<<std::endl;
	}

	opt.includeCR = OP.getBool("includecrystalrotamer");
	if (OP.fail()){
	  opt.includeCR = false;
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
		std::cerr << "WARNING no distanceDielectric specified, using 'true' as default.\n";
		opt.distanceDependentElectrostatics = true;
	}

	opt.vdwScale = OP.getDouble("vdwScale");
	if (OP.fail()){
		opt.vdwScale = 1.0;
	}
	return opt;
}

EnergyTableOptions setupEnergyTableOptions(std::string _confFile){

	// Create the options
	EnergyTableOptions opt;
	

	// Parse the options
	MSL::OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	OP.readFile(_confFile);
	if (OP.fail()) {
		std::string errorMessages = "Cannot read configuration file " + _confFile + "\n";
		std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
		exit(1111);
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
		std::cout << "#     can be CONSTANT, LINEAR, EXPONENTIAL, SIGMOIDAL, SOFT, LINEAR_CYCLES, EXPONENTIAL_CYCLES, SIGMOIDAL\n";
		std::cout << "annealScheduleType EXPONENTIAL\n";
		std::cout << "annealScheduleStartTemp 1000\n";
		std::cout << "annealScheduleEndTemp   1 \n\n";
		std::cout << "numberOfAnnealCycles    1 \n\n";
		std::cout << "# Number of MC Cycles \n";
		std::cout << "numberOfCycles 10000\n\n";
		std::cout << "# Initial Rotamer Configuration set via initializationAlgorithm\n";
		std::cout << "#     can be RANDOM LOWESTSELF QUICKSCAN USERDEF\n";
		std::cout << "initializationAlgorithm LOWESTSELF\n\n";
		std::cout << "# Initial Configuration set by user, RotamerNumberForPosition1:RotamerNumberForPosition2:\n";
		std::cout << "#   Notice the last character is a ':'\n";
		std::cout << "# initializationConfiguration 0:0:0:0:\n\n";
		std::cout << "# Number of Rotamer Configurations to report , top 100 lowest energy scoring found\n";
		std::cout << "numberOfStoredConfigurations 100\n\n";
		std::cout << "# Random Seed to use, -1 means create a time-based seed\n";
		std::cout << "randomSeed 838201\n\n";
		std::cout << "# energyTable configuration file\n";
		std::cout << "energyConfig FILENAME\n";
		std::cout << "# Max number of rejections per cycle\n";
		std::cout << "maxRejections 100 \n";
		std::cout << "# Delta Steps\n";
		std::cout << "deltaSteps 100 \n";
		std::cout << "# minDeltaEnergy\n";
		std::cout << "minDeltaE 0.1 \n";
		exit(0);
        
	}

	opt.minDeltaE = OP.getDouble("minDeltaE");
	if (OP.fail()){
		opt.minDeltaE = 0.1;
	}

	opt.deltaSteps = OP.getInt("deltaSteps");
	if (OP.fail()){
		opt.deltaSteps = opt.numCycles / 10;
	}

	opt.maxRejections = OP.getInt("maxRejections");
	if (OP.fail()){
		opt.maxRejections = opt.numCycles / 10;
	}

	opt.energyConfig = OP.getString("energyConfig");
	if (!OP.fail()){
		opt.energyOpt = setupEnergyTableOptions(opt.energyConfig);   
	} 

	opt.debug = OP.getBool("debug");

	opt.energyTable = OP.getString("energyTable");
	if (OP.fail()){
	  opt.energyTable = "";
	}


	opt.annealType = OP.getString("annealScheduleType");
	if (OP.fail()){
		opt.annealType = "LINEAR";
	}
	map<string,MonteCarloManager::ANNEALTYPES>::iterator it;
	it = opt.annealShapeMap.find(opt.annealType);
	if (it == opt.annealShapeMap.end()){
	  cerr << "ERROR 1111 anneal shape "<<opt.annealType<<" is not known.\n";
	  exit(1111);
	} else {
	  opt.annealShape = it->second;
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
	map<string,MonteCarloOptimization::INITTYPE>::iterator it2;
	it2 = opt.initAlgoMap.find(opt.initAlgorithm);
	if (it2 == opt.initAlgoMap.end()){
	  cerr << "ERROR 1111 init algorithm: "<<opt.initAlgorithm<<" is unknown.\n";
	  exit(1111);
	} else {
	  opt.initAlgo = it2->second;
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


	cout << OP<<std::endl;
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


	cout << OP<<std::endl;
	return opt;
}
AnalysisOptions setupAnalysisOptions(int theArgc, char * theArgv[]) {
	// Create the options
	AnalysisOptions opt;
    

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
	    
		std::cout << "# Options for examineSideChains\n\n";
		std::cout << "structureConfig STRUCT_CONFIG_FILE\n";
		std::cout << "position A,1\n";

		exit(0);
        
	}

	opt.position = OP.getString("position");
	if (OP.fail()){
	  cerr << "ERROR 1111 position not specified\n";
	  exit(1111);
	}

	opt.structureConfig = OP.getString("structureConfig");
	if (!OP.fail()){
		opt.structOpt = setupStructureOptions(opt.structureConfig);
	} else {
	  std::cerr << "ERRROR 1111 structureConfig not specified\n";
	  exit(1111);
	}

	opt.hbondParFile = OP.getString("hbondParFile");
	if (OP.fail()){
	  opt.hbondParFile = SYSENV.getEnv("MSL_HBOND_PAR");
	  if (opt.hbondParFile == "UNDEF"){
	    cerr << "ERROR 1111 hbond parfile not defined\n";
	    exit(1111);
	  }else {
	    cerr << "WARNING hbondParFile defaulted to: "<<opt.hbondParFile<<endl;
	  }
	}


	opt.posToCheck = OP.getString("posToCheck");
	if (OP.fail()){
	  opt.posToCheck = "A_100";
	}
	opt.aaToCheck = OP.getString("aaToCheck");
	if (OP.fail()){
	  opt.aaToCheck = "HSD";
	}
	cout << OP<<std::endl;
	return opt;
}

void createSystem(StructureOptions &_structOpt, MSL::System &_sys,EnergyTableOptions *_eneOpt){


	// Read-in PDB as initialSystem
	MSLOUT_OPT.stream() << "Read in PDB: "<<_structOpt.pdb<<std::endl;
	MSL::System initialSystem;
	initialSystem.readPdb(_structOpt.pdb);

	// Build Map of Variable Positions
	cout << "Build Designable PolymerSequence"<<std::endl;
	std::map<std::string, int> variablePositionMap;
	uint psize = _structOpt.positions.size();

	// Assign the first indices to the variable positions specified explicitly in the file
	for (uint v = 0; v < psize;v++){
		MSLOUT_OPT.stream() << "Position "<<v<<" '"<<_structOpt.positions[v]<<"'"<< std::endl;

		std::string chain = "";
		int resnum = 0;
		std::string icode = "";
		bool OK = MslTools::parsePositionId(_structOpt.positions[v], chain, resnum, icode, 0);

		if (!OK){
			std::cerr << "ERROR 2222: Position "<<v<<" '"<<_structOpt.positions[v]<<"' does not have CHAIN,RESIDUE format\n";
			exit(2222);
		}

		if (initialSystem.positionExists(chain, resnum, icode)){
		        MSLOUT_OPT.stream() << "Adding position: "<<_structOpt.positions[v]<<" to variable position map."<<endl;
			variablePositionMap[_structOpt.positions[v]] = v;					
		}
	}

	if (_structOpt.flexNeighborDistance != -1){
		for (uint v = 0; v < psize;v++){
			MSLOUT_OPT.stream() << "Position "<<v<<" '"<<_structOpt.positions[v]<<"'"<< std::endl;

			std::string chain = "";
			int resnum = 0;
			std::string icode = "";
			bool OK = MslTools::parsePositionId(_structOpt.positions[v], chain, resnum, icode, 0);

			if (!OK){
				std::cerr << "ERROR 2222: Position "<<v<<" '"<<_structOpt.positions[v]<<"' does not have CHAIN,RESIDUE format\n";
				exit(2222);
			}

			if (initialSystem.positionExists(chain, resnum, icode)){

				// If add neighbors automatically
				MSL::Residue &res = initialSystem.getResidue(chain, resnum, icode); // get the current identity

				cout << "Finding Neighbors"<<std::endl;
				std::vector<int> neighbors = res.findNeighbors(_structOpt.flexNeighborDistance);
				cout << "Neighbors are: "<<neighbors.size()<<std::endl;

				
				for (uint n = 0; n < neighbors.size();n++){
					MSL::Residue &neighbor = initialSystem.getResidue(neighbors[n]);
					cout << "\t"<<neighbors[n]<<" "<<neighbor.getResidueNumber()<<std::endl;
				}
				for (uint n = 0; n < neighbors.size();n++){

					MSL::Residue &neighbor = initialSystem.getResidue(neighbors[n]);


					if (neighbor.getResidueName() == "ALA" || neighbor.getResidueName() == "PRO" || neighbor.getResidueName() ==  "GLY"){
						continue;
					}

					// Skip if in std::map already...
					if (variablePositionMap.find(neighbor.getPositionId()) != variablePositionMap.end()){
						continue;
					}

					// Add to _opt.positions, _opt.rotNumbers, _opt.identities and variablePositionMap
				//	char tmpStr[80];
				//	sprintf(tmpStr, "%1s_%d",neighbor.getChainId().c_str(),neighbor.getResidueNumber());
				//	_opt.positions.push_back(tmpStr);
					//cout << "Pushing back " << neighbor.getPositionId() << endl;
					_structOpt.positions.push_back(MslTools::getPositionId(neighbor.getChainId(), neighbor.getResidueNumber(), neighbor.getResidueIcode()));
					_structOpt.identities.push_back(std::vector<std::string>(1,neighbor.getResidueName()));

					if (neighbor.getResidueName() == "SER" || neighbor.getResidueName() == "THR" || neighbor.getResidueName() == "CYS" || neighbor.getResidueName() == "VAL"){
						_structOpt.rotNumbers.push_back(std::vector<int>(1,10));
					} else {
						_structOpt.rotNumbers.push_back(std::vector<int>(1,_structOpt.flexNeighborRotamers));
					}

					
					variablePositionMap[neighbor.getPositionId()] = _structOpt.positions.size() - 1;
					
				} // END neighbors.size()

			} else {
			  cerr << "ERROR 1231 Position: "<<_structOpt.positions[v]<<" does not exist in the PDB."<<endl;
			  exit(1231);
			} // IF positionExists in initialSystem
		} // END _opt.positions.size()
	} // IF flexibleNeighbors


	// Add linked positions, first position in each linkedPositions[INDEX] is the "master" position
	// This is needed to build a proper PolymerSequence..
	initialSystem.setLinkedPositions(_structOpt.linkedPositions);
	
	// Build PolymerSequence.
	cout << "build it"<<std::endl;
	MSL::PolymerSequence pseq(initialSystem, variablePositionMap, _structOpt.identities);
	cout << "Done."<<std::endl;
	cout << "PolymerSequence: "<<std::endl<<pseq.toString()<<std::endl;

	// Build a new system from polymer sequence and create energySet from energy terms
	cout << "Build Charmm System"<<std::endl;
	cout << "Use toppar " << _structOpt.topfile << ", " << _structOpt.parfile << std::endl;
	MSL::CharmmSystemBuilder CSB(_sys,_structOpt.topfile,_structOpt.parfile,_structOpt.solvfile);
	if(_structOpt.solvfile != "") {
		cout << "Use solvation " << _structOpt.solvfile << ", " << _structOpt.solvent << endl;
		CSB.setSolvent(_structOpt.solvent);
	}


	// Check for type of energy calculation...
	if (_eneOpt != NULL){
	  CSB.setDielectricConstant(_eneOpt->dielectric);
	  CSB.setUseRdielectric(_eneOpt->distanceDependentElectrostatics);
	  CSB.setVdwRescalingFactor(_eneOpt->vdwScale);
	}
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).
	//cout << pseq << endl;


	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	int numAssignedAtoms = _sys.assignCoordinates(initialSystem.getAtomPointers(),false);
	fprintf(stdout, "Assigned %8d atoms\n",numAssignedAtoms);

	// Build the all atoms without coordinates (not in initial PDB)
	cout <<"Build Atoms"<<std::endl;
	_sys.buildAllAtoms(); 

	CSB.updateNonBonded(9.0, 10.0, 11.0);

	HydrogenBondBuilder hb;
	if(_structOpt.hbondfile != "") {
		hb.setSystem(_sys);
		if(!hb.readParameters(_structOpt.hbondfile)) {
			cerr << "ERROR 1234 Unable to read hbondfile " << _structOpt.hbondfile << endl;
			exit(1234);
		}
		hb.buildInteractions(-1.0);
	}

	BaselineEnergyBuilder bb;
	if(_structOpt.baselinefile != "") {
		bb.setSystem(_sys);
		if(!bb.readParameters(_structOpt.baselinefile)) {
			std::cerr << "Unable to read baselineenergies " <<  _structOpt.baselinefile << endl;
			exit(0);
		}

		if(!bb.buildInteractions()) {
			std::cerr << "Unable to build baselineInteractions " <<  endl;
			exit(0);
		}
	}

	// Set linked positions in our new sequence-built system
	_sys.setLinkedPositions(_structOpt.linkedPositions);


	// Write initially built file.
	std::string filename = "/tmp/initialBuild.pdb";
	cout << "Write pdb " << filename << std::endl;
	MSL::PDBWriter writer;
	writer.open(filename);
	if (!writer.write(_sys.getAtomPointers())) {
		std::cerr << "Problem writing " << filename << std::endl;
	}
	writer.close();

	
	// Build rotamers
	MSLOUT_OPT.stream() << "Read rotamer library " << _structOpt.rotlib << " and load rotamers"<<std::endl;	
	MSL::SystemRotamerLoader sysRot(_sys, _structOpt.rotlib);
	std::map<std::string,int >::iterator it1;
	for (uint i = 0; i < _sys.positionSize();i++){
		MSL::Position &pos = _sys.getPosition(i);

		std::string chainId = pos.getChainId();
		int resNum          = pos.getResidueNumber();
		string resIcode     = pos.getResidueIcode();

		// If a "slave" position, then take its "master" positions chainId and residue number.
		string posId = pos.getPositionId();
		if (pos.getLinkedPositionType() == MSL::Position::SLAVE){
			chainId    = pos.getLinkedPositions()[0]->getChainId();
			resNum     = pos.getLinkedPositions()[0]->getResidueNumber();
			resIcode   = pos.getLinkedPositions()[0]->getResidueIcode();
			cout << " FOUND LINKED POS: "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" to "<<chainId<<" " <<resNum<<std::endl;
			
			// redefine posid to master posid...
			posId = MslTools::getPositionId(chainId, resNum, resIcode);
		}	
		
		int index = -1;
		
		it1 = variablePositionMap.find(posId);
		if (it1 != variablePositionMap.end()){
		        index = it1->second;
		}
		
		if (index != -1){
			for (uint j = 0; j < _structOpt.identities[index].size();j++){
				MSLOUT_OPT.stream() << "Loading "<<_structOpt.rotNumbers[index][j]<<" rotamers of type "<<_structOpt.identities[index][j]<<" at postion "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName();
				if (pos.getLinkedPositionType() == MSL::Position::SLAVE){
					cout << " LINKED TO "<< chainId<<" "<<resNum;
				}
				cout <<std::endl;

				// Test for position specific library...
				char posSpecificLib[80];
				sprintf(posSpecificLib,"POS_%1s_%04d",pos.getChainId().c_str(),pos.getResidueNumber());

				// Load rotamers
				if (sysRot.getRotamerLibrary()->libraryExists((std::string)posSpecificLib)) {

					int lastRotamerIndex = _structOpt.rotNumbers[index][j];
					if (sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_structOpt.identities[index][j]) < _structOpt.rotNumbers[index][j]){
						std::cerr << "WARNING "<<pos.getChainId()<<" "<<pos.getResidueNumber()<<" "<<pos.getResidueName()<<" has a specific rotamer library but not enough rotamers it has: "<<sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_structOpt.identities[index][j])<<std::endl;
						lastRotamerIndex =sysRot.getRotamerLibrary()->size((std::string)posSpecificLib,_structOpt.identities[index][j]);
					}
					//sysRot.loadRotamers(&pos, posSpecificLib, _opt.identities[index][j], 0, lastRotamerIndex); 
					sysRot.loadRotamers(&pos, _structOpt.identities[index][j], lastRotamerIndex, posSpecificLib,_structOpt.includeCR); 
				} else {
					//sysRot.loadRotamers(&pos, "BALANCED-200", _opt.identities[index][j], 0, _opt.rotNumbers[index][j]-1); 
				        sysRot.loadRotamers(&pos, _structOpt.identities[index][j], _structOpt.rotNumbers[index][j], "",_structOpt.includeCR); 
				}
				
				
			}
		}
	}

	cout << "Create system looks like this: "<<endl<< _sys.toString()<<endl;

}


void changeRotamerState(StructureOptions &_opt, MSL::System &_sys, std::vector<int> &_rotamerState){


		for (uint i = 0; i < _opt.positions.size();i++){
			//MSLOUT_OPT.stream() << "Working on absres: "<<_opt.positions[i]<<std::endl;

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
			//cout << "AA -Position: "<<_opt.positions[i]<< " "<<pos.getNumberOfIdentities() << " "<< pos.getTotalNumberOfRotamers()<<" changing to "<<_rotamerState[i];
			pos.setActiveRotamer(_rotamerState[i]);
			//cout << " identity is: "<<pos.getCurrentIdentity().getResidueName()<<endl;
		}
}
