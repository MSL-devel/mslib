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
#include <string>
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "Timer.h"
#include "CharmmEnergyCalculator.h"
#include "metalRotamers.h"

using namespace std;

using namespace MSL;
#include "energyOptimizations.h"


// Global objects
Timer t;
double startTime = 0.0;

// MslOut 
static MslOut MSLOUT("metalRotamers");

// SysEnv - ALREADY DEFINED IN energyOptimizations.h
//static SysEnv SYSENV;




int main(int argc, char *argv[]){

	// Option Parser
	cout << "Setup Options"<<endl;
	Options opt = setupOptions(argc, argv);

	StructureOptions structOpt;
	structOpt.pdb        = opt.pdb;
	structOpt.rotlib     = opt.rotlib;
	structOpt.parfile    = opt.parfile;
	structOpt.topfile    = opt.topfile;
	structOpt.positions  = opt.variable_positions;
	structOpt.identities = opt.variable_identities;
	structOpt.rotNumbers = opt.variable_rotamers;
	structOpt.flexNeighborDistance = -1;

	// Create a system from the structural input options
	cout << "Create System"<<endl;
	startTime = t.getWallTime();
	System sys;
	createSystem(structOpt, sys);
	cout << "Built system after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	startTime = t.getWallTime();

	CharmmEnergyCalculator cec(opt.parfile);
	cec.extractBondedInteractions(*sys.getEnergySet());

	/*
	  Fixed Position (HIS) - Zn - Variable Position (HIS)
	    Build all rotamers on Variable Position HIS
	    Check geometry.
	    Keep list of potential_rotamers.
	    Evaluate energies of potential_rotamers, write out top opt.numModels
	  
	 */

	for (uint i = 0; i < opt.variable_positions.size();i++){
	  Position &var_pos = sys.getPosition(opt.variable_positions[i]);
	  Position &fix_pos = sys.getPosition(opt.fixed_positions[i]);
	  Position &zn_pos  = sys.getPosition(opt.zn_positions[i]);
	  if (!zn_pos.atomExists("ZN")){
	    MSLOUT.stream() << "ERROR 5455 ZN position: "<<zn_pos.toString()<<" has no ZN atom"<<endl;
	    continue;
	  }

	  Atom &ZN = zn_pos.getAtom("ZN");
	  priority_queue< pair<double,int>, vector< pair<double,int> >, less<pair<double,int> > > potential_rotamers;


	  for (uint j = 0; j < var_pos.getTotalNumberOfRotamers();j++){
	    var_pos.setActiveRotamer(j);
	    if (! (var_pos.getResidueName() == "HSD" || var_pos.getResidueName() == "HSE") ) {
	      MSLOUT.stream() << "Skipping rotamer "<<j<<" because identity is: "<<var_pos.getResidueName()<<"."<<endl;
	      continue;
	    }
	    
	    MSLOUT.stream() << "Rotamers: "<<j<<endl;

	    
	    if (!checkMetalSiteGeometry(fix_pos,var_pos,ZN)){

	      // ERROR MSG
	      MSLOUT.stream() << "SKIP Position: "<<opt.variable_positions[i]<<" Rotamer: "<<j<<" does not meet geometry criteria"<<endl;
	      continue;

	    }

	    MSLOUT.stream() << "GOOD Position: "<<opt.variable_positions[i]<<" Rotamer: "<<j<<" does meet geometry criteria with an energy of : ";

	    // Check energy of this rotamer
	    double energy = cec.calculateBackgroundEnergy(sys, i, j);
	    MSLOUT.stream()<< energy<<endl;

	    // Add to priority_queue for this position (energy,index)
	    potential_rotamers.push(pair<double,int>(energy,j));
	  }

	  // If priority-queue for this position is empty ERROR
	  if (potential_rotamers.empty()){
	    cerr << "ERROR 4444 No rotamers found to satisfy conditions."<<endl;
	    exit(4444);
	  }

	  // Write out models for the top numModels from priority_queue here
	  int numModels = 0;
	  int lowestEnergyRotamer = potential_rotamers.top().second;
	  while (!potential_rotamers.empty()){
	    numModels++;
	    if (numModels == opt.numModels) break;


	    pair<double,int> model = potential_rotamers.top();
	    potential_rotamers.pop();

	    var_pos.setActiveRotamer(model.second);

	    // Write out file
	    string fname = MslTools::stringf("%s_%s_%04d_%05d.pdb",MslTools::getFileName(opt.pdb).c_str(),var_pos.getChainId().c_str(),var_pos.getResidueNumber(),numModels);
	    MSLOUT.stream() << "Writing: "<<fname<<endl;
	    sys.writePdb(fname);
	  }

	  // Set lowest-energy rotamer for this position in the priority_queue
	  var_pos.setActiveRotamer(lowestEnergyRotamer);
	}


	MSLOUT.stream() << "Done."<<endl;


}

/*

	       NE2-chelating:
	          NE2-Zn 2.0 - 2.6
		  NE2-Zn-NE2 80-100
		  ND1-Zn-ND1 70-110

	       ND1-chelating:
	          ND1-Zn 2.0 - 2.6
		  ND1-Zn-ND1 80-100
		  NE2-Zn-NE2 70-110		  

	       Mixed-chelating1 (first = fixed [ND1-chelating], second = variable [NE2-chelating] ):
	          ND1-Zn; NE2-Zn 2.0 - 2.6
		  ND1-Zn-NE2 80-100
		  NE2-Zn-ND1 70-110		  

	       Mixed-chelating2 (first = fixed [NE2-chelating], second = variable [ND1-chelating] ):
	          NE2-Zn; ND1-Zn 2.0 - 2.6
		  NE2-Zn-ND1 80-100
		  ND1-Zn-NE2 70-110		  

 */
bool checkMetalSiteGeometry(Position &_fix, Position &_var, Atom &_zn){
  
  if (!(_fix.atomExists("NE2") && _fix.atomExists("ND1"))){
    MSLOUT.stream() << "ERROR fix position: "<<_fix.toString()<<" doesn't have NE2 or ND1"<<endl;
    return false;
  }

  if (!(_var.atomExists("NE2") && _var.atomExists("ND1"))){
    MSLOUT.stream() << "ERROR var position: "<<_var.toString()<<" doesn't have NE2 or ND1"<<endl;
    return false;
  }

  // Check two angles to ZN from var.
  double ring_angle1 = _var.getAtom("NE2").angle(_var.getAtom("ND1"),_zn);
  double ring_angle2 = _var.getAtom("ND1").angle(_var.getAtom("NE2"),_zn);
  if ( !   
      (
       (ring_angle1 > 120 && ring_angle1 < 240) ||
       (ring_angle2 > 120 && ring_angle2 < 240) 
       )
    ){
    MSLOUT.stream() << "ERROR skipping because of bad ring angles: "<<ring_angle1<<","<<ring_angle2<<endl;
    return false;
  }
      
	 
  // NE2-chelating
  double dist1  = _var.getAtom("NE2").distance(_zn);
  double dist2  = _fix.getAtom("NE2").distance(_zn);
  double angle1 = _fix.getAtom("NE2").angle(_zn,_var.getAtom("NE2"));
  double angle2 = _fix.getAtom("ND1").angle(_zn,_var.getAtom("ND1"));
  if (
      dist1   > 2.0 && 
      dist1   < 2.6 &&
      dist2   > 2.0 && 
      dist2   < 2.6 &&
      angle1 > 80 &&
      angle1 < 100 &&
      angle2 > 70 &&
      angle2 < 110
      ){
    MSLOUT.stream() << "Metal site is NE2-chelating "<<_var.toString()<<" and "<<_fix.toString()<<endl;
    return true;
  }

  // ND1-chelating
  dist1  = _var.getAtom("ND1").distance(_zn);
  dist2  = _fix.getAtom("ND1").distance(_zn);
  angle1 = _fix.getAtom("ND1").angle(_zn,_var.getAtom("ND1"));
  angle2 = _fix.getAtom("NE2").angle(_zn,_var.getAtom("NE2"));
  if (
      dist1   > 2.0 && 
      dist1   < 2.6 &&
      dist2   > 2.0 && 
      dist2   < 2.6 &&
      angle1 > 80 &&
      angle1 < 100 &&
      angle2 > 70 &&
      angle2 < 110
      ){
    MSLOUT.stream() << "Metal site is ND1-chelating "<<_var.toString()<<" and "<<_fix.toString()<<endl;
    return true;
  }  

  // Mixed-chelating1
  dist1  = _fix.getAtom("ND1").distance(_zn);
  dist2  = _var.getAtom("NE2").distance(_zn);
  angle1 = _fix.getAtom("ND1").angle(_zn,_var.getAtom("NE2"));
  angle2 = _fix.getAtom("NE2").angle(_zn,_var.getAtom("ND1"));
  if (
      dist1   > 2.0 && 
      dist1   < 2.6 &&
      dist2   > 2.0 && 
      dist2   < 2.6 &&
      angle1 > 80 &&
      angle1 < 100 &&
      angle2 > 70 &&
      angle2 < 110
      ){
    MSLOUT.stream() << "Metal site is Mixed-chelating1 "<<_var.toString()<<" and "<<_fix.toString()<<endl;
    return true;
  }  

  // Mixed-chelating2
  dist1  = _fix.getAtom("NE2").distance(_zn);
  dist2  = _var.getAtom("ND1").distance(_zn);
  angle1 = _fix.getAtom("NE2").angle(_zn,_var.getAtom("ND1"));
  angle2 = _fix.getAtom("ND1").angle(_zn,_var.getAtom("NE2"));

  if (
      dist1   > 2.0 && 
      dist1   < 2.6 &&
      dist2   > 2.0 && 
      dist2   < 2.6 &&
      angle1 > 80 &&
      angle1 < 100 &&
      angle2 > 70 &&
      angle2 < 110
      ){
    MSLOUT.stream() << "Metal site is Mixed-chelating2 "<<_var.toString()<<" and "<<_fix.toString()<<endl;
    return true;
  }  

  return false;
}



void cleanExit(int sig) {

	cout << "*************SIGNAL CAUGHT*****************\n";
	cout << "\tSIG"<<sig<<endl;


	cout << "GoodBye."<<endl;
	exit(0);
}





Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	// Parse the options
	MSL::OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.readArgv(theArgc, theArgv);
	opt.configFile = OP.getString("configfile");
	if (OP.fail()){
	  cerr << "YOU DID NOT USE A CONFIG FILE\n";
	}
	

	OP.readFile(opt.configFile);
	if (OP.fail()) {
		std::string errorMessages = "Cannot read configuration file '" + opt.configFile + "'\n";
		std::cerr << "ERROR 1111 "<<errorMessages<<std::endl;
		exit(1111);
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
		std::cout << "variablePosition     CHAIN_RESIDUE\n";
		std::cout << "variableRotamers     NUM\n";
		std::cout << "fixedPosition        CHAIN_RESIDUE\n";
		std::cout << "znPosition        CHAIN_RESIDUE\n";
		std::cout << "numModels INT\n";
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
		std::string position  = OP.getString("variablePosition",index);
		if (OP.fail()) {
			break;
		}
		int rotamers          = OP.getInt("variableRotamer", index);
		if (OP.fail()) {
			break;
		}

		std::string fixed = OP.getString("fixedPosition",index);
		if (OP.fail()) {
			break;
		}
		std::string zn    = OP.getString("znPosition",index);
		if (OP.fail()) {
			break;
		}
		opt.variable_positions.push_back(position);
		opt.fixed_positions.push_back(fixed);
		opt.zn_positions.push_back(zn);

		vector<string> ids;
		ids.push_back("HSD");
		ids.push_back("HSE");
		opt.variable_identities.push_back(ids);

		vector<int> rots;
		rots.push_back(rotamers);
		rots.push_back(rotamers);

		opt.variable_rotamers.push_back(rots);

		index++;
	}
	if (opt.variable_positions.size() == 0) {
		std::cerr << "The number of entries for option 'position' needs to be at least 1\n";
	}
	
	if (opt.variable_positions.size() != opt.variable_identities.size()) {
		std::cerr << "The number of entries for options positions and identities need to be identical. Since its not, the program will take WT Identities\n";
	}

	if (opt.variable_positions.size() != opt.variable_rotamers.size()) {
		std::cerr << "The number of entries for options positions and rotNumbers need to be identical\n";
	} else {
		for (unsigned i=0; i<opt.variable_identities.size(); i++) {
			if (opt.variable_rotamers[i].size() != opt.variable_identities[i].size()) {
				std::cerr<< "The series of number of rotamers \"rotNumbers\" needs to be identical to the series of identities at position " + opt.variable_positions[i] + "\n";
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

	opt.numModels = OP.getInt("numModels");
	if (OP.fail()){
	  opt.numModels = 5;
	}
	return opt;
}
