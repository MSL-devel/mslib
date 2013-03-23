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

#include "SysEnv.h"

using namespace MSL;
using namespace std;




SysEnv::SysEnv(){
	setup();
}

SysEnv::~SysEnv(){
}


void SysEnv::setup(){
	// Defaults should be set here
	defaultString = "";
	undefinedString = "UNDEF";

	string mslHome = getEnv("MSL_DIR");
	if (mslHome == undefinedString || mslHome == defaultString) {
		mslHome = getEnv("HOME");
		if (mslHome == undefinedString || mslHome == defaultString) {
			mslHome = "/usr/local/mslib";
		} else {
			mslHome = mslHome+"/mslib";
		}
	}
	env["MSL_DIR"]         = mslHome;
	
	// default locations for commonly used input files 

	// charmm toppar (standard charmm 22)
	//env["MSL_CHARMM_TOP"] = env["MSL_DIR"]+"/toppar/charmm/top_all27_prot_lipid.inp";
	//env["MSL_CHARMM_PAR"] = env["MSL_DIR"]+"/toppar/charmm/par_all27_prot_lipid.inp";
	env["MSL_CHARMM_TOP"] = env["MSL_DIR"]+"/toppar/charmm22.top";
	env["MSL_CHARMM_PAR"] = env["MSL_DIR"]+"/toppar/charmm22.par";

	// scwrl 4 hydrogen bond (canonical, no CA hbond)
	env["MSL_HBOND_PAR"]   = env["MSL_DIR"]+"/toppar/scwrl4hb/par_hbond_1.txt";

	// scwrl 4 hydrogen bond (canonical + CA hbond)
	env["MSL_HBOND_CA_PAR"]   = env["MSL_DIR"]+"/toppar/scwrl4hb/par_hbond_CA_1.txt";
	
	// balanced rotamer library
	env["MSL_PDB_TOP"]     = env["MSL_DIR"]+"/toppar/top_pdb2.3_noH.inp";
	env["MSL_PDBH_TOP"]     = env["MSL_DIR"]+"/toppar/top_pdb2.3_H.inp";

	// the default is now the Energy-Based rotamer library
	env["MSL_ROTLIB"]        = env["MSL_DIR"]+"/library/EBL_11-2011_CHARMM22.txt";
//	env["MSL_ROTLIB"]     = env["MSL_DIR"]+"/rotlib/balanced/rotlib-balanced-200.txt";
//	env["MSL_EBL"]        = env["MSL_DIR"]+"/library/EBL_11-2011_CHARMM22.txt";

	// location of example files
	env["MSL_EXAMPLE_FILE_DIR"]     = env["MSL_DIR"]+"/exampleFiles";

	env["MSL_BBQ_TABLE"]         = env["MSL_DIR"]+"/tables/PiscesBBQTable.txt";
	env["MSL_PDBFRAG_TABLE_LINUX32"]     = env["MSL_DIR"]+"/tables/nr1000.fragdb";
	env["MSL_PDBFRAG_TABLE_MAC32"]     = env["MSL_DIR"]+"/tables/nr1000.fragdb";

	env["MSL_PHIPSI_TABLE"] = env["MSL_DIR"]+"/tables/phiPsiCounts.txt";

	// location of degree of freedom definition files
	env["MSL_PDB_2_3_DOF"] = env["MSL_DIR"]+"/toppar/pdb_2.3_DegOfFreedoms.txt";
	env["MSL_CHARMM_22_DOF"] = env["MSL_DIR"]+"/toppar/CHARMM_22_DegOfFreedoms.txt";

	// location of salt bridge propensity table
	env["MSL_SB_PROP_TABLE"] = env["MSL_DIR"]+"/tables/sb_prop_table_12_14_12.txt";

	// location of prosite data
	env["MSL_PROSITE"] = env["MSL_DIR"]+"/tables/prosite.dat";
}

bool SysEnv::setEnv(string &_var, string &_value){
	env[_var] = _value;
	return true;
}

bool SysEnv::addEnvVariable(string &_var){
	env[_var] = defaultString;
	return true;
}

bool SysEnv::isDefined(string &_var){
	envIt = env.find(_var);
	char *envVar = NULL;
	envVar = getenv(_var.c_str());

	// if stored in SysEnv or in users environment, return true
	return (envIt != env.end() || envVar != NULL);
}

std::string SysEnv::getEnv(std::string const &_var){
	
	// First look for variable in the users environment..
	char *envVar = NULL;
	envVar = getenv(_var.c_str());

	// It not found in the environment, use our defaults..
	if (envVar == NULL){

		// Make sure the variable has been defined
		envIt = env.find(_var);

		// Return default string when variable is not in users environment or we have not defined a default
		if (envIt == env.end()){
			return undefinedString;
		} else {
			return envIt->second;
		}

	} else {
		// Store environment variable.
		env[_var] = (string)envVar;
	}

	return env[_var];
}
