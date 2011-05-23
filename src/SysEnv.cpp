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

	std::string mslHome = getEnv("MSLHOME");
	if (mslHome == undefinedString || mslHome == defaultString){
		char *envVar = NULL;
		envVar = getenv("HOME");
		if (envVar == NULL){
			mslHome = "/usr/local/mslib";
		} else {
			mslHome = ((string)envVar)+"/mslib";
		}
	}
	
	env["MSLHOME"]      = mslHome;

	env["HBONDPARDIR"]   = env["MSLHOME"]+"/hbondPar";
	env["HBONDPAR"]      = env["HBONDPARDIR"]+"/canonical_hbond_1.txt";

	env["CHARMMDIR"]      = env["MSLHOME"]+"/charmmTopPar";
	env["CHARMMPAR"]      = env["CHARMMDIR"]+"/par_all27_prot_lipid.inp";
	env["CHARMMTOP"]      = env["CHARMMDIR"]+"/top_all27_prot_lipid.inp";
	env["CHARMMTOP_EEF1"] = env["CHARMMDIR"]+"/toph19_eef1.1.inp";
	env["CHARMMPAR_EEF1"] = env["CHARMMDIR"]+"/param19_eef1.1.nowildcards.inp";
	env["CHARMMSOL"]      = env["CHARMMDIR"]+"/solvpar.inp";

	env["ROTLIBDIR"]    = env["MSLHOME"]+"/rotlib";



	env["ROTLIB"]       = env["ROTLIBDIR"]+"/balanced/rotlib-balanced-200.txt";

	
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
