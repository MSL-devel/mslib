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

#ifndef CHARMMPARAMETERREADER_H
#define CHARMMPARAMETERREADER_H
// STL Includes
#include <vector>
#include <map>
#include <iostream>
#include <math.h>

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

using namespace std;

class CharmmParameterReader : public Reader {

	public:
		CharmmParameterReader();
		CharmmParameterReader(const string & _filename);
		CharmmParameterReader(const CharmmParameterReader & _par);
		~CharmmParameterReader();

		void operator=(const CharmmParameterReader & _par);

		bool read();
/*
		double getBondMinimumDistance(string _type1, string _type2);
		double getBondSpringConstant(string _type1, string _type2);
*/
		void reset();

		//VdwParam will return a vector with 4 values (Eps,Rmin,Esp14,Rmin14)
		vector<double> vdwParam(string type) const;


		/*******************************************************************
		 *  VdwParam will return a vector with 4 values 
		 * (sqrt(Eps1*Eps2),Rmin1+Rmin2,sqrt(Esp14_1+Eps14_2),Rmin14_1+Rmin14_2)
		 * NOTE: note created automatically, need to call createVdwParamPairs() first
		 *******************************************************************/
		vector<double> vdwParamPair(string type1, string type2) const;
		void createVdwParamPairs();

		//bondparam will return a vector with 2 values Kb and B0
		vector<double> bondParam(string type1, string type2) const;
	
		//Angle params will return a vector with 2 values (Ktheta, Theta0)
		vector<double> angleParam(string type1, string type2, string type3) const;

		//Ureybradley param will return a vector with 2 values (Kub, S0)
		vector<double> ureyBradleyParam(string type1, string type2, string type3) const;

		//Combined angle and ureybradley param will return a vector with 4 values (Ktheta, Theta0, Kub, S0)
		vector<double> angleAndUreyBradleyParam(string type1, string type2, string type3) const;

		//Dihedral Params will return a vector of vectors with 3 values each(Kchi, N, Delta). The vector is for each line in the Dihedral block which contains a match for these four types
		vector<vector<double> > dihedralParam(string type1, string type2, string type3, string type4) const;
		
		//Improper Params will return a vector with 2 values (Kpsi,Psi0) 
		vector<double> improperParam(string type1, string type2, string type3, string type4) const;

		//To be implemented
		vector<double> EEF1Param(string type) const; // default solvent WATER
		vector<double> EEF1Param(string type, string solvent) const;
		vector<double> IMM1Param(string type) const;
		vector<double> IMM1Param(string _type, string _solvent1, string _solvent2) const;

	private:
		void addBond(string type1, string type2, double Kb, double B0);
		void addAngle(string type1, string type2, string type3, double Ktheta, double Theta0, double Kub=0.0, double S0=0.0);
		void addDihedral(string type1, string type2, string type3, string type4, double kchi, double N, double Delta);
		void addImproper(string type1, string type2, string type3, string type4, double Kpsi, double Psi0);
		void addVdw(string type1, double Eps, double Rmin, double Esp14, double Rmin14);
		void addEEF1(string typeType, string solvent, double V, double Gref, double Gfree, double Href, double CPref, double Sigw);
		
		void setup();
		void copy(const CharmmParameterReader & _par);
		
		//Bond params will contain a vector with 2 values Kb and B0
		map<string, map<string, vector<double> > > bondParamMap;
		
		//Angle params will contain a vector with 4 values (Ktheta, Theta0, Kub, S0)
		map<string , map<string, map<string, vector<double> >  >  > angleParamMap;
		
		//Dihedral Params will contain a vector of vectors with 3 values each(Kchi, N, Delta). The vector is for each line in the Dihedral block which contains a match for these four types
		map< string, map<string , map<string, map<string, vector<vector<double> > >  >  >  > dihedralParamMap;
		

		//Improper Params will contain a vector with 2 values (Kpsi,Psi0) 
		map< string, map<string , map<string, map<string, vector<double> >  >  >  > improperParamMap;
		
		//VdwParams will contain a vector with 4 values (Eps,Rmin,Esp14,Rmin14)
		map<string,vector<double> > vdwParamMap;
		map<string,map<string, vector<double> > > vdwParamPairMap;

};

#endif

