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

#ifndef PHIPSISTATISTICS
#define PHIPSISTATISTICS


//MSL Includes
#include "Residue.h"
#include "RandomNumberGenerator.h"

// STL Includes
#include <stdio.h>

namespace MSL { 
class PhiPsiStatistics {
	public:
		PhiPsiStatistics();
		PhiPsiStatistics(const PhiPsiStatistics &_phiPsiStat);
		~PhiPsiStatistics();

		void operator=(const PhiPsiStatistics &_phiPsiStat);
		void addStatisitics(std::string _residueType, std::string _phiBin, std::string _psiBin,int _count);

		int operator()(std::string _key);
		void setGridSize(double _gridSize);
		double getGridSize();
		int getCounts(std::string resName, double phi, double psi);
		int getCounts(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getProbability(std::string &resName, double phi, double psi);
		double getProbability(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getProbabilityAll(double phi, double psi);
		double getProbabilityAll(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getPropensity(std::string &resName, double phi, double psi);
		double getPropensity(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getFreqInBin(std::string _resName, double phi, double psi);

		static double getPhi(const Residue &nMinus1, const Residue &n);
		static double getPsi(const Residue &n, const Residue &nPlus1);
		static double getOmega(const Residue &n, const Residue &nPlus1);
		void computeTotalCounts();

		/*
			Functions to get a random phi or psi , given the probablity distribution that has been read in.
		        The function gets a random pair of phi,psi angles.
		 */
		
		std::pair<double,double> getRandomPhiPsi(std::string _resType);

		std::map<std::string,int>  getPhiPsiCounts() const { return phiPsiTable; }
		double getPhiPsiBin(double _angle);
	private:
     

		void copy(const PhiPsiStatistics &_phiPsiStat);
                double gridSize;

		std::map<std::string,int> phiPsiTable;

		struct PhiPsiRNG {
		    RandomNumberGenerator rng;
		    std::vector<double> counts;
		    std::vector<std::pair<double,double> > phiPsiValues;

		    /* NO LONGER NEEDED, DEFAULTED IN THE OBJECT
		    PhiPsiRNG(){
			rng.setRNGType("knuthran2002");
			rng.setRNGTimeBasedSeed();
		    }
		    */

		    std::pair<double,double>& getRandomAngles(){
			    int randomIndex = rng.getRandomDiscreteIndex();
			    if (randomIndex < 0 || randomIndex > phiPsiValues.size()){
				    std::cerr << "ERROR 3432 PhiPsiStatistics::PhiPsiRNG::getRandomAngles(), index is bigger than data size..."<<randomIndex<<" size: "<<phiPsiValues.size()<<"\n";
				exit(3432);
			    }
			    return phiPsiValues[randomIndex];
		    }
		};

		std::map<std::string,PhiPsiRNG *> phiPsiRandom;


};

inline	void   PhiPsiStatistics::setGridSize(double _gridSize){
  gridSize = _gridSize;
}
inline 	double PhiPsiStatistics::getGridSize(){
  return gridSize;
}	

}

#endif
