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

#ifndef SCMF_H
#define SCMF_H

#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sys/stat.h>

#include "MslTools.h"
#include "MonteCarloManager.h"

using namespace std;
using namespace MSL;

/*! \brief Self Consitent Mean Field class
 */
namespace MSL {

class SelfConsistentMeanField {
	public:
		SelfConsistentMeanField();
		SelfConsistentMeanField(vector<vector<double> > & selfEnergies, vector<vector<vector<vector<double> > > > & pairEnergies);
		SelfConsistentMeanField(double & fixEnergy, vector<vector<double> > & selfEnergies, vector<vector<vector<vector<double> > > > & pairEnergies);
		SelfConsistentMeanField(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		SelfConsistentMeanField(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		~SelfConsistentMeanField();

		void setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies);
		void setEnergyTables(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies);
		void setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		void setEnergyTables(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		void setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines);

		//void getExternalRNG(RandomNumberGenerator * _pExternalRNG);

		void cycle();
		std::vector<unsigned int> runMC(double _startingTemperature, double _endingTemperature, int _scheduleCycles, int _scheduleShape, int _maxRejectionsNumber, int _convergedSteps, double _convergedE); // run an SCMF Biased Monte Carlo search and returns the state with the best Energy
		void setLambda(double theLambda);
		double getLambda() const;
		vector<vector<double> > & getP();
		vector<vector<double> > & getE();

		double getStateP(vector<unsigned int> states) const; /*! get the probability of a rotameric state */

		double getPvariation() const;
		double getNumberOfCycles() const;
		void initialize();
		void randomize();
		
		void setMask(vector<vector<bool> > theMask);
		vector<vector<bool> > getMask() const;
		void readMaskFromFile(string _fileName);

		vector<unsigned int> getMostProbableState(); /*! return the most probable state */
		double getAverageEnergy();
		double getStateEnergy(vector<unsigned int> _state) const;

		vector<unsigned int> getRandomState(); /*! return a state selected randomly based on the probabilities */
		vector<unsigned int> moveRandomState(); /*! change randomly just one position of the random state */
		vector<unsigned int> moveRandomState(unsigned int _numOfSteps); /*! change randomly the specified number of position of the random state */
		void setCurrentState(vector<unsigned int> _state);
		
		void setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG);
		RandomNumberGenerator * getRandomNumberGenerator() const;

		void seed(unsigned int _seed); // use 0 for time based seed
		unsigned int getSeed() const; // get back the seed (useful to get the time based one)

		void setT(double _T);
		double getT() const;

		void setVerbose(bool _verbose = true);

				
	private:
		
		void initializeProbs();
		void initializeMask();
		int selectRandomStateAtPosition(int _position) const;
		void setup();
		void deletePointers();
		bool deleteRng;

		RandomNumberGenerator * pRng;

		double * pFixed;
		vector<vector<double> > * pSelfE;
		vector<vector<vector<vector<double> > > > * pPairE;
		vector<vector<double> > * pBaseLines;
		

		vector<vector<double> > selfConsE;
		vector<vector<double> > p;
		vector<vector<double> > pPrev;
		// the mask specified rotamers that should not be considered
		vector<vector<bool> > mask;

		vector<unsigned int> currentState; 
		
		int cycleCounter;
		double sumSquare;
		double lambda;
		double T;
		double RT;

		bool verbose;

		
};

inline void SelfConsistentMeanField::setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}

inline void SelfConsistentMeanField::setVerbose(bool _verbose) {
	verbose = _verbose;
}
inline RandomNumberGenerator * SelfConsistentMeanField::getRandomNumberGenerator() const {
	return pRng;
}

inline void SelfConsistentMeanField::seed(unsigned int _seed){
	if (_seed == 0) {
		pRng->setTimeBasedSeed();
	}
	else {
		pRng->setSeed(_seed);
	}
}


inline unsigned int SelfConsistentMeanField::getSeed() const {
	return pRng->getSeed();
}

}

#endif

