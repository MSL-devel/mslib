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

#ifndef MONTECARLOMANAGER_H
#define MONTECARLOMANAGER_H

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <vector>

// MSL includes
#include "MslExceptions.h"
#include "Real.h"
#include "release.h"

#include "MslTools.h"


using namespace std;
using namespace MSL;


/*! \brief MonteCarloManager object */

namespace MSL {



class MonteCarloManager {
	public:

                enum ANNEALTYPES {CONSTANT, LINEAR, EXPONENTIAL, SIGMOIDAL, SOFT};

		MonteCarloManager();
		MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber);
		MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber, int convergedSteps, double convergedE);
		~MonteCarloManager();

		void reset();
		void reset(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber);
		void reset(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber, int convergedSteps, double convergedE);

		void addStep(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber);

		bool accept(double energy);
		bool accept(double _energy, double _factor); // this is used for maitaining detailed balanced with SCMF/MC
		bool getComplete();
		void setEner(double ener);
		double getEner() const;
		int getCurrentCycle() const;
		int getTotalCycles() const;
		double getCurrentT() const;
		//double getStartT() const;
		//double getEndT() const;
		
		//int getShape() const;

		void setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG);
		RandomNumberGenerator * getRandomNumberGenerator() const;

		void seed(unsigned int _seed); // use 0 for time based seed
		unsigned int getSeed() const; // get back the seed (useful to get the time based one)
		//void seed(int _seed);
		//void getExternalRNG(RandomNumberGenerator * _pExternalRNG);

		double getRandomDouble() const {return -((double)rand() / (double)(RAND_MAX+1.0));};
		double getRandomDouble(double _max) const {return _max * -((double)rand() / (double)(RAND_MAX+1.0));};
		double getRandomDouble(double _max, double _min) const {return ((_max-_min) * -((double)rand() / (double)(RAND_MAX+1.0))) + _min;};

		int getRandomInt(int _max) const {return (int)(((double)(_max+1) * -((double)rand() / (double)(RAND_MAX+1.0))));};
		int getRandomInt(int _max, int _min) const {return (int)(((double)(_max+1-_min) * -((double)rand() / (double)(RAND_MAX+1.0)))) + _min;};

		void setConvergencyStop(int convergedSteps, double convergedE);
		
		string getReasonCompleted() const {return reasonCompleted;};
		
	private:
		void cool();
		void setup();
		void deletePointers();

		bool deleteRng;

		RandomNumberGenerator * pRng;
		
		double T;
		vector<double> startT;
		vector<double> endT;
		vector<int> shape;
		vector<int> totalCycles;
		vector<int> currCycle;
		vector<int> maxRejections;
		int currStep;
		double lastE;

		// quit if in for minDeltaInSteps acceptances the
		// energy decreases by less of delta 
		unsigned int minDeltaSteps;
		double minDeltaE;
		unsigned int deltaCounter;
		double deltaLast;

		double alpha;
		double beta;
		double gamma;
		int rejectionCounter;
		int totalAccepted;
		int totalRejected;
		string reasonCompleted;
};
inline void MonteCarloManager::setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}

inline RandomNumberGenerator * MonteCarloManager::getRandomNumberGenerator() const {
	return pRng;
}

inline void MonteCarloManager::seed(unsigned int _seed){
	if (_seed == 0) {
		pRng->setTimeBasedSeed();
	}
	else {
		pRng->setSeed(_seed);
	}
}


inline unsigned int MonteCarloManager::getSeed() const {
	return pRng->getSeed();
}

}

#endif
