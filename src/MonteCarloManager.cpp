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

#include "MonteCarloManager.h"

using namespace std;
using namespace MSL;


MonteCarloManager::MonteCarloManager() {
	reset(298, 298, 0, CONSTANT, 0, 0, 0.0);
	setup();
}

MonteCarloManager::MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber) {
	// default with no convergency stop
	reset(startingTemperature, endingTemperature, scheduleCycles, scheduleShape, maxRejectionsNumber, 0, 0.0);
	setup();
}

MonteCarloManager::MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber, int convergedSteps, double convergedE) {
	reset(startingTemperature, endingTemperature, scheduleCycles, scheduleShape, maxRejectionsNumber, convergedSteps, convergedE);
	setup();
}

MonteCarloManager::~MonteCarloManager() {
	deletePointers();
}

void MonteCarloManager::setup() {
	deleteRng = true;
	pRng = new RandomNumberGenerator;
}

void MonteCarloManager::deletePointers() {
	if (deleteRng == true) {
		delete pRng;
	}
}

void MonteCarloManager::reset() {
	// reset counter and temperature to the original setting
	T = 0; 
	lastE = +1e+100;
	currStep = 0;
	rejectionCounter = 0;
	totalAccepted = 0;
	totalRejected = 0;

	minDeltaSteps = 0;
	minDeltaE = 0.0;
	deltaCounter = 0;
	deltaLast = lastE;

	startT.clear();
	endT.clear();
	shape.clear();
	totalCycles.clear();
	currCycle.clear();
	maxRejections.clear();

	reasonCompleted = "";
}

void MonteCarloManager::reset(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber) {
	// default with no convergency stop
	reset(startingTemperature, endingTemperature, scheduleCycles, scheduleShape, maxRejectionsNumber, 0, 0.0);
}

void MonteCarloManager::reset(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber, int convergedSteps, double convergedE) {
	reset();
	T = startingTemperature; 
	minDeltaSteps = convergedSteps;
	minDeltaE = convergedE;
	addStep(startingTemperature, endingTemperature, scheduleCycles, scheduleShape, maxRejectionsNumber);
}

void MonteCarloManager::addStep(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber) {

	if (scheduleShape < 0 || scheduleShape > 4) {
		cerr << "ERROR 4656: unknown value " << scheduleShape << " for montecarlo curve shape (CONSTANT = 0, LINEAR = 1, EXPONENTIAL1 = 2, SIGMOIDAL = 3, SOFT = 4) in MonteCarloManager::MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape)" << endl;
		exit(4656);
	}
	if (endingTemperature < 0.0) {
		cerr << "ERROR 4657: ending temperature " << endingTemperature << "<0.0 in MonteCarloManager::MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape)" << endl;
		exit(4657);
	}
	if (startingTemperature < endingTemperature) {
		cerr << "ERROR 4658: starting temperature " << startingTemperature << "< of ending temperature " << endingTemperature << " in MonteCarloManager::MonteCarloManager(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape)" << endl;
		exit(4658);
	}
	/*
	if (scheduleShape == EXPONENTIAL && endingTemperature == 0.0) {
		cerr << "ERROR 4659: ending temperature of 0 with exponential (2) cooling in void MonteCarloManager::reset(double startingTemperature, double endingTemperature, int scheduleCycles, int scheduleShape, int maxRejectionsNumber)" << endl;
		exit(4659);
	}
	*/

	startT.push_back(startingTemperature);
	endT.push_back(endingTemperature);
	shape.push_back(scheduleShape);
	totalCycles.push_back(scheduleCycles);
	currCycle.push_back(0);
	maxRejections.push_back(maxRejectionsNumber);

	reasonCompleted = "";

}

/*
void MonteCarloManager::getExternalRNG(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}
*/

bool MonteCarloManager::accept(double _energy) {
	return accept(_energy, 1.0);
}

bool MonteCarloManager::accept(double _energy, double _factor) {

	/***********************************
	 *  Reject if we are done with all steps
	 ***********************************/
	if (currStep >= currCycle.size()) {
		return false;
	}


	bool accept_flag = false;

	/***********************************
	 *  Metropolis criterion
	 ***********************************/
	double exp = _factor * pow(M_E,(lastE - _energy) / (MslTools::R * T));
	if (exp > 1) {
		lastE = _energy;
		accept_flag = true;
	} else {

		double randomN = pRng->getRandomDouble();

		if (exp > randomN) {
			lastE = _energy;
			accept_flag = true;
		}
	}

	/*
	if (_energy < lastE) {
		lastE = _energy;
		accept_flag = true;
	} else {
		double diff = _energy - lastE;
		double randomN = -(double)rand()/(double)(RAND_MAX+1.0);

		double exp = pow(M_E,(-(diff)) / (R * T));
		//cout << "Exp: " << exp << "  RAND: " << randomN << endl;
		if (exp > randomN) {
			lastE = _energy;
			accept_flag = true;
		}
	}
	*/

	if (!accept_flag) {
		rejectionCounter++;
		totalRejected++;
	} else {
		totalAccepted++;
	}
	if (accept_flag || rejectionCounter > maxRejections[currStep]) {
		/***************************************
		 *  If the move is accepted cool and
		 *  reset the rejection counter
		 ***************************************/
		currCycle[currStep]++;
		rejectionCounter = 0;
		if (currCycle[currStep] < totalCycles[currStep]) {
			cool();
		} else {
			currStep++;
			if (currStep < currCycle.size()) {
				T = startT[currStep];
			}
		}
		if (deltaLast - _energy > minDeltaE) {
			deltaCounter = 0;
			deltaLast = _energy;
		} else {
			deltaCounter++;
		}
	}
	return accept_flag;
}

bool MonteCarloManager::getComplete() {
	if (currCycle.size() == 0) {
		reasonCompleted = "No cycles set";
		return true;
	}
	if (currStep >= currCycle.size()) {
		reasonCompleted = "Completed all cycles. Accepted " + MslTools::intToString(totalAccepted) + "/" + MslTools::intToString(totalAccepted+totalRejected);
		return true;
	}
	if (currStep == currCycle.size() - 1 && currCycle[currStep] >= totalCycles[currStep]) {
		reasonCompleted = "Completed all cycles. Accepted " + MslTools::intToString(totalAccepted) + "/" + MslTools::intToString(totalAccepted+totalRejected);
		return true;
	}
	if (minDeltaSteps != 0 && minDeltaE != 0.0) {
		// stop if the energy has not lowered of the 
		// required delta over the required number of
		// steps
		if (deltaCounter > minDeltaSteps) {
			reasonCompleted = "Reached " + MslTools::intToString(deltaCounter) + " cycles without improvement of " + MslTools::doubleToString(minDeltaE) + ". Accepted " + MslTools::intToString(totalAccepted) + "/" + MslTools::intToString(totalAccepted+totalRejected);
			return true;
		}
	}
	return false;
}

void MonteCarloManager::setConvergencyStop(int convergedSteps, double convergedE) {
	minDeltaSteps =  convergedSteps;
	minDeltaE = convergedE;
}

void MonteCarloManager::setEner(double ener) {
	lastE = ener;
}

double MonteCarloManager::getEner() const {
	return lastE;
}

int MonteCarloManager::getCurrentCycle() const {
	// check this function
	int sum = 0;
	for (int i=0; i< currCycle.size(); i++) {
		sum += currCycle[currStep];
	}
	return sum;
}

int MonteCarloManager::getTotalCycles() const {
	int sum = 0;
	for (int i=0; i< totalCycles.size(); i++) {
		sum += totalCycles[currStep];
	}
	return sum;
}

double MonteCarloManager::getCurrentT() const {
	return T;
}

/*
double MonteCarloManager::getStartT() const {
	return startT;
}

double MonteCarloManager::getEndT() const {
	return endT;
}

int MonteCarloManager::getShape() const {
	return shape;
}
*/


void MonteCarloManager::cool(){

	switch (shape[currStep]) {
		case CONSTANT:
			/*************************************
			 *  #  0
			 *************************************/
			break;
		case LINEAR:
			/*************************************
			 *  #  1
			 *************************************/
			T = (startT[currStep] - (currCycle[currStep] * ((startT[currStep] - endT[currStep]) / (double)totalCycles[currStep])));
			break;
		case EXPONENTIAL:
			/*********************************************
			 *  #  2
			 * exponential decay constant for the formula
			 * of the temperature of the i-th step
			 *
			 *                          startT-endT
			 * Ti = (startT-endT) * (-----------------)^(i/totalCycles) + endT
			 *                       500*(startT-endT)
			 *********************************************/
			/*
			alpha = (double)(startT[currStep] - endT[currStep]);
			beta = alpha/500.0;
			T = alpha * pow((beta/alpha), ((double)currCycle[currStep]/(double)totalCycles[currStep])) + endT[currStep];
			*/
			
			/***********************************************************
			*   Ti = startT * 10 ^ (-i * currCycle)
			*   where i = - log (endT/startT) / totalCycles
			*
			************************************************************/
			T = startT[currStep]  * pow(10,currCycle[currStep] * log10(endT[currStep]/startT[currStep])/totalCycles[currStep]);
			break;
		case SIGMOIDAL:
			/*********************************************
			 *  #  3
			 * sigmoidal curve, flat start, steep middle
			 * discend, flat end
			 *
			 *                     startT - endT
			 * Ti = ------------------------------------------ + endT
			 *      1 + e^(12/totalCycles * (i - (totalCycles/2)))

			 * The constant 12 controls the steepness in the middle region
			 * - higher the value steeper the slope
			 *********************************************/
			/*
			alpha = (double)(startT[currStep] - endT[currStep]);
			beta = (double)totalCycles[currStep] / 2.0;
			gamma = 12.0/(double)totalCycles[currStep];
			T = endT[currStep] + (alpha / (1.0 + exp(gamma * ((double)currCycle[currStep] - beta ))));
			*/
			
			/*********************************************
			 * sigmoidal curve, flat start, steep middle
			 * discend, flat end
			 *
			 *      (startT - endT)    
			 * Ti = --------------------  * (R(i) - R(endT)) + endT 
			 *      (R(startT)- R(endT))           
			 *	
			 *	           		1
                         * R(i) =   -----------------------------------------------
                         *        (1 + e^(12/totalCycles * (i - (totalCycles/2))))
                         *
			 * The constant 12 controls the steepness in the middle region
			 * - higher the value steeper the slope
			 *********************************************/
			alpha = startT[currStep] - endT[currStep];
			beta = 1.0 / (1.0 + exp(6.0)); // R(endT)
			gamma = 1.0 / (1.0 + exp(-6.0)) - beta; // R(startT) - R(endT)
			T = (alpha/gamma) * ( 1.0/ (1.0 + exp(12.0/totalCycles[currStep] * (currCycle[currStep]- (totalCycles[currStep]/2.0)))) - beta)	+ endT[currStep]; 
			break;
		case SOFT:
			/*********************************************
			 *  #  4
			 *  Soft initial decay, steeper midrange cooling
			 *  and slow convergency
			 *      startT - endT             12 * currCycle
			 * Ti = ------------- * (1 + cosh(------------- - 5) + endT
			 *            2                        N
			 *********************************************/
			alpha = 0.5 * (double)(startT[currStep]- endT[currStep]);
			T = (((double)(startT[currStep] - endT[currStep])) / cosh( (((double)currCycle[currStep] * 10.0)/ (double)totalCycles[currStep]))) + endT[currStep];
			break;
	}


}

/*
void MonteCarloManager::seed(int _seed) {
	if (_seed == 0)  {
		pRng->setRNGTimeBasedSeed();
	}
	else {
		pRng->setRNGSeed(_seed);
	}
}
*/


