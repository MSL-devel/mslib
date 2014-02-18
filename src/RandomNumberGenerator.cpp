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


#include "RandomNumberGenerator.h"
#include "Timer.h"

using namespace MSL;
using namespace std;



RandomNumberGenerator::RandomNumberGenerator(bool includeUpperLimit){

	// If random int includes upper limit, set to one, otherwise, to zero
	upperLimitOffset = 1;
	if (!includeUpperLimit) {
		upperLimitOffset = 0;
	}

	randSeed = 0; // time based is the default

#ifdef __GSL__
	// Using this , by default is set to use:
	//    GSL_RNG_TYPE="taus"
	//    GSL_RNG_SEED=123
	gsl_rng_env_setup(); 
	Type   = gsl_rng_default;
	rngObj = gsl_rng_alloc(Type);
	randType = "knuthran2";

	gsl_discrete = NULL;

	// default to knuthran2 (this will also seed it)
	setRNGType("knuthran2");

#else
	// seed it with time by default
	setTimeBasedSeed();
	randType = "";
#endif

}

RandomNumberGenerator::~RandomNumberGenerator(){
#ifdef __GSL__
	gsl_rng_free(rngObj);
	if (gsl_discrete != NULL){
		gsl_ran_discrete_free(gsl_discrete);
	}
	Type = NULL;
	rngObj = NULL;
#endif
}

void RandomNumberGenerator::setRNGType(string _type){

#ifndef __GSL__
	cerr << "Functionvoid RandomNumberGenerator::setRNGType(string _type) not available if not compiled with GLS" << endl; 
#else
	randType = _type;

	// Convert from string to gsl_rng_type; (Comments from GSL
	// documentation)

	// UNIX: the BSD 'rand' generator
	if (randType == "rand")   { Type = gsl_rng_rand;    }


	// MT19937: MT19937 generator of Makoto Matsumoto and Takuji Nishimura
	//    is a variant of the twisted generalized feedback shift-register algorithm, 
	//    and is known as the $(B!H(BMersenne Twister$(B!I(B generator.
	if (randType == "mt19937"){ Type = gsl_rng_mt19937; }

	// ranlx : second-generation version of the ranlux algorithm
	//    of LN|scher,which produces $(B!H(Bluxury random numbers$(B!I(B
	if (randType == "ranlxs0"){ Type = gsl_rng_ranlxs0; }
	if (randType == "ranlxs1"){ Type = gsl_rng_ranlxs1; }
	if (randType == "ranlxs2"){ Type = gsl_rng_ranlxs2; }
	if (randType == "ranlxd1"){ Type = gsl_rng_ranlxd1; }
	if (randType == "ranlxd2"){ Type = gsl_rng_ranlxd2; }

	// CMRG: combined multiple recursive generator by L'Ecuyer
	if (randType == "cmrg")   { Type = gsl_rng_cmrg;    }

	// Taus,Taus2: maximally equidistributed combined Tausworthe generator by L'Ecuyer
	if (randType == "taus")   { Type = gsl_rng_taus;    }
	if (randType == "taus2")  { Type = gsl_rng_taus2;   }

	// gfsr4: generator is like a lagged-fibonacci generator
	//   and produces each number as an xor'd sum of four previous values
	if (randType == "gfsr4")  { Type = gsl_rng_gfsr4;   } 
	
	
	// Knuth: second-order multiple recursive generator described
	// by Knuth in Seminumerical Algorithms, 3rd Ed., Section
	// 3.6. Knuth provides its C code. 
	if (randType == "knuthran")      { Type = gsl_rng_knuthran; }
	if (randType == "knuthran2002")  { Type = gsl_rng_knuthran2002; }
	if (randType == "knuthran2")     { Type = gsl_rng_knuthran2; }

	// Other algorithms from Kunth's Seminumerical Algortihms.
	if (randType == "borosh13")      { Type = gsl_rng_borosh13; }
	if (randType == "fishman18")     { Type = gsl_rng_fishman18; }
	if (randType == "fishman20")     { Type = gsl_rng_fishman20; }
	if (randType == "lecuyer21")     { Type = gsl_rng_lecuyer21; }
	if (randType == "waterman14")    { Type = gsl_rng_waterman14; }
	if (randType == "fishman2x")     { Type = gsl_rng_fishman2x; }
	if (randType == "coveyou")       { Type = gsl_rng_coveyou; }



	// Free previous randomNumber object and create new one with
	// this type.
	gsl_rng_free(rngObj);
        rngObj = gsl_rng_alloc(Type);
	setSeed(randSeed);
#endif

}
string RandomNumberGenerator::getRNGType(){
	return randType;
}
string RandomNumberGenerator::getRNGTypeGSL(){
#ifndef __GSL__
	cerr << "Functionvoid RandomNumberGenerator::setRNGType(string _type) not available if not compiled with GLS" << endl; 
	return "";
#else 
	return gsl_rng_name(rngObj);
#endif
}

void RandomNumberGenerator::setTimeBasedSeed(){
	/*
	Timer t;
	randSeed = (int)t.getWallTime();		
	gsl_rng_set(rngObj,randSeed);
	*/
	setSeed(0);
}

void RandomNumberGenerator::setSeed(int _seed){
	if (_seed == 0) {
		_seed = (unsigned int)time((time_t *)NULL);
	}
	randSeed = _seed;

#ifndef __GSL__
	srand(_seed);
#else
	gsl_rng_set(rngObj, _seed);
#endif
}

double RandomNumberGenerator::getRandomDouble(){
#ifndef __GSL__
	return (double)rand() / (double)RAND_MAX;
#else
	return gsl_rng_uniform(rngObj); 
#endif
}

double RandomNumberGenerator::getRandomDouble(double _upperLimit) {
#ifndef __GSL__
	return _upperLimit * (double)rand() / (double)RAND_MAX;
#else
	return getRandomDouble() * _upperLimit;
#endif
}
double RandomNumberGenerator::getRandomDouble(double _lowerLimit, double _upperLimit) {
#ifndef __GSL__
	return ((_upperLimit - _lowerLimit) * (double)rand() / (double)RAND_MAX) + _lowerLimit;
#else
	return getRandomDouble() * (_upperLimit - _lowerLimit) + _lowerLimit;
#endif
}

unsigned long int RandomNumberGenerator::getRandomInt(){
#ifndef __GSL__
	return rand();
#else
	return gsl_rng_get(rngObj); 
#endif
}


unsigned long int RandomNumberGenerator::getRandomInt(unsigned long int _upperLimit){
#ifndef __GSL__
	return rand() % (_upperLimit + upperLimitOffset);
#else
	return gsl_rng_uniform_int(rngObj,_upperLimit+upperLimitOffset); 
#endif
}

long int RandomNumberGenerator::getRandomInt(long int _lowerLimit, long int _upperLimit){
#ifndef __GSL__
	return _lowerLimit + rand() % (_upperLimit + upperLimitOffset - _lowerLimit);
#else
	return gsl_rng_uniform_int(rngObj,_upperLimit+upperLimitOffset-_lowerLimit) + _lowerLimit; 
#endif
}

unsigned long int RandomNumberGenerator::getRandomIntLimit(int _upperLimit){
	cerr << "WARNING: deprecated function unsigned long int RandomNumberGenerator::getRandomIntLimit(int _upperLimit)" << endl;
	return getRandomInt(_upperLimit);
}

std::vector <unsigned int> RandomNumberGenerator::getRandomOrder (uint _size) {
	unsigned int start = 0;
	unsigned int end = _size - 1;
	return getRandomOrder (start, end); 
}

std::vector <unsigned int> RandomNumberGenerator::getRandomOrder (uint _start, uint _end) {
	unsigned int size = (_end - _start) + 1;
	std::vector <unsigned int> ordered (size);
	std::vector <unsigned int> random (size);
	for (unsigned int i = _start; i <= _end; i++) {
		ordered[i] = i;
	}
	for (unsigned int j = 0; j<size; j++) {
		unsigned int randInt = getRandomDouble() * ordered.size();
		random[j] = ordered[randInt];
		ordered.erase(ordered.begin() + randInt);
	}
	return random;
}

void RandomNumberGenerator::printAvailableRNGAlgorithms(){

#ifndef __GSL__
	cout << "None (program compiled without GLS)" << endl;
#else
	cout << "RNG Algorithms:"<<endl;
	cout << "\trand"<<endl;
	cout << "\tmt19937"<<endl;
	cout << "\tranlxs0"<<endl;
	cout << "\tranlxs1"<<endl;
	cout << "\tranlxs2"<<endl;
	cout << "\tranlxd1"<<endl;
	cout << "\tranlxd2"<<endl;
	cout << "\tcmrg"<<endl;
	cout << "\ttaus"<<endl;
	cout << "\ttaus2"<<endl;
	cout << "\tgfsr4"<<endl;
	cout << "\tknuthran"<<endl;
	cout << "\tknuthran2002"<<endl;
	cout << "\tknuthran2"<<endl;
	cout << "\tborosh13"<<endl;
	cout << "\tfishman18"<<endl;
	cout << "\tfishman20"<<endl;
	cout << "\tfishman2x"<<endl;
	cout << "\tlecuyer21"<<endl;
	cout << "\twaterman14"<<endl;
	cout << "\tcoveyou"<<endl;
#endif

}


void RandomNumberGenerator::setDiscreteProb(const vector<double> _prob){
	if (_prob.size() == 0) {
		return;
	}
#ifndef __GSL__

	// create a comulative probability distrubution
	cumulProb.clear();
	double sum = 0.0;
	for (unsigned int i=0; i<_prob.size(); i++) {
		sum += _prob[i];
	}
	if (sum == 0.0) {
		return;
	}
	for (unsigned int i=0; i<_prob.size(); i++) {
		cumulProb.push_back(_prob[i]/sum);
		if (i != 0) {
			cumulProb[i] += cumulProb[i-1];
		}
	}
#else
        double *_probPtr = (double *)&_prob[0];

	// make sure we dealloc the space
	if(gsl_discrete) {
		gsl_ran_discrete_free(gsl_discrete);
	}
	gsl_discrete = gsl_ran_discrete_preproc(_prob.size(), _probPtr);
#endif
}

/*
void RandomNumberGenerator::setDiscreteProb(const double *_prob, int _size){

	gsl_discrete = gsl_ran_discrete_preproc(_size,_prob);

}
*/


int RandomNumberGenerator::getRandomDiscreteIndex(){

#ifndef __GSL__
	double n = getRandomDouble(); // get a number 0-1
	
	for (unsigned int i=0; i<cumulProb.size(); i++) {
		if (n < cumulProb[i]) {
			return i;
		}
	}
	return cumulProb.size() - 1;
		
#else
	if (gsl_discrete == NULL || rngObj == NULL){
		return -1;
	}

	return (gsl_ran_discrete(rngObj,gsl_discrete));
#endif
}
