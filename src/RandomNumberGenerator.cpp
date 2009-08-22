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


#include "RandomNumberGenerator.h"
#include "Timer.h"


RandomNumberGenerator::RandomNumberGenerator(){

	// Using this , by default is set to use:
	//    GSL_RNG_TYPE="taus"
	//    GSL_RNG_SEED=123
	gsl_rng_env_setup(); 
	Type   = gsl_rng_default;
	rngObj = gsl_rng_alloc(Type);
	randSeed = gsl_rng_default_seed;
	randType = gsl_rng_name(rngObj);
}
RandomNumberGenerator::~RandomNumberGenerator(){
	gsl_rng_free(rngObj);
	Type = NULL;
	rngObj = NULL;
}

void RandomNumberGenerator::setRNGType(string _type){

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

}
string RandomNumberGenerator::getRNGType(){
 
	return randType;
}
string RandomNumberGenerator::getRNGTypeGSL(){
	return gsl_rng_name(rngObj);
}

void RandomNumberGenerator::setRNGTimeBasedSeed(){
	Timer t;
	randSeed = (int)t.getWallTime();		
	gsl_rng_set(rngObj,randSeed);
}

void RandomNumberGenerator::setRNGSeed(int _seed){
	randSeed = _seed;
	gsl_rng_set(rngObj, randSeed);
}
int RandomNumberGenerator::getRNGSeed(){
	return randSeed;
}

double RandomNumberGenerator::getRandomDouble(){
	return gsl_rng_uniform(rngObj); 
}


unsigned long int RandomNumberGenerator::getRandomInt(){
	return gsl_rng_get(rngObj); 
}


unsigned long int RandomNumberGenerator::getRandomIntLimit(int _upperLimit){
	return gsl_rng_uniform_int(rngObj,_upperLimit); 
}


void RandomNumberGenerator::printAvailableRNGAlgorithms(){

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

}

