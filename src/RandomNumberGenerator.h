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

/*
  Random Number Generation Object
 */

#ifndef RANDOM_H
#define RANDOM_H

// STL Includes
#include <iostream> 
#include <vector>

// GSL Includes
#ifdef __GSL__
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

// MSL Includes
#include "MslTools.h"


namespace MSL { 
class RandomNumberGenerator {
	
	public:
		RandomNumberGenerator(bool includeUpperLimit=true);
		~RandomNumberGenerator();

		int operator()(int _upperLimit); // get a random int

		void setRNGType(std::string _type);
		std::string getRNGType();    // Stored in random 
		std::string getRNGTypeGSL(); // Directly from GSL

		void setSeed(int _seed);
		void setTimeBasedSeed();
		unsigned int getSeed() const;

		// get random double (NOTE LOWER limits *INCLUDED*, UPPER limit *NOT* included!)
		double getRandomDouble(); // between 0 and 1 
		double getRandomDouble(double _upperLimit);  // between 0 and _upperLimit
		double getRandomDouble(double _lowerLimit, double _upperLimit); // between _lowerLimit and _upperLimit

		// get random int (NOTE LOWER and UPPER limits *INCLUDED*, unless includeUpperLimit set to false)
		unsigned long int    getRandomInt(); // between 0 and RAND_MAX
		unsigned long int    getRandomInt(unsigned long int _upperLimit); // between 0 and _upperLimit (included by default)
		long int             getRandomInt(long int _lowerLimit, long int _upperLimit); // between _lowerLimit and _upperLimit (both included by default)
		unsigned long int    getRandomIntLimit(int _upperLimit); // DEPRECATED

		std::vector <unsigned int> getRandomOrder (uint _size);			//returns vector between 0 and _size-1 numbered in a random order
		std::vector <unsigned int> getRandomOrder (uint _start, uint _end);


		/* The following takes a vector of probabilities and return a biased
		   ramdom indes. For example (0.25, 0.5, 0.25) is twice as likely to
		   return 1 than 0 or 2                                              */
		void setDiscreteProb(const std::vector<double> _prob);
		//void setDiscreteProb(const double *_prob,int _size); // old C style array function
		int getRandomDiscreteIndex();

		void printAvailableRNGAlgorithms();
	private:		

		int randSeed;
		int upperLimitOffset;
		std::string randType;
		

#ifndef __GSL__
		std::vector<double> cumulProb;
#else
		const gsl_rng_type *Type;
		gsl_rng *rngObj;
		gsl_ran_discrete_t *gsl_discrete;
#endif
};

inline int RandomNumberGenerator::operator()(int _upperLimit) { return getRandomInt(_upperLimit); }

inline unsigned int RandomNumberGenerator::getSeed() const {
	return randSeed;
}

}

#endif


