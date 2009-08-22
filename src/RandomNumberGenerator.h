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

/*
  Random Number Generation Object
 */

#ifndef RANDOM_H
#define RANDOM_H

#ifndef __GSL__
// This gives an error, but not the kind I was thinking...
#error message("GSL libraries not defined, can't use RandomNumberGenerator without them.")
#endif
// STL Includes
#include <iostream> 

// GSL Includes
#include <gsl/gsl_rng.h>

// MSL Includes
#include "MslTools.h"


class RandomNumberGenerator {
	
	public:
		RandomNumberGenerator();
		~RandomNumberGenerator();

		void setRNGType(string _type);
		string getRNGType();    // Stored in random 
		string getRNGTypeGSL(); // Directly from GSL

		void setRNGSeed(int _seed);
		int getRNGSeed();
		void setRNGTimeBasedSeed();

		double getRandomDouble();
		unsigned long int    getRandomInt();
		unsigned long int    getRandomIntLimit(int _upperLimit);

		void printAvailableRNGAlgorithms();
	private:		

		int randSeed;
		string randType;
#ifdef __GSL__
		const gsl_rng_type *Type;
		gsl_rng *rngObj;
#endif

};
#endif


