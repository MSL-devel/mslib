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

#include <string>
#include <stdlib.h>
#include <iostream>
#include "RandomNumberGenerator.h"

using namespace std;

using namespace MSL;



int main(int argc, char **argv) {


	RandomNumberGenerator rng;


	// Print availble algorithms
	rng.printAvailableRNGAlgorithms();



	double p[] = { 5,1,7,100,1,1,10};


	rng.setDiscreteProb(&p[0],7);

	vector<int> counts;
	counts.resize(7);

	for (uint i = 0; i < 1000; i++){
		int index = rng.getRandomDiscreteIndex();

		//fprintf(stdout, "Random number from discrete prob: %d %8.3f\n",index,p[index]);
		counts[index]++;
		
	}

	fprintf(stdout, " Summary after selecting 1000 times from this probability distribution, does prob match counts?\n");
	for (uint i = 0; i < counts.size();i++){
	    fprintf(stdout, "Index %d , prob %-3.0f, counts %-4d, count-based prob %8.3f\n",i,p[i],counts[i],((double)counts[i])/10.0);
	}



}
