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

#include <string>
#include <stdlib.h>
#include <iostream>
#include "RandomNumberGenerator.h"
#include "MslTools.h"
#include "Timer.h"

using namespace std;

using namespace MSL;



int main(int argc, char **argv) {


	RandomNumberGenerator rng;

	if (argc > 1) {
		rng.setSeed(MslTools::toInt((string)argv[1]));
		cout << "Seeded with " << argv[1] << endl;
	} else {
		rng.setTimeBasedSeed();
	}

	// Print availble algorithms
	rng.printAvailableRNGAlgorithms();

	cout << "Get random double" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomDouble() << endl;;
	}
	cout << "========================" << endl;

	cout << "Get random double between 0.0 and 15.0" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomDouble(15.0) << endl;;
	}
	cout << "========================" << endl;

	cout << "Get random double between -7.0 and 8.0" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomDouble(-7.0, 8.0) << endl;;
	}

	cout << "Get random ints" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomInt() << endl;;
	}
	cout << "========================" << endl;

	cout << "Get random ints between 0 and 15" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomInt(15) << endl;;
	}
	cout << "========================" << endl;

	cout << "Get random ints between -7 and 8" << endl;
	for (unsigned int i=0; i<100; i++) {
		cout << rng.getRandomInt(-7, 8) << endl;;
	}


//	double p[] = { 5,1,7,100,1,1,10};
	vector<double> p(7, 0.0);
	p[0] = 5;
	p[1] = 1;
	p[2] = 7;
	p[3] = 100;
	p[4] = 1;
	p[5] = 1;
	p[6] = 10;
	double sum = 0.0;
	for (unsigned int i=0; i<p.size(); i++) {
		sum += p[i];
	}
	vector<double> relP;
	for (unsigned int i=0; i<p.size(); i++) {
		relP.push_back(p[i]/sum);
		cout << "* " << relP[i] << endl;
	}
	rng.setDiscreteProb(p);

	vector<int> counts;
	counts.resize(7);

	for (uint i = 0; i < 1000; i++){
		int index = rng.getRandomDiscreteIndex();

		//fprintf(stdout, "Random number from discrete prob: %d %8.3f\n",index,p[index]);
		counts[index]++;
		
	}

	fprintf(stdout, " Summary after selecting 1000 times from this probability distribution, does prob match counts?\n");
	for (uint i = 0; i < counts.size();i++){
	    fprintf(stdout, "Index %d , prob %6.4f, counts %-4d, count-based prob %8.3f\n",i,relP[i],counts[i],((double)counts[i])/10.0);
	}



}
