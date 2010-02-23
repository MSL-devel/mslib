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
#include <math.h>

#include "CartesianGeometry.h"
#include "CartesianPoint.h"
#include "RandomNumberGenerator.h"

using namespace MSL;
using namespace std;


bool randomSample(vector<CartesianPoint *> _points, double _boxSize);
void printPartials(vector<double> &_derivatives);

RandomNumberGenerator rng;

int main(){
	rng.setRNGTimeBasedSeed();
	fprintf(stdout, "Running with random seed: %10d\n",rng.getRNGSeed());


	// Define 4 points
	CartesianPoint *p1 = new CartesianPoint(0.0,0.0,0.0);
	CartesianPoint *p2 = new CartesianPoint(0.0,0.0,0.0);
	CartesianPoint *p3 = new CartesianPoint(0.0,0.0,0.0);
	CartesianPoint *p4 = new CartesianPoint(0.0,0.0,0.0);

	// Store points in a list
	vector<CartesianPoint *> points;
	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);

	/*
	  Random positions within a box test... good for testing distance-related problems
	 */
	string testType = "RANDOM BOX TEST";
	fprintf(stdout, "********** %25s  **********\n\n",testType.c_str());
	int iterations = 5;
	int i = 0;
	while (i++ < iterations){

		// Get 4 random points in a 5 by 5 by 5 Angstrom box
		randomSample(points, pow(10.0,-5));


		vector<double> numericDistPartials     = CartesianGeometry::instance()->distanceNumericalDerivative(*p1,*p2);
		//vector<double> numericAnglePartials    = CartesianGeometry::instance()->angleNumericalDerivative(*p1,*p2,*p3,pow(10,-15));
		vector<double> numericDihedralPartials = CartesianGeometry::instance()->dihedralNumericalCosDerivative(*p1,*p2,*p3,*p4);


		// Compare here to analytic...
		vector<double> analyticDistPartials     = CartesianGeometry::instance()->distanceDerivative(*p1,*p2);
		//vector<double> analyticAnglePartials    = CartesianGeometry::instance()->angleDerivative(*p1,*p2,*p3);
		vector<double> analyticDihedralPartials = CartesianGeometry::instance()->dihedralCosDerivative(*p1,*p2,*p3,*p4);
		

		cout << "Numeric Distances: "<<endl;
		printPartials(numericDistPartials);
		
		cout << "Analytic Distances: "<<endl;
		printPartials(analyticDistPartials);

		cout << "Numeric Dihedrals: "<<endl;
		printPartials(numericDihedralPartials);
		
		cout << "Analytic Dihedrals: "<<endl;
		printPartials(analyticDihedralPartials);

		cout << "======="<<endl;


	}



	/*
	  3 points on line test 

	   A. close to 180 degrees
	 */

	map<int,double> errorMap;
	for (int d = -2; d < 20; d++){
		double deltaOffLine = exp(-d);

		testType = "3 co-linear near 180";
		fprintf(stdout, "********** %25s (%8.10f delta) **********\n\n",testType.c_str(),deltaOffLine);

		p1->setCoor(0,0,0);
		p2->setCoor(1,0,0);
		p3->setCoor(2,0,0);

	
		double aveError = 0.0;
		iterations = 5;
		i = 0;
		while (i++ < iterations){

			// Get random Y,Z values close to 0, biggest = 0.25
			p3->setY(rng.getRandomDouble()*deltaOffLine);
			p3->setZ(rng.getRandomDouble()*deltaOffLine);

			//vector<double> numericDistPartials     = CartesianGeometry::instance()->distanceNumericalDerivative(*p1,*p2);
			vector<double> numericAnglePartials    = CartesianGeometry::instance()->angleNumericalDerivative(*p1,*p2,*p3,0.01);
			//vector<double> numericDihedralPartials = CartesianGeometry::instance()->dihedralNumericalDerivative(*p1,*p2,*p3,*p4);


			// Compare here to analytic...
			//vector<double> analyticDistPartials     = CartesianGeometry::instance()->distanceDerivative(*p1,*p2);
			vector<double> analyticAnglePartials    = CartesianGeometry::instance()->angleDerivative(*p1,*p2,*p3);
			//vector<double> analyticDihedralPartials = CartesianGeometry::instance()->dihedralDerivative(*p1,*p2,*p3,*p4);
		

			cout << "Numeric Angles: "<<endl;
			printPartials(numericAnglePartials);
		
			cout << "Analytic Angles: "<<endl;
			printPartials(analyticAnglePartials);

			cout << "======="<<endl;



			for (uint j = 0; j < numericAnglePartials.size();j++){
				aveError += abs(numericAnglePartials[j] - analyticAnglePartials[j]);
			}

		}	
		aveError /= iterations;


		errorMap[d] = aveError;

		
	}

	/*
	  3 points on line test 

	   B. close to 0 degrees
	 */

	ofstream fout;
	fout.open("/tmp/angleDerError.txt");
	fout << "ERR180 ERR0 DELTA"<<endl;

	for (int d = -2; d < 20; d++){
		double deltaOffLine = exp(-d);
		testType = "3 co-linear near 0";
		fprintf(stdout, "********** %25s (%8.10f delta) **********\n\n",testType.c_str(),deltaOffLine);

		double aveError = 0.0;
		iterations = 5;
		i = 0;
		while (i++ < iterations){

			// Get random X value between 0 and 1
			p3->setX(0);

			// Get random Y,Z values close to 0, biggest = 0.25
			p3->setY(rng.getRandomDouble()*deltaOffLine);
			p3->setZ(rng.getRandomDouble()*deltaOffLine);

			//vector<double> numericDistPartials     = CartesianGeometry::instance()->distanceNumericalDerivative(*p1,*p2);
			vector<double> numericAnglePartials    = CartesianGeometry::instance()->angleNumericalDerivative(*p1,*p2,*p3,0.01);
			//vector<double> numericDihedralPartials = CartesianGeometry::instance()->dihedralNumericalDerivative(*p1,*p2,*p3,*p4);


			// Compare here to analytic...
			//vector<double> analyticDistPartials     = CartesianGeometry::instance()->distanceDerivative(*p1,*p2);
			vector<double> analyticAnglePartials    = CartesianGeometry::instance()->angleDerivative(*p1,*p2,*p3);
			//vector<double> analyticDihedralPartials = CartesianGeometry::instance()->dihedralDerivative(*p1,*p2,*p3,*p4);
		

			cout << "Numeric Angles: "<<endl;
			printPartials(numericAnglePartials);
		
			cout << "Analytic Angles: "<<endl;
			printPartials(analyticAnglePartials);

			cout << "======="<<endl;

			for (uint j = 0; j < numericAnglePartials.size();j++){
				aveError += abs(numericAnglePartials[j] - analyticAnglePartials[j]);
			}
		}	

		aveError /= iterations;

		char outstr[80];
		fout << errorMap[d]<<" "<<aveError<<" "<<deltaOffLine<<endl;
		//fout << sprintf(outstr, "%8.3f %8.3f %10.10f",errorMap[d],aveError,deltaOffLine)<<endl;
	}
	
	fout.close();



};



bool randomSample(vector<CartesianPoint *> _points, double _boxSize) {
	
	for (uint i = 0; i < _points.size();i++){
		_points[i]->setX(rng.getRandomDouble()*_boxSize);
		_points[i]->setY(rng.getRandomDouble()*_boxSize);
		_points[i]->setZ(rng.getRandomDouble()*_boxSize);
	}

	return true;
}

void printPartials(vector<double> &_derivatives){

	
	for (uint i = 0; i < _derivatives.size();i+=3){
		fprintf(stdout, "[ %8.3f %8.3f %8.3f ]\n",_derivatives[i],_derivatives[i+1],_derivatives[i+2]);
	}
}
