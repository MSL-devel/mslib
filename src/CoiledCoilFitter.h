/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
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

#ifndef COILEDCOILFITTER_H
#define COILEDCOILFITTER_H


#ifndef __GSL__
#error message("CoiledCoilFitter can't compile unless GSL libraries are present")
#endif

#include "AtomContainer.h"
#include <gsl/gsl_multimin.h>

using namespace std;

namespace MSL {
class CoiledCoilFitter  {

	public:

                CoiledCoilFitter();  
		~CoiledCoilFitter();

		/*
		  Must add helices in order. A,B,C,D means A neighbors D and B.  B neighbors A and C. C neighbors B and D.  
		  // Add N,Ca,C,Cb only
		 */
		void addNextHelix(AtomPointerVector *_av);

		// Minimize
		//TODO: LOOK INTO CREATING ENERGY SUBSETS
		bool fit();

		// Get, Sets
		void setSystem(System& _sys);

		void setStepSize(double _stepsize);
		double getStepSize();
		
		void setTolerance(double _tol);
		double getTolerance();

		void setMaxIterations(int _maxIter);
		int getMaxIterations();

		void setMinimizeAlgorithm(int _minAlgo);
		int getMinimizeAlgorithm();

		void setSymmetry(string _sym);
		string getSymmetry();

		vector<double> getMinimizedParameters();

		//void printData();
		enum MinimizingAlgorithms { NELDERMEAD1=1, NELDERMEAD2=2,NELDERMEAD2RAND=3};


	private:

		AtomContainer pAtoms;

		double stepsize;
		double tolerance;
		int    maxIterations;
		int    minimizeAlgorithm;    
		string symmetry_;
		int numHelices_;
		int numResidues_;
		int parameterSize;
		int minIndex;
		vector<double> params;

		// Defining function pointers
		double  my_f   (const gsl_vector *xvec_ptr, void *params);
		static double  my_static_f   (const gsl_vector *xvec_ptr, void *params);


};
#endif

//INLINES

inline void CoiledCoilFitter::setStepSize(double _stepsize){ stepsize = _stepsize; }
inline double CoiledCoilFitter::getStepSize(){ return stepsize;}


inline void CoiledCoilFitter::setTolerance(double _tol){ tolerance = _tol; }
inline double CoiledCoilFitter::getTolerance(){ return tolerance;}

inline void CoiledCoilFitter::setMaxIterations(int _maxIter){ maxIterations = _maxIter; }
inline int CoiledCoilFitter::getMaxIterations(){ return maxIterations;}


inline void CoiledCoilFitter::setMinimizeAlgorithm(int _minAlgo){ minimizeAlgorithm = _minAlgo; }
inline int CoiledCoilFitter::getMinimizeAlgorithm(){ return minimizeAlgorithm;}

inline void CoiledCoilFitter::setSymmetry(string _sym) { symmetry_ = _sym; }
inline string  CoiledCoilFitter::getSymmetry() { return symmetry_;}

 inline vector<double> CoiledCoilFitter::getMinimizedParameters() { return params; }
};
