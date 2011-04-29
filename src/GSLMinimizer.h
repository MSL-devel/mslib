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

#ifndef GSLMIN_H
#define GSLMIN_H


#include "EnergySet.h"
#include "AtomPointerVector.h"
#include "AtomicPairwiseEnergy.h"

#ifndef __GSL__
#error message("GSLMinimizer can't compile unless GSL libraries are present")
#endif

#include <gsl/gsl_multimin.h>

using namespace std;

namespace MSL {
class GSLMinimizer  {

	public:

		GSLMinimizer();  
		GSLMinimizer(EnergySet *_pEs);  
		GSLMinimizer(EnergySet *_pEs,AtomPointerVector *_atoms);  
		~GSLMinimizer();

		// Minimize
		void Minimize();	

		// Get, Sets
		void setEnergySet(EnergySet *_pEs) { pEset = _pEs;}

		void setStepSize(double _stepsize);
		double getStepSize();
		
		void setTolerance(double _tol);
		double getTolerance();

		void setMaxIterations(int _maxIter);
		int getMaxIterations();

		void setMinimizeAlgorithm(int _minAlgo);
		int getMinimizeAlgorithm();

		void addAtoms(AtomPointerVector *_av);
		//AtomPointerVector& getAtoms();     


		void setSystem(System *_sys);
		System& getSystem();
		
		void setPosition(int _position);
		int getPosition();
		
		// Restrict Minimization when using atoms...
		//void freezeAtoms(AtomPointerVector &_av, double _springConstant);
		

		//void printData();
		enum MinimizingAlgorithms { NELDERMEAD1=1, NELDERMEAD2=2,NELDERMEAD2RAND=3,CG_FLETCHER_REEVES=4,CG_POLAK_RIBIERE=5,BFGS=6,STEEPEST_DESCENT=7};


	private:

		EnergySet* pEset;
		System* sys;
		int position;

		double stepsize;
		double tolerance;
		int    maxIterations;
		int    minimizeAlgorithm;    

		AtomPointerVector *atoms;

		AtomPointerVector springControlledAtoms;
		EnergySet springEnergy;

		//  Not using yet..
		static map<string,int> algorithmList;



		void setup(EnergySet *_pEset, AtomPointerVector *_pAv);
		void resetCoordinates(const gsl_vector *xvec_ptr);

		// Defining function pointers
		double  my_f   (const gsl_vector *xvec_ptr, void *params);
		void    my_df  (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);
		void    my_fdf (const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr);
		static double  my_static_f   (const gsl_vector *xvec_ptr, void *params);
		static void    my_static_df  (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);
		static void    my_static_fdf (const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr);

		// Multidimensional Minimizable Functions in GSL
		gsl_multimin_function f;
		gsl_multimin_function_fdf fdf;

		// Mulit-dimensional minmizer objects
		gsl_multimin_fminimizer     *s1;
		gsl_multimin_fdfminimizer   *s2;

		// GSL Constants for Minimization Algorithm Types
		const gsl_multimin_fminimizer_type    *R;
		const gsl_multimin_fdfminimizer_type  *F;

		// Double data in GSL
		gsl_vector *gslData;

		// Step size vector in GSL
		gsl_vector *ss;


	
};
#endif

//INLINES

inline void GSLMinimizer::setStepSize(double _stepsize){ stepsize = _stepsize; }
inline double GSLMinimizer::getStepSize(){ return stepsize;}


inline void GSLMinimizer::setTolerance(double _tol){ tolerance = _tol; }
inline double GSLMinimizer::getTolerance(){ return tolerance;}

inline void GSLMinimizer::setMaxIterations(int _maxIter){ maxIterations = _maxIter; }
inline int GSLMinimizer::getMaxIterations(){ return maxIterations;}


inline void GSLMinimizer::setMinimizeAlgorithm(int _minAlgo){ minimizeAlgorithm = _minAlgo; }
inline int GSLMinimizer::getMinimizeAlgorithm(){ return minimizeAlgorithm;}


inline void GSLMinimizer::addAtoms(AtomPointerVector *_av){ atoms = _av;} 


inline void GSLMinimizer::setSystem(System *_sys) { sys = _sys;}
inline System& GSLMinimizer::getSystem() { return (*sys); }

inline void GSLMinimizer::setPosition(int _pos) { position = _pos; }
inline int GSLMinimizer::getPosition() { return position; }


};
