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

#ifndef GSLMIN_H
#define GSLMIN_H


#include "EnergySet.h"
#include "AtomPointerVector.h"
#include "AtomSelection.h"
#include "CharmmEnergyCalculator.h"
#include "SpringConstraintInteraction.h"

#ifndef __GSL__
#error message("GSLMinimizer can't compile unless GSL libraries are present")
#endif

#include <gsl/gsl_multimin.h>

using namespace std;

namespace MSL {
class GSLMinimizer  {

	public:

		GSLMinimizer();  
		GSLMinimizer(System& _sys);  
		GSLMinimizer(EnergySet* _es, AtomPointerVector* _av);
		~GSLMinimizer();

		// Minimize
		//TODO: LOOK INTO CREATING ENERGY SUBSETS
		bool minimize();	

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

		/* Restrict Minimization when using pAtoms...
		* Adds CONSTRAINT interactions to the system's energySet to simulate pAtoms on a spring
		* So once minimization is done, a single call to removeConstraints is necessary
		*/
		void setConstraintForce(string _selection, double _springConstant); // overwrites springConstant of existing interactions if selected
		void setConstraintForce(AtomPointerVector &_av, double _springConstant); // overwrites springConstant of existing interactions if selected

		// exclude pAtoms from minimization - sets their minimization index to -1
		void fixAtoms(AtomPointerVector &_av);
		void fixAtoms(string _selection);

		void resetConstraints(); // resets all constraint terms to 0.0
		void removeConstraints(); // deletes all constraint energy terms
		

		//void printData();
		enum MinimizingAlgorithms { NELDERMEAD1=1, NELDERMEAD2=2,NELDERMEAD2RAND=3,CG_FLETCHER_REEVES=4,CG_POLAK_RIBIERE=5,BFGS=6,STEEPEST_DESCENT=7};


	private:

		AtomPointerVector* pAtoms;
		EnergySet* pEset;
		map<Atom*,SpringConstraintInteraction*> springInteractionPointers;

		double stepsize;
		double tolerance;
		int    maxIterations;
		int    minimizeAlgorithm;    

		void addSpringInteraction(Atom* a1, double _springConstant);

		//  Not using yet..
		//static map<string,int> algorithmList;

		void setup(EnergySet* _es, AtomPointerVector* _av);
		void resetCoordinates(const gsl_vector *xvec_ptr);

		// Defining function pointers
		double  my_f   (const gsl_vector *xvec_ptr, void *params);
		void    my_df  (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);
		void    my_fdf (const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr);
		static double  my_static_f   (const gsl_vector *xvec_ptr, void *params);
		static void    my_static_df  (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);
		static void    my_static_fdf (const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr);

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

};
