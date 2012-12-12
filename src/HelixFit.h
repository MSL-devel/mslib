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

#ifndef HELIXFIT_H
#define HELIXFIT_H

#ifndef __GSL__
#error message("GSLMinimizer can't compile unless GSL libraries are present")
#endif

#include <gsl/gsl_multimin.h>

// MSL Includes
#include "Matrix.h"
#include "CartesianPoint.h"
#include "Chain.h"
#include "AtomVector.h"


// Namespaces
using namespace std;

class HelixFit {
	public:

		HelixFit();
		HelixFit(Chain &ch);
		HelixFit(AtomVector &caOnly);
		~HelixFit();



		bool fit(int minimizationAlgorithm=NELDERMEAD1);

		void setChain(Chain &ch);
		void setAtoms(AtomVector &av);
		void setParameters(vector<double> _params);
		void setStepSize(double _stepSize);
		void setTolerance(double _tolerance);
		void setNumberSteps(int _numSteps);
		void setRigidMovementSet(int _moveSet);
		
		double getMinimizedValue() { return minimizedValue;}
		double getLastRMSD() { return lastFitRmsd; }
		double getLastSSQ() { return lastFitSSQ; }
		double getLastConstraint() { return lastFitConstraints; }

		vector<double> &getFittedParameters() { return parameters;}
		AtomVector& fittedHelix(vector<double> &_parameters);
		

		enum RigidMovements { NONE=-1,C2FIT=1,ALLFIT=2};
		enum MinimizingAlgorithms { NELDERMEAD1=1, NELDERMEAD2=2,NELDERMEAD2RAND=3,CG_FLETCHER_REEVES=4,CG_POLAK_RIBIERE=5,BFGS=6,STEEPEST_DESCENT=7};

		void setVerbose(bool _verboseOn=true);

	private:

		Chain *ch;
		AtomVector *caOnly;

		double stepSize;
		double tolerance;
		int numSteps;
		int moveSet;

		vector<double> parameters;
		vector<double> startingParameters;
		vector<int> angleParmeters;
		vector<pair<int,double> > constrainParameters; // pair = index of param, force constant

		bool verbose;

		// Transformation storage
		Matrix rotationMatrix;
		Matrix symmetryMatrix;
		CartesianPoint translationVector;

		// Last value minimized
		double minimizedValue;

		// After fittedHelix, these will be set.
		AtomVector lastFitHelix;		
		double lastFitRmsd;
		double lastFitSSQ;
		double lastFitConstraints;
		

		void setup();
		void createTransformation(vector<double> parameters);
		CartesianPoint createIdealHelixResidue(vector<double> _parameters, double index);

		// Defining function pointers
		double  RigidBodyHelixFit(const gsl_vector *xvec_ptr);
		void    RigidBodyHelixFitDerivatives(const gsl_vector *xvec_ptr, gsl_vector *df_ptr);

		static double  my_static_f   (const gsl_vector *xvec_ptr, void *params);
		static void    my_static_df  (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);
		static void    my_static_fdf (const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr);		

};
#endif
