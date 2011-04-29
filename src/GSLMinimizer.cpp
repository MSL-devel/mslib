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

#include "GSLMinimizer.h"

using namespace MSL;

#include "MslTools.h"
#include "Timer.h"

GSLMinimizer::GSLMinimizer(){	
	setup(NULL,NULL);
}

GSLMinimizer::GSLMinimizer(EnergySet *_pEs){
	setup(_pEs,NULL);
}

GSLMinimizer::GSLMinimizer(EnergySet *_pEs,AtomPointerVector *_pAv){
	setup(_pEs,_pAv);
}

GSLMinimizer::~GSLMinimizer(){
}

void GSLMinimizer::setup(EnergySet *_pEset, AtomPointerVector *_pAv){
	pEset = _pEset;
	atoms = _pAv;
	
	stepsize = 0.1;
	tolerance = 0.01;
	maxIterations = 200;
	minimizeAlgorithm = STEEPEST_DESCENT;
}
void GSLMinimizer::Minimize(){

	if (atoms == NULL  || pEset == NULL){
		cerr << "ERROR GSLMinimizer::Minimize() either atoms or energySet is NULL.\n";
		exit(3);
	}


	// Variable declaration
	int retval,iter,status;
	double size;
	size = retval = iter = status = 0;
	int coordinateSize = atoms->size()*3;

	// Compute the initial value
	double initialValue, minimizedValue, deltaValue;
	minimizedValue = deltaValue = MslTools::doubleMax;
	initialValue   = pEset->calcEnergy();

	cout << "Setting up, initial value: "<<initialValue<<","<<coordinateSize<<endl;

	// Setup step size vector
	ss = NULL;
	ss = gsl_vector_alloc(coordinateSize);
	gsl_vector_set_all(ss,stepsize);  // Set all step sizes to .5


	// Create gsl version of our data..
	gslData = gsl_vector_alloc( coordinateSize);
	uint i = 0;
	for (uint a = 0; a < atoms->size();a++){
		gsl_vector_set(gslData,i++,(*atoms)(a).getCoor()[0]);
		gsl_vector_set(gslData,i++,(*atoms)(a).getCoor()[1]);
		gsl_vector_set(gslData,i++,(*atoms)(a).getCoor()[2]);
	}


	// Setup GSL Minimizer objects
	gsl_multimin_function f;
	f.f = &GSLMinimizer::my_static_f;        // the function itself
	f.n = coordinateSize;                 // ASSUMES cartieasan minimization
	f.params = (this);                     // MUST UPDATE THIS PARAMETER!!!

	gsl_multimin_function_fdf fdf;
	fdf.f      = &GSLMinimizer::my_static_f;      // the function itself
	fdf.df     = &GSLMinimizer::my_static_df;    // the gradient
	fdf.fdf    = &GSLMinimizer::my_static_fdf;  // combined 
	fdf.n      = coordinateSize;
	fdf.params = (this);                   // MUST UPDATE THIS PARAMETER!!!

	R = NULL;
	F = NULL;
	s1 = NULL;
	s2 = NULL;
	// set up minimizer
	switch (minimizeAlgorithm) {
	      case NELDERMEAD1:
		      R      = gsl_multimin_fminimizer_nmsimplex;
		      s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
		      retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
		      break;
//	      case NELDERMEAD2:
//		      R      = gsl_multimin_fminimizer_nmsimplex2;
//		      s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
//		      retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
//		      break;
//	      case NELDERMEAD2RAND:
//		      R      = gsl_multimin_fminimizer_nmsimplex2rand;
//		      s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
//		      retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
//		      break;
	      case CG_FLETCHER_REEVES:
		      F      = gsl_multimin_fdfminimizer_conjugate_fr;
		      s2     = gsl_multimin_fdfminimizer_alloc(F,coordinateSize); // initalize minimizer
		      retval = gsl_multimin_fdfminimizer_set(s2, &fdf, gslData, stepsize, tolerance);
		      break;
	      case CG_POLAK_RIBIERE:
		      F      = gsl_multimin_fdfminimizer_conjugate_pr;
		      s2     = gsl_multimin_fdfminimizer_alloc(F,coordinateSize); // initalize minimizer
		      retval = gsl_multimin_fdfminimizer_set(s2, &fdf, gslData, stepsize, tolerance);
		      break;
	      case BFGS:
		      F      = gsl_multimin_fdfminimizer_vector_bfgs2;
		      s2     = gsl_multimin_fdfminimizer_alloc(F,coordinateSize); // initalize minimizer
		      retval = gsl_multimin_fdfminimizer_set(s2, &fdf, gslData, stepsize, tolerance);
		      break;
	      case STEEPEST_DESCENT:
		      F      = gsl_multimin_fdfminimizer_steepest_descent;
		      s2     = gsl_multimin_fdfminimizer_alloc(F,coordinateSize); // initalize minimizer
		      retval = gsl_multimin_fdfminimizer_set(s2, &fdf, gslData, stepsize, tolerance);
		      break;

	      default:
		      cerr << "GSLMinimizer::Minimize() undefined minimization algorithm: "<<minimizeAlgorithm<<endl;
		      return;
	}


	if (s1 == NULL && s2 == NULL){
		cerr << "GSLMinimizer::Minimize() problem no minimizer object set up\n";
		return;
	}


	//  Store positions before minimization .... ????
	bool converged = false;
	bool derivateMinimization;
	if (s1 != NULL) {
		derivateMinimization = false;
	}
	if (s2 != NULL){
		derivateMinimization = true;
	}
	
	cout << "RUN MINIMIZATION"<<endl;
	do {
		iter++;

		// Perform a minimization step
		if (derivateMinimization){         status = gsl_multimin_fdfminimizer_iterate(s2); }
		else {                             status = gsl_multimin_fminimizer_iterate(s1);   }

		// Check for errors, handling should be better
		if (status) {

	 		printf("We have errors %d\n",status);
			break;  
		}

		// Get size info, used to see where we are in the minimization
		if (derivateMinimization) { 
			status = gsl_multimin_test_gradient (s2->gradient, tolerance);
		} else {
			size = gsl_multimin_fminimizer_size(s1);  
			status = gsl_multimin_test_size(size, tolerance);
		}

		// Check status to see if we have minimized (converged)
		if (status == GSL_SUCCESS) {  
			printf ("converged to minimum\n"); 
			converged = true;
		}           

		
	} while (status == GSL_CONTINUE && iter < maxIterations);


	// Print out values..
	minimizedValue  = 0.0;
	deltaValue     = 0.0;
	if (derivateMinimization){
		minimizedValue = gsl_multimin_fdfminimizer_minimum(s2);
		deltaValue     = (minimizedValue - initialValue);
		fprintf(stdout,"==> Minimized Value: %8.3f.  Value before minimization: %8.3f\n",minimizedValue,initialValue);
	} else {
		minimizedValue = gsl_multimin_fminimizer_minimum(s1);
		deltaValue     = (minimizedValue - initialValue);
		fprintf(stdout,"==> Minimized Value: %8.3f.  Value before minimization: %8.3f\n",minimizedValue,initialValue);
	}

	// I may have to retrieve lowest energy coordinates and set the atoms here...

	gsl_vector_free(ss);
	gsl_vector_free(gslData);

	// clean up and free memory
	if (derivateMinimization){
		gsl_multimin_fdfminimizer_free(s2);   
	} else {
		gsl_multimin_fminimizer_free(s1);      
	}

	R = NULL;
	F = NULL;
	
}


double GSLMinimizer::my_static_f(const gsl_vector *xvec_ptr, void *params){
	double temp;
	((GSLMinimizer *)params)->resetCoordinates(xvec_ptr);
	temp = ((GSLMinimizer *)params)->my_f(xvec_ptr,NULL);
	return temp;
}
void GSLMinimizer::my_static_df(const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr){
	// When commenting this out make sure my_df is defined
	((GSLMinimizer *)params)->resetCoordinates(xvec_ptr);
	((GSLMinimizer *)params)->my_df(xvec_ptr,params,df_ptr);
}
void GSLMinimizer::my_static_fdf(const gsl_vector *xvec_ptr, void *params,double *f_ptr, gsl_vector *df_ptr){
	// When commenting this out make sure my_fdf is defined
	((GSLMinimizer *)params)->resetCoordinates(xvec_ptr);
	((GSLMinimizer *)params)->my_fdf(xvec_ptr,params,f_ptr,df_ptr);
}



 double GSLMinimizer::my_f(const gsl_vector *xvec_ptr, void *params){

	double ans = 0.0;
	//cout << "HERE my_f "<<(*atoms).size()<<","<<(*atoms).size()*3<<endl;

	// extract new atom data and update coords
	//  if we stored gsl_vector in CartesianPoint as an alternative
	//  we could get rid of this loop...
	//cout << "F(X) = F(";

	

	// Call our funcion pointer..
	//Timer t;
	//double start = t.getWallTime();
	//ans = pEset->calcEnergy("active", "all");
	ans = pEset->calcEnergy();
	//fprintf(stdout," Ewall: %8.6f\n",t.getWallTime() - start);
	//ans = pEset->getAtomicPairwiseEnergy()->calculateTemplateEnergy(params->getSystem(),params->getPosition());
	

	double springE = 0.0;
	if ((*this).springControlledAtoms.size() > 0){
		springE = (*this).springEnergy.calcEnergy();
		ans += springE;
	}

	//fprintf(stdout, ") = %8.3f , %8.3f\n",ans,springE);
	
	/*
	if (ans < minValue) {
      
		minValue = ans;
		// save coor for each atom to "minCoor" or something like that..
	}
	*/

	return ans;
}

void GSLMinimizer::my_df(const gsl_vector *xvec_ptr, void *params, gsl_vector *df) {

	//cout << "HERE my_df "<<(*atoms).size()<<","<<(*atoms).size()*3<<endl;
	// Extract new atom data then compute gradient
	
	vector<double> gradient((*atoms).size()*3,0.0);
	//Timer t;
	//double start = t.getWallTime();
	pEset->calcEnergyGradient(gradient);
	//fprintf(stdout," Gwall: %8.6f\n",t.getWallTime() - start);


	// Do I have to set *df with gradient now? ... probably
	for (uint i=0; i < gradient.size();i+=1){
		gsl_vector_set(df, i, gradient[i]);
	}
}

void GSLMinimizer::my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
	//cout << "my_fdf"<<endl;

	// New way compute Energy and Gradient with single EnergySet Function
	vector<double> gradient((*atoms).size()*3,0.0);
	
	*f = pEset->calcEnergyAndEnergyGradient(gradient);

	for (uint i=0; i < gradient.size();i+=1){
		gsl_vector_set(df, i, gradient[i]);
	}


	// Old way, this works but is less efficient.
	/*
       *f = my_f(x, params); 
       my_df(x, params, df);
	*/
}





void GSLMinimizer::resetCoordinates(const gsl_vector *xvec_ptr){
	
	//cout << "Reset Coords "<<(*atoms).size()<<","<<(*atoms).size()*3<<endl;
	uint a = 0;
	for (uint i=0;i<atoms->size()*3;i+=3) {

		//*((*this).data)[i] = gsl_vector_get(xvec_ptr,i);
		(*atoms)(a).getCoor().setX(gsl_vector_get(xvec_ptr,i));
		(*atoms)(a).getCoor().setY(gsl_vector_get(xvec_ptr,i+1));
		(*atoms)(a).getCoor().setZ(gsl_vector_get(xvec_ptr,i+2));

		//fprintf(stdout, "%s",(*atoms)(a).toString().c_str());
		a++;
	}


}
