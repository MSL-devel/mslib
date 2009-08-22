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

template<class T>GSLMinimizer<T>::GSLMinimizer(){
}

template<class T>GSLMinimizer<T>::~GSLMinimizer(){
}

template<class T>void GSLMinimizer<T>::Minimize(){

	// Variable declaration
	int retval,iter,status;
	double size;
	size = retval = iter = status = 0;

	// Compute the initial value
	double initialValue, minimizedValue, deltaValue;
	minimizedValue = deltaValue = doubleMax;
	initialValue   = (*this).CallFunc();

	cout << "Setting up, initial value: "<<initialValue<<endl;
	// Setup step size vector
	ss = gsl_vector_alloc((*this).data.size());
	gsl_vector_set_all(ss, (*this).stepsize);  // Set all step sizes to .5A


	// Create gsl version of our data..
	gslData = gsl_vector_alloc( (*this).data.size());
	for (uint i = 0; i < (*this).data.size();i++){
		gsl_vector_set(gslData,i,*((*this).data)[i]);
	}


	// Setup GSL Minimizer objects
	GSLMinimizer *this_ptr;
	this_ptr = (this);

	f.f = &GSLMinimizer::my_static_f;        // the function itself
	f.n = (*this).data.size();                 // ASSUMES cartieasan minimization
	f.params = this_ptr;                     // MUST UPDATE THIS PARAMETER!!!
    
	fdf.f      = &GSLMinimizer::my_static_f;      // the function itself
	fdf.df     = &GSLMinimizer::my_static_df;    // the gradient
	fdf.fdf    = &GSLMinimizer::my_static_fdf;  // combined 
	fdf.n      = (*this).data.size();                     // ASSUMES cartieasan minimization
	fdf.params = this_ptr;                   // MUST UPDATE THIS PARAMETER!!!




	// set up minimizer
	if ( (*this).minimizeAlgorithm == 1 ) {
		R = gsl_multimin_fminimizer_nmsimplex;
		s1        = gsl_multimin_fminimizer_alloc(R,(*this).data.size()); // initalize minimizer
		retval    = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);

	}
	else {    
		s2 = gsl_multimin_fdfminimizer_alloc(F,(*this).data.size()); // initalize minimizer
		retval      = gsl_multimin_fdfminimizer_set(s2, &fdf, gslData, (*this).stepsize, (*this).tolerance);
	}


	//  Store positions before minimization .... ????


	bool converged = false;
	do {
		iter++;

		// Perform a minimization step
		if ((*this).minimizeAlgorithm == 1) { status = gsl_multimin_fminimizer_iterate(s1);   }
		else {                              status = gsl_multimin_fdfminimizer_iterate(s2); }

		// Check for errors, handling should be better
		if (status) {

	 		printf("We have errors %d\n",status);
			cout << "Status is bad: "<<status<<endl;
			break;  
		}

		// Get size info, used to see where we are in the minimization
		if ((*this).minimizeAlgorithm == 1) { 
			size = gsl_multimin_fminimizer_size(s1);  
			status = gsl_multimin_test_size(size, (*this).tolerance);
		} else {
			status = gsl_multimin_test_gradient (s2->gradient, (*this).tolerance);
		}

		// Check status to see if we have minimized (converged)
		if (status == GSL_SUCCESS) {  
			printf ("converged to minimum at\n"); 
			converged = true;
		}           

		
	} while (status == GSL_CONTINUE && iter < (*this).maxIterations);


	// Print out values..
	minimizedValue  = 0.0;
	deltaValue     = 0.0;
	if ((*this).minimizeAlgorithm == 1) {
		minimizedValue = gsl_multimin_fminimizer_minimum(s1);
		deltaValue     = (minimizedValue - initialValue);
		cout << "==> Minimized Value: "<< deltaValue<<".  Value before minimization: " << initialValue<<".  Value after minimization: "<<minimizedValue<< endl;
	}


	//ss = gsl_vector_alloc(data.size());   // Initial vertex size vector for all cartiseans
	gsl_vector_free(ss);
	gsl_vector_free(gslData);

	// clean up and free memory
	if ((*this).minimizeAlgorithm == 1 ) {
		gsl_multimin_fminimizer_free(s1);      
	} else {
		gsl_multimin_fdfminimizer_free(s2);   
	}
	
}


template<class T>double GSLMinimizer<T>::my_static_f(const gsl_vector *xvec_ptr, void *params){
	double temp;
	temp = ((GSLMinimizer *)params)->my_f(xvec_ptr,NULL);
	return temp;
}
template<class T>void GSLMinimizer<T>::my_static_df(const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr){
	// When commenting this out make sure my_df is defined
	//((GSLMinimizer *)params)->my_df(xvec_ptr,params,df_ptr);
}
template<class T>void GSLMinimizer<T>::my_static_fdf(const gsl_vector *xvec_ptr, void *params_ptr,double *f_ptr, gsl_vector *df_ptr){
	// When commenting this out make sure my_fdf is defined
	//((GSLMinimizer *)params_ptr)->my_fdf(xvec_ptr,params_ptr,f_ptr,df_ptr);
}



template<class T> double GSLMinimizer<T>::my_f(const gsl_vector *xvec_ptr, void *params){

	double ans = 0.0;


	// extract new atom data and update coords
	//  if we stored gsl_vector in CartesianPoint as an alternative
	//  we could get rid of this loop...
	cout << "F(X) = F(";
	for (uint i=0;i<(*this).data.size();i++) {
		*((*this).data)[i] = gsl_vector_get(xvec_ptr,i);
		fprintf(stdout, "%8.3f, ",*(*this).data[i]);
	}
	

	// Call our funcion pointer..
	ans = (*this).CallFunc();

	double springE = 0.0;
	if ((*this).springControlledAtoms.size() > 0){
		springE = (*this).springEnergy.calcEnergy();
		ans += springE;
	}

	fprintf(stdout, ") = %8.3f , %8.3f\n",ans,springE);
	
	/*
	if (ans < minValue) {
      
		minValue = ans;
		// save coor for each atom to "minCoor" or something like that..
	}
	*/

	return ans;
}

class logTest;
template class GSLMinimizer<logTest>;


class EnergySet;
template class GSLMinimizer<EnergySet>;
