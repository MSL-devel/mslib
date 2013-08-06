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

#include "GSLMinimizer.h"

using namespace MSL;

#include "MslTools.h"
#include "Timer.h"

GSLMinimizer::GSLMinimizer(){	
	setup(NULL,NULL);
}

GSLMinimizer::GSLMinimizer(System& _sys){
	setup(_sys.getEnergySet(),&(_sys.getAtomPointers()));
}

GSLMinimizer::GSLMinimizer(EnergySet* _es, AtomPointerVector* _av){
	setup(_es,_av);
}
	
GSLMinimizer::~GSLMinimizer(){
	removeConstraints();
}

void GSLMinimizer::setSystem(System& _sys) {
	setup(_sys.getEnergySet(),&(_sys.getAtomPointers()));
}

void GSLMinimizer::setup(EnergySet* _es, AtomPointerVector* _av){
	pEset = _es; 
	pAtoms = _av;
	// By default minimize all pAtoms
	if(pAtoms) {
		for(int i = 0; i < pAtoms->size(); i++) {
			// The minimization index has to start from 1
			(*pAtoms)[i]->setMinimizationIndex(i+1);
		}
	}
	
	stepsize = 0.1;
	tolerance = 0.01;
	maxIterations = 200;
	minimizeAlgorithm = BFGS;
}

void GSLMinimizer::resetConstraints() {

	// reset each interaction to zero energy
	for(map<Atom*,SpringConstraintInteraction*>::iterator it = springInteractionPointers.begin(); it != springInteractionPointers.end(); it++) {
		(it->second)->reset();
	}
}

void GSLMinimizer::fixAtoms(string _selection) {
	AtomSelection sel(*pAtoms);
	fixAtoms(sel.select(_selection));
}

void GSLMinimizer::fixAtoms(AtomPointerVector &_av) {
	//cout << "Number of pAtoms to fix " << _av.size() << endl;
	for(AtomPointerVector::iterator it = _av.begin(); it != _av.end(); it++) {
		if(*it) {
			(*it)->setMinimizationIndex(-1);
		}
	}
}


void GSLMinimizer::setConstraintForce(string _selection, double _springConstant) {
	AtomSelection sel(*pAtoms);
	setConstraintForce(sel.select(_selection),_springConstant);

}
void GSLMinimizer::setConstraintForce(AtomPointerVector &_av, double _springConstant) {
	
	//cout << "Number of pAtoms to constrain " << _av.size() << endl;
	if(pEset->getTotalNumberOfInteractions("SPRING_CONSTRAINT") == 0) {
		// no constraint interactions exist 
		for(AtomPointerVector::iterator it = _av.begin(); it != _av.end(); it++) {
			addSpringInteraction(*it,_springConstant);
		}
	} else {
		// if interactions already exist for these pAtoms just change the spring constant
		// else create new interactions
		
		for(AtomPointerVector::iterator it = _av.begin(); it != _av.end(); it++) {
			// check if interaction exists
			if(springInteractionPointers.find(*it) == springInteractionPointers.end()) {
				addSpringInteraction(*it,_springConstant);
			} else {
				SpringConstraintInteraction* scI = springInteractionPointers[*it];
				scI->setSpringConstant(_springConstant);
			}
		}
	}
	    
}

void GSLMinimizer::addSpringInteraction(Atom* _a1, double _springConstant) {
	SpringConstraintInteraction* scI = new SpringConstraintInteraction(*_a1,_springConstant,0.0);
	springInteractionPointers[_a1] = scI;
	pEset->addInteraction(scI);
}

void GSLMinimizer::removeConstraints() {
	pEset->eraseTerm("SPRING_CONSTRAINT");
	springInteractionPointers.clear();
}

bool GSLMinimizer::minimize(){

	if (pAtoms == NULL  || pEset == NULL){
		cerr << "ERROR GSLMinimizer::Minimize() either pAtoms or energySet is NULL.\n";
		return false;	 
	}

	// Variable declaration
	int retval,iter,status;
	double size;
	size = retval = iter = status = 0;
	int coordinateSize = pAtoms->size()*3;

	// Compute the initial value
	double initialValue, minimizedValue, deltaValue;
	minimizedValue = deltaValue = MslTools::doubleMax;
	initialValue   = pEset->calcEnergy();

	//cout << "Setting up, initial value: "<<initialValue<<","<<coordinateSize<<endl;

	// Double data in GSL
	gsl_vector *gslData = NULL;
	// Create gsl version of our data..
	gslData = gsl_vector_alloc( coordinateSize);
	uint i = 0;
	for (uint a = 0; a < pAtoms->size();a++){
		gsl_vector_set(gslData,i++,(*pAtoms)[a]->getX());
		gsl_vector_set(gslData,i++,(*pAtoms)[a]->getY());
		gsl_vector_set(gslData,i++,(*pAtoms)[a]->getZ());
	}


	// Step size vector in GSL
	gsl_vector *ss =  NULL;
	ss = gsl_vector_alloc(coordinateSize);
	gsl_vector_set_all(ss,stepsize);  // Set all step sizes to .5

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

	// GSL Constants for Minimization Algorithm Types
	const gsl_multimin_fminimizer_type    *R = NULL;
	const gsl_multimin_fdfminimizer_type  *F = NULL;

	// Mulit-dimensional minimizer objects
	gsl_multimin_fminimizer     *s1 = NULL;
	gsl_multimin_fdfminimizer   *s2 = NULL;

	// set up minimizer
	switch (minimizeAlgorithm) {
		case NELDERMEAD1:
			R      = gsl_multimin_fminimizer_nmsimplex;
			s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
#ifndef __GSL_OLD__
		// the following two functions compile only with GSL V.1.14 or above
		// set a $GSL_OLD = T environmental variable to avoid compiling them
		case NELDERMEAD2:
			R      = gsl_multimin_fminimizer_nmsimplex2;
			s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
		case NELDERMEAD2RAND:
			R      = gsl_multimin_fminimizer_nmsimplex2rand;
			s1     = gsl_multimin_fminimizer_alloc(R,coordinateSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
#endif
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
			return false;
	}


	if (s1 == NULL && s2 == NULL){
		cerr << "GSLMinimizer::Minimize() problem no minimizer object set up\n";
		return false;
	}


	//  Store positions before minimization .... ????
	//bool converged = false;
	bool derivateMinimization;
	if (s1 != NULL) {
		derivateMinimization = false;
	}
	if (s2 != NULL){
		derivateMinimization = true;
	}
	
	//cout << "RUN MINIMIZATION"<<endl;
	do {
		iter++;

		// Perform a minimization step
		if (derivateMinimization){    
			status = gsl_multimin_fdfminimizer_iterate(s2);
		} else {
			status = gsl_multimin_fminimizer_iterate(s1);   
		}

		// Check for errors, handling should be better
		if (status != GSL_SUCCESS && status != GSL_CONTINUE) {

	 		cerr << "Error 13425: " << gsl_strerror(status) << endl;
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
		//if (status == GSL_SUCCESS) {  
			//cout << "converged to minimum" << endl; 
			//converged = true;
		//}           

		
	} while (status == GSL_CONTINUE && iter < maxIterations);


	// Print out values..
	minimizedValue  = 0.0;
	deltaValue     = 0.0;
	if (derivateMinimization){
		minimizedValue = gsl_multimin_fdfminimizer_minimum(s2);
		deltaValue     = (minimizedValue - initialValue);
		//fprintf(stdout,"==> Minimized Value: %8.3f.  Value before minimization: %8.3f\n",minimizedValue,initialValue);
	} else {
		minimizedValue = gsl_multimin_fminimizer_minimum(s1);
		deltaValue     = (minimizedValue - initialValue);
		//fprintf(stdout,"==> Minimized Value: %8.3f.  Value before minimization: %8.3f\n",minimizedValue,initialValue);
	}
	gsl_vector_free(ss);
	gsl_vector_free(gslData);

	// clean up and free memory
	if (derivateMinimization){
		gsl_multimin_fdfminimizer_free(s2);   
	} else {
		gsl_multimin_fminimizer_free(s1);      
	}
	return true;
	
}


double GSLMinimizer::my_static_f(const gsl_vector *_xvec_ptr, void *_params){
	double temp;
	((GSLMinimizer *)_params)->resetCoordinates(_xvec_ptr);
	temp = ((GSLMinimizer *)_params)->my_f(_xvec_ptr,NULL);
	return temp;
}
void GSLMinimizer::my_static_df(const gsl_vector *_xvec_ptr, void *_params, gsl_vector *_df_ptr){
	// When commenting this out make sure my_df is defined
	((GSLMinimizer *)_params)->resetCoordinates(_xvec_ptr);
	((GSLMinimizer *)_params)->my_df(_xvec_ptr,_params,_df_ptr);
}
void GSLMinimizer::my_static_fdf(const gsl_vector *_xvec_ptr, void *_params,double *_f_ptr, gsl_vector *_df_ptr){
	// When commenting this out make sure my_fdf is defined
	((GSLMinimizer *)_params)->resetCoordinates(_xvec_ptr);
	((GSLMinimizer *)_params)->my_fdf(_xvec_ptr,_params,_f_ptr,_df_ptr);
}



double GSLMinimizer::my_f(const gsl_vector *_xvec_ptr, void *_params){
	// Calculate energy without the switching function.
	return pEset->calcEnergyWithoutSwitchingFunction();
}

void GSLMinimizer::my_df(const gsl_vector *_xvec_ptr, void *_params, gsl_vector *_df) {

	vector<double> gradient(pAtoms->size()*3,0.0);
	pEset->calcEnergyGradient(gradient);

	for (uint i=0; i < gradient.size();i+=1){
		gsl_vector_set(_df, i, gradient[i]);
	}
}

void GSLMinimizer::my_fdf(const gsl_vector *_x, void *_params, double *_f, gsl_vector *_df) {
	// New way compute Energy and Gradient with single EnergySet Function
	vector<double> gradient(pAtoms->size() * 3,0.0);
	
	*_f = pEset->calcEnergyAndEnergyGradient(gradient);
	//cout << "Energy: " << *_f << endl;

	for (uint i=0; i < gradient.size();i+=1){
		gsl_vector_set(_df, i, gradient[i]);
	}
	// Old way, this works but is less efficient.
	/*
       *f = my_f(x, params); 
       my_df(x, params, df);
	*/
}


void GSLMinimizer::resetCoordinates(const gsl_vector *_xvec_ptr){
	
	//cout << "Reset Coords "<<(*pAtoms).size()<<","<<(*pAtoms).size()*3<<endl;
	uint a = 0;
	for (uint i=0;i<pAtoms->size()*3;i+=3) {
		(*pAtoms)[a]->setCoor(gsl_vector_get(_xvec_ptr,i),gsl_vector_get(_xvec_ptr,i+1),gsl_vector_get(_xvec_ptr,i+2));
		a++;
	}

}
