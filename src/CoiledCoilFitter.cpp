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


#include "CoiledCoilFitter.h"

using namespace MSL;

#include "MslTools.h"
#include "CoiledCoils.h"
#include "Timer.h"

CoiledCoilFitter::CoiledCoilFitter(){	

  numHelices_ = 0;
  numResidues_ = 0;
  stepsize = 0.1;
  tolerance = 0.01;
  maxIterations = 200;
  symmetry_ = "C";
  parameterSize = 3;
  minimizeAlgorithm = NELDERMEAD1;
  minIndex = 1;
}

	
CoiledCoilFitter::~CoiledCoilFitter(){
}

void CoiledCoilFitter::addNextHelix(AtomPointerVector *_helix){
     pAtoms.addAtoms(*_helix);

     if (numHelices_ == 0){
       numResidues_ = _helix->size() / 4;
     } else {
       if (_helix->size() / 4 != numResidues_){
	 cerr << "ERROR 23423 CoiledCoilFitter::addNextHelix, last helix had "<<numResidues_<<" residues and this one has "<<_helix->size()/4<<" each helix needs the same number of residues"<<endl;
	 exit(23423);
       }
     }
     numHelices_++;
}

bool CoiledCoilFitter::fit(){

      if (pAtoms.size() == 0){
		cerr << "ERROR CoiledCoilFitter::fit() no atoms to minimize.\n";
		return false;	 
	}

	// Variable declaration
	int retval,iter,status;
	double size;
	size = retval = iter = status = 0;
	int parameterSize = 3;

	// Compute the initial value
	double initialValue, minimizedValue, deltaValue;
	minimizedValue = deltaValue = initialValue = MslTools::doubleMax;


	//cout << "Setting up, initial value: "<<initialValue<<","<<coordinateSize<<endl;

	// Double data in GSL
	gsl_vector *gslData = NULL;

	// Create gsl version of our data..
	gslData = gsl_vector_alloc( parameterSize );
	gsl_vector_set(gslData,0,5.0);                           // superhelical radius = r0 = 5
	gsl_vector_set(gslData,1,(2*M_PI*5)/tan(M_PI*12/180)); // super-helical pitch (12 degrees converted into Angstroms = distance between superhelical repeats)
	gsl_vector_set(gslData,2,0.0);                           // phi1 = alpha helical phase = 0
	


	// Step size vector in GSL
	gsl_vector *ss =  NULL;
	ss = gsl_vector_alloc(parameterSize);
	gsl_vector_set_all(ss,stepsize);  // Set all step sizes to .5

	// Setup GSL Minimizer objects
	gsl_multimin_function f;
	f.f = &CoiledCoilFitter::my_static_f;   // the function itself
	f.n = parameterSize;                     // ASSUMES cartieasan minimization
	f.params = (this);                     // MUST UPDATE THIS PARAMETER!!!

	// GSL Constants for Minimization Algorithm Types
	const gsl_multimin_fminimizer_type    *R = NULL;

	// Mulit-dimensional minimizer objects
	gsl_multimin_fminimizer     *s1 = NULL;

	// set up minimizer
	switch (minimizeAlgorithm) {
		case NELDERMEAD1:
			R      = gsl_multimin_fminimizer_nmsimplex;
			s1     = gsl_multimin_fminimizer_alloc(R,parameterSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
#ifndef __GSL_OLD__
		// the following two functions compile only with GSL V.1.14 or above
		// set a $GSL_OLD = T environmental variable to avoid compiling them
		case NELDERMEAD2:
			R      = gsl_multimin_fminimizer_nmsimplex2;
			s1     = gsl_multimin_fminimizer_alloc(R,parameterSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
		case NELDERMEAD2RAND:
			R      = gsl_multimin_fminimizer_nmsimplex2rand;
			s1     = gsl_multimin_fminimizer_alloc(R,parameterSize); // initalize minimizer
			retval = gsl_multimin_fminimizer_set(s1, &f, gslData, ss);
			break;
#endif
		default:
			cerr << "CoiledCoilFitter::Minimize() undefined minimization algorithm: "<<minimizeAlgorithm<<endl;
			return false;
	}


	if (s1 == NULL){
		cerr << "CoiledCoilFitter::Minimize() problem no minimizer object set up\n";
		return false;
	}


	//  Store positions before minimization .... ????
	bool converged = false;
	
	//cout << "RUN MINIMIZATION"<<endl;
	do {
		iter++;

		// Perform a minimization step
		status = gsl_multimin_fminimizer_iterate(s1);   

		// Check for errors, handling should be better
		if (status != GSL_SUCCESS && status != GSL_CONTINUE) {

	 		cerr << "Error 13425: " << gsl_strerror(status) << endl;
			break;  
		}

		// Get size info, used to see where we are in the minimization
		size = gsl_multimin_fminimizer_size(s1);  
		status = gsl_multimin_test_size(size, tolerance);

		// Check status to see if we have minimized (converged)
		if (status == GSL_SUCCESS) {  
		  converged = true;
		  break;
		}           

		
	} while (status == GSL_CONTINUE && iter < maxIterations);


	// Print out values..
	minimizedValue  = 0.0;
	deltaValue     = 0.0;
	minimizedValue = gsl_multimin_fminimizer_minimum(s1);
	deltaValue     = (minimizedValue - initialValue);

	if (converged){
	  CoiledCoils cc;
	  AtomPointerVector &ats = cc.getCoiledCoilBundle(
					    gsl_vector_get(s1->x,0), //r0
					    1.51, //risePerRes
					    gsl_vector_get(s1->x,1), //pitch
					    2.26, //r1
					    102.8, //w1
					    gsl_vector_get(s1->x,2), //phi1
					    0.0, //dZ
					    numResidues_,
					    symmetry_,
					    numHelices_);
  
  
	  fprintf(stdout,"%8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f\n",
					    gsl_vector_get(s1->x,0), //r0
					    1.51, //risePerRes
					    gsl_vector_get(s1->x,1), //pitch
					    2.26, //r1
					    102.8, //w1
					    gsl_vector_get(s1->x,2), //phi1
					    0.0 //dZ
		  );



	  Transforms t;
	  t.rmsdAlignment(ats,pAtoms.getAtomPointers());

	  PDBWriter pout;
	  pout.open("/tmp/ccf_min.pdb");
	  pout.write(ats);
	  pout.close();

	  // Compute and return rmsd
	  double rmsd = pAtoms.getAtomPointers().rmsd(ats);
	  fprintf(stdout, "RSMD Final = %8.3f\n",rmsd);

	  // Store these parameters
	  params.push_back(gsl_vector_get(s1->x,0));
	  params.push_back(1.51); //risePerRes
	  params.push_back(gsl_vector_get(s1->x,1)); //pitch
	  params.push_back(2.26); //r1
	  params.push_back(102.8); //w1
	  params.push_back(gsl_vector_get(s1->x,2)); //phi1
	  params.push_back(0.0); //dZ


	} else {
	  cerr << "ERROR CoiledCoilFitter did not converge!\n";
	}



	gsl_vector_free(ss);
	gsl_vector_free(gslData);

	// clean up and free memory
	gsl_multimin_fminimizer_free(s1);      

	return true;
	
}


double CoiledCoilFitter::my_static_f(const gsl_vector *_xvec_ptr, void *_params){
	double temp;
	//((CoiledCoilFitter *)_params)->resetCoordinates(_xvec_ptr);// why?
	temp = ((CoiledCoilFitter *)_params)->my_f(_xvec_ptr,NULL);
	return temp;
}

double CoiledCoilFitter::my_f(const gsl_vector *_xvec_ptr, void *_params){

  // Create coiled coil given gsl_vector parameters and number of helices
  CoiledCoils cc;
  AtomPointerVector &ats = cc.getCoiledCoilBundle(
					    gsl_vector_get(_xvec_ptr,0), //super-helical radius
					    1.51, //risePerRes
					    gsl_vector_get(_xvec_ptr,1), //super-helical pitch
					    2.26, //r1
					    102.8, //w1
					    gsl_vector_get(_xvec_ptr,2), //alpha-helical phase
					    0.0, //dZ
					    numResidues_,
					    symmetry_,
					    numHelices_);

  PDBWriter pout;
  pout.open(MslTools::stringf("/tmp/min_%025d.pdb",minIndex++));
  pout.write(ats);
  pout.close();

  Transforms t;
  t.rmsdAlignment(ats,pAtoms.getAtomPointers());
  

  // Compute and return rmsd
  double rmsd = pAtoms.getAtomPointers().rmsd(ats);
  fprintf(stdout, "RSMD[%20d] = %8.3f\n",minIndex,rmsd);
  return rmsd;

}


