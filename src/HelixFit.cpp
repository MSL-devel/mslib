#include "HelixFit.h"


HelixFit::HelixFit(){
	moveSet = ALLFIT;
	setup();
}
HelixFit::HelixFit(Chain &_ch){
	ch = &_ch;
	moveSet = ALLFIT;
	setup();
}

HelixFit::HelixFit(AtomVector &_caOnly){
	caOnly = &_caOnly;
	moveSet = ALLFIT;
	setup();
}
HelixFit::~HelixFit(){
	lastFitHelix.deletePointers();

}
void HelixFit::setStepSize(double _stepSize) { 
	stepSize = _stepSize;
}
void HelixFit::setTolerance(double _tolerance){
	tolerance = _tolerance;
}
void HelixFit::setNumberSteps(int _numSteps){
	numSteps = _numSteps;
}

void HelixFit::setVerbose(bool _verbose){
	verbose = _verbose;
}
void HelixFit::setup(){
	minimizedValue = MslTools::doubleMax;
	tolerance = 0.0001;
	stepSize  = 1.0;
	numSteps  = 1000;
	verbose   = false;
	lastFitRmsd = 0.0;
	lastFitSSQ = 0.0;
	lastFitConstraints = 0.0;
	
	rotationMatrix = Matrix(3,3,0.0);
	symmetryMatrix = Matrix(3,3,0.0);
	translationVector = CartesianPoint(0.0,0.0,0.0);

	// Default symmetry matrix to identity.

	symmetryMatrix[0][0] = 1.0;
	symmetryMatrix[1][1] = 1.0;
	symmetryMatrix[2][2] = 1.0;

	if (moveSet == C2FIT){
		parameters.resize(10);
		parameters[0] = 5.0*M_PI/180;
		parameters[1] = 5.0*M_PI/180;
		parameters[2] = 5.0*M_PI/180;
		parameters[3] = 5.0;
		parameters[4] = 5.0;
		parameters[5] = 5.0;
		parameters[6] = 2.0;
		parameters[7] = 1.0;
		parameters[8] = 176.937*M_PI/180;
		parameters[9] = -47.203*M_PI/180;

		angleParmeters.push_back(0);
		angleParmeters.push_back(1);
		angleParmeters.push_back(2);
		angleParmeters.push_back(8);
		angleParmeters.push_back(9);


	} 

	if (moveSet == ALLFIT){
		parameters.resize(6);
		parameters[0] = 5.0*M_PI/180;
		parameters[1] = 5.0*M_PI/180;
		parameters[2] = 5.0*M_PI/180;
		parameters[3] = 5.0;
		parameters[4] = 5.0;
		parameters[5] = 5.0*M_PI/180;


		//angleParmeters.push_back(0);
		//angleParmeters.push_back(1);
		//angleParmeters.push_back(2);
		//angleParmeters.push_back(5);

		constrainParameters.push_back(pair<int,double>(3,10));
		constrainParameters.push_back(pair<int,double>(0,0.75));


	}


}
void HelixFit::setChain(Chain &_ch){
	ch = &_ch;
}

void HelixFit::setAtoms(AtomVector &_caOnly){
	caOnly = &_caOnly;
}

void HelixFit::setRigidMovementSet(int _moveSet){
	moveSet = _moveSet;
}
double HelixFit::my_static_f(const gsl_vector *xvec_ptr, void *params){
	double temp;
	temp = ((HelixFit *)params)->RigidBodyHelixFit(xvec_ptr);
	return temp;
}
void HelixFit::my_static_df(const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr){

	((HelixFit *)params)->RigidBodyHelixFitDerivatives(xvec_ptr,df_ptr);

}

void HelixFit::my_static_fdf(const gsl_vector *xvec_ptr, void *params,double *f_ptr, gsl_vector *df_ptr){

	*f_ptr = ((HelixFit *)params)->RigidBodyHelixFit(xvec_ptr);
	((HelixFit *)params)->RigidBodyHelixFitDerivatives(xvec_ptr,df_ptr);
}


void HelixFit::setParameters(vector<double> _params){
	parameters = _params;
}
bool HelixFit::fit(int minimizeAlgorithm){


	startingParameters = parameters;

	// Multidimensional Minimizable Functions in GSL
	gsl_multimin_function f;
	gsl_multimin_function_fdf fdf;

	// Mulit-dimensional minmizer objects
	gsl_multimin_fminimizer     *s1;
	gsl_multimin_fdfminimizer   *s2;

	// GSL Constants for Minimization Algorithm Types
	const gsl_multimin_fminimizer_type    *R;
	const gsl_multimin_fdfminimizer_type  *T;

	gsl_vector *ss;
	gsl_vector *x;
	int i; 
	size_t iter = 0;
	int status;
	double size;
	int retval;
	retval = 0;

     
	/**********************************
	 * Initialize the starting vector. You should note that you will need
	 * to set these values each time you come into this function.* 
	 **********************************/
	x = gsl_vector_alloc(parameters.size());
	for(i=0;i<parameters.size();i++){
		gsl_vector_set (x, i, parameters[i]); 
	}

	/*********************************
	 * Initialize method and iterate *
	 *********************************/
	f.f      = &HelixFit::my_static_f;
	f.n      = parameters.size();
	f.params = (this);


	fdf.f      = &HelixFit::my_static_f;
	fdf.df     = &HelixFit::my_static_df;
	fdf.fdf    = &HelixFit::my_static_fdf;
	fdf.n      = parameters.size();
	fdf.params = (this);

	//gsl_multimin_fdfminimizer_set(s,&function,x, 0.01,1e-3);


	R = NULL;
	T = NULL;
	s1 = NULL;
	s2 = NULL;
	ss = NULL;

	// set up minimizer
	switch (minimizeAlgorithm) {
	case NELDERMEAD1:
		ss = gsl_vector_alloc (parameters.size());
		gsl_vector_set_all (ss, stepSize);  //step size
		R      = gsl_multimin_fminimizer_nmsimplex;
		s1     = gsl_multimin_fminimizer_alloc(R,parameters.size()); // initalize minimizer
		retval = gsl_multimin_fminimizer_set(s1, &f, x, ss);
		break;
		//	      case NELDERMEAD2:
		//		      R      = gsl_multimin_fminimizer_nmsimplex2;
		//		      s1     = gsl_multimin_fminimizer_alloc(R,parameters.size()); // initalize minimizer
		//		      retval = gsl_multimin_fminimizer_set(s1, &f, x, ss);
		//		      break;
		//	      case NELDERMEAD2RAND:
		//		      R      = gsl_multimin_fminimizer_nmsimplex2rand;
		//		      s1     = gsl_multimin_fminimizer_alloc(R,parameters.size()); // initalize minimizer
		//		      retval = gsl_multimin_fminimizer_set(s1, &f, x, ss);
		//		      break;
	case CG_FLETCHER_REEVES:
		T      = gsl_multimin_fdfminimizer_conjugate_fr;
		s2     = gsl_multimin_fdfminimizer_alloc(T,parameters.size()); // initalize minimizer
		retval = gsl_multimin_fdfminimizer_set(s2, &fdf, x, stepSize, tolerance);
		break;
	case CG_POLAK_RIBIERE:
		T      = gsl_multimin_fdfminimizer_conjugate_pr;
		s2     = gsl_multimin_fdfminimizer_alloc(T,parameters.size()); // initalize minimizer
		retval = gsl_multimin_fdfminimizer_set(s2, &fdf, x, stepSize, tolerance);
		break;
	case BFGS:
		T      = gsl_multimin_fdfminimizer_vector_bfgs2;
		s2     = gsl_multimin_fdfminimizer_alloc(T,parameters.size()); // initalize minimizer
		retval = gsl_multimin_fdfminimizer_set(s2, &fdf, x, stepSize, tolerance);
		break;
	case STEEPEST_DESCENT:
		T      = gsl_multimin_fdfminimizer_steepest_descent;
		s2     = gsl_multimin_fdfminimizer_alloc(T,parameters.size()); // initalize minimizer
		retval = gsl_multimin_fdfminimizer_set(s2, &fdf, x, stepSize, tolerance);
		break;

	default:
		cerr << "GSLMinimizer::Minimize() undefined minimization algorithm: "<<minimizeAlgorithm<<endl;
		return false;
	}


	if (s1 == NULL && s2 == NULL){
		cerr << "GSLMinimizer::Minimize() problem no minimizer object set up\n";
		return false;
	}

	bool derivateMinimization = false;
	if (s2 != NULL){
		derivateMinimization = true;
	}



	bool converged = false;
	do{
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
			if (verbose){
				printf ("converged to minimum\n"); 
			}
			converged = true;
		}           

		
		if (verbose){
			fprintf(stderr,"Step: %4d\t",iter);
			if (derivateMinimization) { 

				/*
				  Zrot   = RotTheta
				  Ztrans = Vz
				  Xtrans = Vx
				  Xrot   = RotPhi
				*/

				if (moveSet == C2FIT){
					fprintf (stderr, "RotZ=%6.3f RotY=%6.3f RotX=%6.3f Tx=%6.3f Ty=%6.3f Tz=%6.3f Vx=%6.3f Vz=%6.3f RotTheta=%6.3f RotPhi=%6.3f SSQ= %.3f\n",
					 gsl_vector_get (s2->x, 0)*180/M_PI,
					 gsl_vector_get (s2->x, 1)*180/M_PI,
					 gsl_vector_get (s2->x, 2)*180/M_PI,
					 gsl_vector_get (s2->x, 3),
					 gsl_vector_get (s2->x, 4),
					 gsl_vector_get (s2->x, 5),
					 gsl_vector_get (s2->x, 6),
					 gsl_vector_get (s2->x, 7),
					 gsl_vector_get (s2->x, 8)*180/M_PI,
					 gsl_vector_get (s2->x, 9)*180/M_PI,s2->f);
				} else{
					
					fprintf (stderr, "Zrot(Gamma)=%6.3f Xrot(Alpha)=%6.3f  Yrot(Beta)=%6.3f Tx=%6.3f Tz=%6.3f RotTheta=%6.3f SSQ= %.3f\n",
					 gsl_vector_get (s2->x, 2)*180/M_PI,
					 gsl_vector_get (s2->x, 0)*180/M_PI,
					 gsl_vector_get (s2->x, 1)*180/M_PI,
					 gsl_vector_get (s2->x, 3),
					 gsl_vector_get (s2->x, 4),
					 gsl_vector_get (s2->x, 5)*180/M_PI,
					 s2->f);

				}

			} else {

				if (moveSet == C2FIT){
					fprintf (stderr, "RotZ=%6.3f RotY=%6.3f  RotX=%6.3f Tx=%6.3f Ty=%6.3f Tz=%6.3f Vx=%6.3f Vz=%6.3f RotTheta=%6.3f RotPhi=%6.3f VAL= %.3f\n",
					 gsl_vector_get (s1->x, 0)*180/M_PI,
					 gsl_vector_get (s1->x, 1)*180/M_PI,
					 gsl_vector_get (s1->x, 2)*180/M_PI,
					 gsl_vector_get (s1->x, 3),
					 gsl_vector_get (s1->x, 4),
					 gsl_vector_get (s1->x, 5),
					 gsl_vector_get (s1->x, 6),
					 gsl_vector_get (s1->x, 7),
					 gsl_vector_get (s1->x, 8)*180/M_PI,
					 gsl_vector_get (s1->x, 9)*180/M_PI,
					 s1->fval);
				} else {

					fprintf (stderr, "Zrot(Gamma)=%6.3f Xrot(Alpha)=%6.3f  Yrot(Beta)=%6.3f Tx=%6.3f Tz=%6.3f RotTheta=%6.3f VAL= %.3f\n",
					 gsl_vector_get (s1->x, 2)*180/M_PI,
					 gsl_vector_get (s1->x, 0)*180/M_PI,
					 gsl_vector_get (s1->x, 1)*180/M_PI,
					 gsl_vector_get (s1->x, 3),
					 gsl_vector_get (s1->x, 4),
					 gsl_vector_get (s1->x, 5)*180/M_PI,
					 s1->fval);
				}

			}
		}

	}while (status == GSL_CONTINUE && iter < numSteps);
       
	if (derivateMinimization){
			minimizedValue = gsl_multimin_fdfminimizer_minimum(s2);
			gsl_vector *data = gsl_multimin_fdfminimizer_x(s2);
			for (uint i =0;i<parameters.size();i++){
				parameters[i] = gsl_vector_get(data,i);
			}
			data = NULL;
	} else {
			minimizedValue = gsl_multimin_fminimizer_minimum(s1);
			gsl_vector *data = gsl_multimin_fminimizer_x(s1);
			for (uint i =0;i<parameters.size();i++){
				parameters[i] = gsl_vector_get(data,i);
			}
			data = NULL;
	}

	
	// Print out values..
	if (verbose){
		if (derivateMinimization){
			fprintf(stdout,"==> Minimized Value: %8.3f\n",minimizedValue);
		} else {
			fprintf(stdout,"==> Minimized Value: %8.3f.\n",minimizedValue);
		}
	}

	// I may have to retrieve lowest energy coordinates and set the atoms here...

	gsl_vector_free(x);
	if (ss != NULL){
		gsl_vector_free(ss);
	}

	// clean up and free memory
	if (derivateMinimization){
		gsl_multimin_fdfminimizer_free(s2);   
	} else {
		gsl_multimin_fminimizer_free(s1);      
	}

	R = NULL;
	T = NULL;

	return converged;

}


double HelixFit::RigidBodyHelixFit(const gsl_vector *_variablesToOptimize){
	
	for (uint i = 0; i < parameters.size();i++){
		parameters[i] = gsl_vector_get(_variablesToOptimize,i);
	}


	// If defined angle parameters, insure 0 to 360.
	for (uint i = 0; i < angleParmeters.size();i++){
		if (parameters[angleParmeters[i]] < 0){
			parameters[angleParmeters[i]] += 2*M_PI;
			gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),angleParmeters[i],parameters[angleParmeters[i]]);
		}

		if (parameters[angleParmeters[i]] > 2*M_PI){
			parameters[angleParmeters[i]] -= 2*M_PI;
			gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),angleParmeters[i],parameters[angleParmeters[i]]);
		}
	}

	// Create rotatationMatrix, translationVector and symmetryMatrix
	createTransformation(parameters);


	int numResidues = 0;
	if (ch != NULL){
		numResidues = ch->size();
	}

	if (caOnly != NULL){
		numResidues = caOnly->size();
	}


	double SSQ = 0.0;

	// Add in constraints
	for (uint i = 0; i < constrainParameters.size();i++){
		int index = constrainParameters[i].first;
		double force = constrainParameters[i].second;
		double con = force*(parameters[index]-startingParameters[index])*(parameters[index]-startingParameters[index]);
		SSQ += con;
	}


	// Add in distance-squares
	CartesianPoint tmp(0.0,0.0,0.0);
	for(uint i=0; i<numResidues;i++){

		CartesianPoint ca;
		if (ch != NULL){
			ca = ch->getResidueByIndex(i).getAtom("CA").getCoor();
		}

		if (caOnly != NULL){
			ca = (*caOnly)(i).getCoor();
		}

		/***********************
		 * Equations for helix *
		 **********************/
		CartesianPoint ideal = createIdealHelixResidue(parameters,(double)i);
		
		
		tmp = ideal*rotationMatrix + translationVector;
		tmp *= symmetryMatrix;


		double diffSq = ((tmp.getX()-ca.getX())*(tmp.getX()-ca.getX())) +
			        ((tmp.getY()-ca.getY())*(tmp.getY()-ca.getY())) +
			        ((tmp.getZ()-ca.getZ())*(tmp.getZ()-ca.getZ()));

		SSQ += diffSq ;

	}


	return(SSQ);
	

	

}


void HelixFit::RigidBodyHelixFitDerivatives(const gsl_vector *_variablesToOptimize, gsl_vector *df){
	if (moveSet != C2FIT){
		cerr << "HelixFit::RigidBodyHelixFitDerivatives() ONLY C2FIT DERIVATIVES ARE SETUP!!!\n";
		exit(0);
	}

	double A, B, C;
	double X, Y, Z;
	double Alpha, Beta, Gamma;
	double Theta, Phi, Vx, Vz;
	double Tx, Ty, Tz;
	double HX, HY, HZ;
	double dPhi=0.0, dTheta=0.0, dVx=0.0, dVz=0.0;
	double dAlpha=0.0, dBeta=0.0, dGamma=0.0;
	double dTx=0.0, dTy=0.0, dTz=0.0;
	int i,N;


	parameters[0] = Alpha = gsl_vector_get(_variablesToOptimize, 0);
	parameters[1] = Beta  = gsl_vector_get(_variablesToOptimize, 1);
	parameters[2] = Gamma = gsl_vector_get(_variablesToOptimize, 2);
	parameters[3] = Tx = gsl_vector_get(_variablesToOptimize, 3);
	parameters[4] = Ty = gsl_vector_get(_variablesToOptimize, 4);
	parameters[5] = Tz = gsl_vector_get(_variablesToOptimize, 5);
	parameters[6] = Vx = gsl_vector_get(_variablesToOptimize, 6);
	parameters[7] = Vz = gsl_vector_get(_variablesToOptimize, 7);
	parameters[8] = Theta = gsl_vector_get(_variablesToOptimize, 8);
	parameters[9] = Phi = gsl_vector_get(_variablesToOptimize, 9);
	N = parameters.size();

	//How many atoms do you have in your array?

	int numResidues = 0;
	if (ch != NULL){
		numResidues = ch->size();
	}

	if (caOnly != NULL){
		numResidues = caOnly->size();
	}
	for(uint i=0; i<numResidues;i++){

		CartesianPoint ca;
		if (ch != NULL){
			ca = ch->getResidueByIndex(i).getAtom("CA").getCoor();
		}

		if (caOnly != NULL){
			ca = (*caOnly)(i).getCoor();
		}

		A = ca.getX();
		B = ca.getY();
		C = ca.getZ();

		/***********************
		 * Equations for helix *
		 **********************/
		HX= 2.203*cos(98.75*M_PI/180*(double)i);
		HY= 2.203*sin(98.75*M_PI/180*(double)i);
		HZ= 1.53*(double)i;
		/******************************************
		 * Equations for rotation of helix using  *
		 * Bill's sampling method.                *
		 ******************************************/
		X = Vx + HX*cos(Phi)-HY*sin(Phi);
		Y = HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi);
		Z = HZ*cos(Theta) + Vz*cos(Theta)+HY*sin(Theta)*cos(Phi)+HX*sin(Theta)*sin(Phi);
		/************************************
		 * Derivative with respect to Gamma *
		 ************************************/

		dGamma += 2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			       Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			(Z*(cos(Gamma)*sin(Alpha)-sin(Gamma)*cos(Alpha)*sin(Beta))-
			 Y*(cos(Gamma)*cos(Alpha)+sin(Gamma)*sin(Beta)*sin(Alpha))-X*sin(Gamma)*cos(Beta))             +
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			(Z*(sin(Gamma)*sin(Alpha)+cos(Gamma)*cos(Alpha)*sin(Beta)) +
			 Y*(cos(Gamma)*sin(Beta)*sin(Alpha)-sin(Gamma)*cos(Alpha))+X*cos(Gamma)*cos(Beta));


		/***********************************
		 * Derivative with respect to Beta *
		 ***********************************/
 
		dBeta += 2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			      Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			(Y*cos(Gamma)*cos(Beta)*sin(Alpha) - X*cos(Gamma)*sin(Beta) + Z*cos(Gamma)*cos(Beta)*cos(Alpha))  +
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			(Y*sin(Gamma)*cos(Beta)*sin(Alpha)-X*sin(Gamma)*sin(Beta)+Z*sin(Gamma)*cos(Beta)*cos(Alpha))     +
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			(-X*cos(Beta)-Y*sin(Beta)*sin(Alpha)-Z*cos(Alpha)*sin(Beta));

		/************************************
		 * Derivative with respect to Alpha *
		 ************************************/
  
		dAlpha +=   2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
				 Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			(Y*(sin(Gamma)*sin(Alpha)+cos(Gamma)*cos(Alpha)*sin(Beta))+Z*(sin(Gamma)*cos(Alpha)-cos(Gamma)*sin(Beta)*sin(Alpha))) +
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			(-Y*(cos(Gamma)*sin(Alpha)-sin(Gamma)*cos(Alpha)*sin(Beta))-Z*cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha))+
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			(Y*cos(Beta)*cos(Alpha)-Z*cos(Beta)*sin(Alpha));

		/*********************************
		 * Derivative with respect to Tx *
		 *********************************/
		dTx += 2.0*(Tx - A - Y*(sin(Gamma)*cos(Alpha) - cos(Gamma)*sin(Beta)*sin(Alpha)) + 
			    Z*(sin(Gamma)*sin(Alpha)+cos(Gamma)*cos(Alpha)*sin(Beta))+X*(cos(Gamma)*cos(Beta) ));

		/*********************************
		 * Derivative with respect to Ty *
		 *********************************/ 
		dTy += 2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) + 
			    Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)));

		/*********************************
		 * Derivative with respect to Tz *
		 *********************************/
		dTz += 2.0*(Tz-C-X*sin(Beta)+Y*cos(Beta)*sin(Alpha)+Z*cos(Beta)*cos(Alpha));
  
		/**********************************
		 * Derivative with respect to Vx  *
		 **********************************/
		dVx +=   2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			      Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			(1.0*cos(Gamma)*cos(Beta)) +
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			(1.0*sin(Gamma)*cos(Beta))+
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			(-1.0*sin(Beta));

		/**********************************
		 * Derivative with respect to Vz  *
		 **********************************/
		dVz +=   2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			      Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			(sin(Gamma)*sin(Alpha+Theta) + cos(Gamma)*sin(Beta)*cos(Alpha+Theta))
			+
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			(sin(Gamma)*sin(Beta)*cos(Alpha+Theta) - cos(Gamma)*sin(Alpha+Theta))
			+
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			( cos(Beta)*cos(Alpha+Theta));

		/*************************************
		 * Derivative with respect to Theta  *
		 *************************************/
		dTheta +=   2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
				 Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			((-HZ*cos(Theta)-Vz*cos(Theta)-HY*cos(Phi)*sin(Theta)-HX*sin(Theta)*sin(Phi))*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			 (HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi))*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))
			+
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			((-HZ*cos(Theta)-Vz*cos(Theta)-HY*cos(Phi)*sin(Theta)-HX*sin(Theta)*sin(Phi))*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			 (HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi))*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))
			+
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			( (-HZ*cos(Theta)-Vz*cos(Theta)-HY*cos(Phi)*sin(Theta)-HX*sin(Theta)*sin(Phi))*(cos(Beta)*sin(Alpha)) +
			  (HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi))*(cos(Beta)*cos(Alpha)));
		/***********************************
		 * Derivative with respect to Phi  *
		 ***********************************/
		dPhi +=   2.0*(Tx - A + X*cos(Gamma)*cos(Beta) + Y*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			       Z*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))*
			((-HY*cos(Phi)-HX*sin(Phi))*(cos(Gamma)*cos(Beta)) +
			 (HX*cos(Theta)*cos(Phi)-HY*cos(Theta)*sin(Phi))*(cos(Gamma)*sin(Beta)*sin(Alpha) - sin(Gamma)*cos(Alpha)) +
			 (HX*cos(Phi)*sin(Theta)-HY*sin(Theta)*sin(Phi))*(sin(Gamma)*sin(Alpha) + cos(Gamma)*cos(Alpha)*sin(Beta)))

			+
			2.0*(Ty - B + X*sin(Gamma)*cos(Beta) + Y*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			     Z*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))*
			((-HY*cos(Phi)-HX*sin(Phi))*(sin(Gamma)*cos(Beta)) +
			 (HX*cos(Theta)*cos(Phi)-HY*cos(Theta)*sin(Phi))*(cos(Gamma)*cos(Alpha) + sin(Gamma)*sin(Beta)*sin(Alpha)) +
			 (HX*cos(Phi)*sin(Theta)-HY*sin(Theta)*sin(Phi))*(sin(Gamma)*cos(Alpha)*sin(Beta) - cos(Gamma)*sin(Alpha)))
			+
			2.0*(Tz - C - X*sin(Beta) + Y*cos(Beta)*sin(Alpha) + Z*cos(Beta)*cos(Alpha))*
			((-HY*cos(Phi)-HX*sin(Phi))*(sin(Beta)) +
			 (HX*cos(Theta)*cos(Phi)-HY*cos(Theta)*sin(Phi))*(cos(Beta)*sin(Alpha)) +
			 (HX*cos(Phi)*sin(Theta)-HY*sin(Theta)*sin(Phi))*(cos(Beta)*cos(Alpha)));
		/*fprintf(stderr,"j=%d\n",j);
		  fprintf(stderr,"A=%8.3lf\tB=%8.3lf\tC=%8.3lf\n",A, B, C);
		  fprintf(stderr,"Tx=%8.3lf\tTy=%8.3lf\tTz=%8.3lf\n",Tx, Ty, Tz);
		  fprintf(stderr,"Alpha=%8.3lf\tBeta=%8.3lf\tGamma=%8.3lf\n",Alpha*180/M_PI, Beta*180/M_PI, Gamma*180/M_PI);
		  fprintf(stderr,"Vx=%8.3lf\tVz=%8.3lf\tTheta=%8.3lf\tPhi=%8.3lf\n",Vx,Vz,Theta*180/M_PI,Phi*180/M_PI);
		  fprintf(stderr,"dAlpha=%8.3lf\tdBeta=%8.3lf\tdGamma=%8.3lf\n",dAlpha,dBeta,dGamma);
		  fprintf(stderr,"dTx=%8.3lf\tdTy=%8.3lf\tdTz=%8.3lf\n",dTx,dTy,dTz);
		  fprintf(stderr,"dVx=%8.3lf\tdVz=%8.3lf\tdTheta=%8.3lf\tdPhi=%8.3lf\n",dVx,dVz,dTheta,dPhi);
		  exit(1);*/
	}
	gsl_vector_set(df, 0, dAlpha);
	gsl_vector_set(df, 1, dBeta);
	gsl_vector_set(df, 2, dGamma);
	gsl_vector_set(df, 3, dTx);
	gsl_vector_set(df, 4, dTy);
	gsl_vector_set(df, 5, dTz);
	gsl_vector_set(df, 6, dVx);
	gsl_vector_set(df, 7, dVz);
	gsl_vector_set(df, 8, dTheta);
	gsl_vector_set(df, 9, dPhi);
}


AtomVector& HelixFit::fittedHelix(vector<double> &_parameters){
   

	// Create rotatationMatrix, translationVector and symmetryMatrix
	createTransformation(parameters);


	// Delete storage of previous pointers
	lastFitHelix.deletePointers();

	int numResidues = 0;
	if (ch != NULL){
		numResidues = ch->size();
	}

	if (caOnly != NULL){
		numResidues = caOnly->size();
	}

	cout << "Numresidues: "<<numResidues<<endl;

	// Add in constraints
	lastFitConstraints = 0.0;
	if (startingParameters.size() > 0){
		for (uint i = 0; i < constrainParameters.size();i++){
			int index = constrainParameters[i].first;
			double force = constrainParameters[i].second;

			double con = force*(parameters[index]-startingParameters[index])*(parameters[index]-startingParameters[index]);
			lastFitConstraints += con;
		}
	}

	CartesianPoint tmp(0,0,0);

	lastFitSSQ = 0.0;
	for(uint i=0; i<numResidues;i++){

		Atom *ca = NULL;
		if (ch != NULL){
			ca = &ch->getResidueByIndex(i).getAtom("CA");
		}

		if (caOnly != NULL){
			ca = (*caOnly)[i];
		}


		lastFitHelix.push_back(new Atom(*ca));


		/***********************
		 * Equations for helix *
		 **********************/
		CartesianPoint ideal = createIdealHelixResidue(parameters,(double)i);

		//tmp = ideal*rotationMatrix + translationVector;
		//tmp *= symmetryMatrix;
		tmp = ideal;


		(*lastFitHelix.back()).getCoor().setX(tmp.getX());
		(*lastFitHelix.back()).getCoor().setY(tmp.getY());
		(*lastFitHelix.back()).getCoor().setZ(tmp.getZ());


		double diffSq = ((tmp.getX()- ca->getCoor().getX())*(tmp.getX()- ca->getCoor().getX())) +
			        ((tmp.getY()- ca->getCoor().getY())*(tmp.getY()- ca->getCoor().getY())) +
			        ((tmp.getZ()- ca->getCoor().getZ())*(tmp.getZ()- ca->getCoor().getZ()));

		lastFitSSQ += diffSq ;
						  
	}


	double rmsd = 0.0;
	if (ch != NULL){
		lastFitRmsd = ch->getAtoms().rmsd(lastFitHelix);
	}

	if (caOnly != NULL){
		lastFitRmsd = caOnly->rmsd(lastFitHelix);
	}

	

	return lastFitHelix;

}







void HelixFit::createTransformation(vector<double> parameters){



	if (moveSet == C2FIT){

		double Alpha = parameters[0];
		double Beta  = parameters[1];
		double Gamma = parameters[2];
		double Tx    = parameters[3];
		double Ty    = parameters[4];
		double Tz    = parameters[5];
		double Vx    = parameters[6];
		double Vz    = parameters[7];
		double Theta = parameters[8];
		double Phi   = parameters[9];

		rotationMatrix[0][0] = cos(Beta)*cos(Gamma);
		rotationMatrix[1][0] = cos(Beta)*sin(Gamma);
		rotationMatrix[2][0] = -sin(Beta);
		
		rotationMatrix[0][1] = sin(Alpha)*sin(Beta)*cos(Gamma) - cos(Alpha)*sin(Gamma);
		rotationMatrix[1][1] = cos(Alpha)*cos(Gamma) + sin(Alpha)*sin(Beta)*sin(Gamma);
		rotationMatrix[2][1] = sin(Alpha)*cos(Beta);
		
		rotationMatrix[0][2] = sin(Alpha)*sin(Gamma)+cos(Alpha)*sin(Beta)*cos(Gamma);
		rotationMatrix[1][2] = cos(Alpha)*sin(Beta)*sin(Gamma)-sin(Alpha)*cos(Gamma);
		rotationMatrix[2][2] = cos(Alpha)*cos(Beta);


		translationVector[0] = Tx;
		translationVector[1] = Ty;
		translationVector[2] = Tz;

		//symmetryMatrix;
	}

	if (moveSet == ALLFIT){

		double Alpha  = parameters[0];
		double Beta   = parameters[1];
		double Gamma  = parameters[2];
		double Tx     = parameters[3];
		double Tz     = parameters[4];
		double Theta  = parameters[5];

		rotationMatrix[0][0] = cos(Gamma)*cos(Beta)            +  sin(Gamma)*sin(Beta)*sin(Alpha);
		rotationMatrix[1][0] = sin(Gamma)*cos(Alpha);
		rotationMatrix[2][0] = sin(Gamma)*cos(Beta)*sin(Alpha) - cos(Gamma)*sin(Beta);

		rotationMatrix[0][1] = cos(Gamma)*sin(Beta)*sin(Alpha) -  sin(Gamma)*cos(Beta);
		rotationMatrix[1][1] = cos(Gamma)*cos(Alpha);
		rotationMatrix[2][1] = sin(Gamma)*sin(Beta)+cos(Gamma)*cos(Beta)*sin(Alpha);

		rotationMatrix[0][2] = cos(Alpha)*sin(Beta);
		rotationMatrix[1][2] = -sin(Alpha);
		rotationMatrix[2][2] = cos(Beta)*cos(Alpha);


		translationVector[0] = Tx;
		translationVector[1] = 0.0;
		translationVector[2] = Tz;


		symmetryMatrix[0][0] = cos(Theta);
		symmetryMatrix[1][0] = sin(Theta);
		symmetryMatrix[2][0] = 0.0;

		symmetryMatrix[0][1] = -sin(Theta);
		symmetryMatrix[1][1] = cos(Theta);
		symmetryMatrix[2][1] = 0.0;

		symmetryMatrix[0][2] = 0.0;
		symmetryMatrix[1][2] = 0.0;
		symmetryMatrix[2][2] = 1.0;
	}

	
}


CartesianPoint HelixFit::createIdealHelixResidue(vector<double> _parameters, double index){

	CartesianPoint tmp(0.0,0.0,0.0);
	tmp.setX( 2.203*cos(98.75*M_PI/180*index) );
	tmp.setY( 2.203*sin(98.75*M_PI/180*index) );
	tmp.setZ( 1.53*index                      );

	if (moveSet == C2FIT){
		
		CartesianPoint ideal = tmp;

		if (parameters.size() != 10){
			cerr << "HelixFit::createIdeaHelixResidue, not 10 parameters only "<<parameters.size()<<endl;
			exit(0);

		}
		double Vx    = parameters[6];
		double Vz    = parameters[7];
		double Theta = parameters[8];
		double Phi   = parameters[9];


		tmp.setX( (Vx + ideal.getX()*cos(Phi) - ideal.getY()*sin(Phi)));
		tmp.setY( (ideal.getY()*cos(Theta)*cos(Phi)-Vz*sin(Theta)-ideal.getZ()*sin(Theta)+ideal.getX()*cos(Theta)*sin(Phi)));
		tmp.setZ( (ideal.getZ()*cos(Theta)+Vz*cos(Theta)+ideal.getY()*cos(Phi)*sin(Theta)+ideal.getX()*sin(Theta)*sin(Phi)));
		
	}


	return tmp;
}


/*

double HelixFit::RigidBodyHelixFitOld(const gsl_vector *_variablesToOptimize){

	
	Matrix Rotation(3,3,0.0);
	double Alpha, Beta, Gamma;
	double Tx, Ty,Tz;
	double HX, HY, HZ;
	double  X, Y, Z;
	double Vx, Vz, Theta, Phi;
	int i,N;


	parameters[0] = Alpha = gsl_vector_get(_variablesToOptimize, 0);
	parameters[1] = Beta  = gsl_vector_get(_variablesToOptimize, 1);
	parameters[2] = Gamma = gsl_vector_get(_variablesToOptimize, 2);
	parameters[3] = Tx = gsl_vector_get(_variablesToOptimize, 3);
	parameters[4] = Tz = gsl_vector_get(_variablesToOptimize, 4);
	parameters[5] = Theta = gsl_vector_get(_variablesToOptimize, 5);
	N =  parameters.size();


	if (Alpha < 0){
		Alpha += 2*M_PI;
	}

	if (Alpha > 2*M_PI){
		Alpha -= 2*M_PI;
	}

	if (Beta < 0){
		Beta += 2*M_PI;
	}

	if (Beta > 2*M_PI){
		Beta -= 2*M_PI;
	}

	if (Gamma < 0){
		Gamma += 2*M_PI;
	}

	if (Gamma > 2*M_PI){
		Gamma -= 2*M_PI;
	}

	if (Theta < 0){
		Theta += 2*M_PI;
	}

	if (Theta > 2*M_PI){
		Theta -= 2*M_PI;
	}


	gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),0,Alpha);
	gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),1,Beta);
	gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),2,Gamma);
	gsl_vector_set(const_cast<gsl_vector *>(_variablesToOptimize),5,Theta);
	
	parameters[0] = Alpha;
	parameters[1] = Beta;
	parameters[2] = Gamma;
	parameters[5] = Theta;


	Rotation[0][0] = cos(Gamma)*cos(Beta)            +  sin(Gamma)*sin(Beta)*sin(Alpha);
	Rotation[1][0] = sin(Gamma)*cos(Alpha);
	Rotation[2][0] = sin(Gamma)*cos(Beta)*sin(Alpha) - cos(Gamma)*sin(Beta);

	Rotation[0][1] = cos(Gamma)*sin(Beta)*sin(Alpha) -  sin(Gamma)*cos(Beta);
	Rotation[1][1] = cos(Gamma)*cos(Alpha);
	Rotation[2][1] = sin(Gamma)*sin(Beta)+cos(Gamma)*cos(Beta)*sin(Alpha);

	Rotation[0][2] = cos(Alpha)*sin(Beta);
	Rotation[1][2] = -sin(Alpha);
	Rotation[2][2] = cos(Beta)*cos(Alpha);

	//cout << endl<<"Rotation:"<<endl<<Rotation.toString()<<endl;
  
	Matrix sym(3,3,0.0);

	sym[0][0] = cos(Theta);
	sym[1][0] = sin(Theta);

	sym[0][1] = -sin(Theta);
	sym[1][1] = cos(Theta);

	sym[2][2] = 1.0;
	
	CartesianPoint trans(Tx,0,Tz);




	CartesianPoint tmp(0.0,0.0,0.0);

	int numResidues = 0;
	if (ch != NULL){
		numResidues = ch->size();
	}

	if (caOnly != NULL){
		numResidues = caOnly->size();
	}


	double SSQ = 0.0;

	// Add in constraints
	double TxCon = 0.25*(Tx-startingParameters[3])*(Tx-startingParameters[3]);


	SSQ = TxCon;

	// Add in distance-squares
	for(uint i=0; i<numResidues;i++){

		CartesianPoint ca;
		if (ch != NULL){
			ca = ch->getResidueByIndex(i).getAtom("CA").getCoor();
		}

		if (caOnly != NULL){
			ca = (*caOnly)(i).getCoor();
		}

	
		CartesianPoint ideal;
		ideal.setX(2.203*cos(98.75*M_PI/180*(double)i));
		ideal.setY(2.203*sin(98.75*M_PI/180*(double)i));
		ideal.setZ(1.53*(double)i);
		
		// Rotations and translations
		//tmp  = CartesianGeometry::instance()->matrixTimesCartesianPoint(ideal,Rotation) + trans;
		tmp[0]  = Rotation[0][0]*ideal.getX()  + Rotation[0][1]*ideal.getY() + Rotation[0][2]*ideal.getZ() + Tx;
		tmp[1]  = Rotation[1][0]*ideal.getX()  + Rotation[1][1]*ideal.getY() + Rotation[1][2]*ideal.getZ() + 0.0;
		tmp[2]  = Rotation[2][0]*ideal.getX()  + Rotation[2][1]*ideal.getY() + Rotation[2][2]*ideal.getZ() + Tz;

		// Symmtery (only 0,180)
		tmp[0]  = sym[0][0]*tmp.getX()  + sym[0][1]*tmp.getY() + sym[0][2]*tmp.getZ();
		tmp[1]  = sym[1][0]*tmp.getX()  + sym[1][1]*tmp.getY() + sym[1][2]*tmp.getZ();
		tmp[2]  = sym[2][0]*tmp.getX()  + sym[2][1]*tmp.getY() + sym[2][2]*tmp.getZ();



		double diffSq = ((tmp.getX()-ca.getX())*(tmp.getX()-ca.getX())) +
			((tmp.getY()-ca.getY())*(tmp.getY()-ca.getY())) +
			((tmp.getZ()-ca.getZ())*(tmp.getZ()-ca.getZ()));

		//cout << "CA: "<<ca.toString()<<"\t"<<tmp.toString()<<" "<<diffSq<<endl;

		SSQ += diffSq ;

	}


	return(SSQ);

	
}


double HelixFit::RigidBodyHelixFitC2(const gsl_vector *_variablesToOptimize){

	double Rotation[3][3];
	double SSQ=0.00;
	double Alpha, Beta, Gamma;
	double Tx, Ty,Tz;
	double HX, HY, HZ;
	double  X, Y, Z;
	double Vx, Vz, Theta, Phi;
	int i,N;


	parameters[0] = Alpha = gsl_vector_get(_variablesToOptimize, 0);
	parameters[1] = Beta  = gsl_vector_get(_variablesToOptimize, 1);
	parameters[2] = Gamma = gsl_vector_get(_variablesToOptimize, 2);
	parameters[3] = Tx = gsl_vector_get(_variablesToOptimize, 3);
	parameters[4] = Ty = gsl_vector_get(_variablesToOptimize, 4);
	parameters[5] = Tz = gsl_vector_get(_variablesToOptimize, 5);
	parameters[6] = Vx = gsl_vector_get(_variablesToOptimize, 6);
	parameters[7] = Vz = gsl_vector_get(_variablesToOptimize, 7);
	parameters[8] = Theta = gsl_vector_get(_variablesToOptimize, 8);
	parameters[9] = Phi = gsl_vector_get(_variablesToOptimize, 9);
	N =  parameters.size();


	Rotation[0][0] = cos(Beta)*cos(Gamma);
	Rotation[1][0] = cos(Beta)*sin(Gamma);
	Rotation[2][0] = -sin(Beta);

	Rotation[0][1] = sin(Alpha)*sin(Beta)*cos(Gamma) - cos(Alpha)*sin(Gamma);
	Rotation[1][1] = cos(Alpha)*cos(Gamma) + sin(Alpha)*sin(Beta)*sin(Gamma);
	Rotation[2][1] = sin(Alpha)*cos(Beta);

	Rotation[0][2] = sin(Alpha)*sin(Gamma)+cos(Alpha)*sin(Beta)*cos(Gamma);
	Rotation[1][2] = cos(Alpha)*sin(Beta)*sin(Gamma)-sin(Alpha)*cos(Gamma);
	Rotation[2][2] = cos(Alpha)*cos(Beta);
  
	CartesianPoint tmp(0.0,0.0,0.0);

	int numResidues = 0;
	if (ch != NULL){
		numResidues = ch->size();
	}

	if (caOnly != NULL){
		numResidues = caOnly->size();
	}
	for(uint i=0; i<numResidues;i++){

		CartesianPoint ca;
		if (ch != NULL){
			ca = ch->getResidueByIndex(i).getAtom("CA").getCoor();
		}

		if (caOnly != NULL){
			ca = (*caOnly)(i).getCoor();
		}


		HX=2.203*cos(98.75*M_PI/180*(double)i);
		HY=2.203*sin(98.75*M_PI/180*(double)i);
		HZ=1.53*(double)i;


		X = (Vx + HX*cos(Phi) - HY*sin(Phi));
		Y = (HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi));
		Z = (HZ*cos(Theta)+Vz*cos(Theta)+HY*cos(Phi)*sin(Theta)+HX*sin(Theta)*sin(Phi));

		tmp[0]  = Rotation[0][0]*X  + Rotation[0][1]*Y + Rotation[0][2]*Z + Tx;
		tmp[1]  = Rotation[1][0]*X  + Rotation[1][1]*Y + Rotation[1][2]*Z + Ty;
		tmp[2]  = Rotation[2][0]*X  + Rotation[2][1]*Y + Rotation[2][2]*Z + Tz;

		SSQ +=  ((tmp.getX()-ca.getX())*(tmp.getX()-ca.getX())) +
			((tmp.getY()-ca.getY())*(tmp.getY()-ca.getY())) +
			((tmp.getZ()-ca.getZ())*(tmp.getZ()-ca.getZ()));

	}


	return(SSQ);
}


AtomVector& HelixFit::fittedHelix(vector<double> &_parameters){
   
    double Temp[3] = {0.00,0.00,0.00};
    Matrix Rotation(3,3,0.0);
    double SSQ=0.00;
    double Alpha, Beta, Gamma;
	double Tx, Ty,Tz;
	double HX, HY, HZ;
	double Vx, Vz, Theta, Phi;
    int i,j,N;

    Alpha = parameters[0]; 
    Beta =  parameters[1];  
    Gamma = parameters[2]; 
    Tx = parameters[3]; 
    Tz = parameters[4]; 
    Theta = parameters[5]; 


    //You should have a variable to tell you how many atoms you have.
    Rotation[0][0] = cos(Gamma)*cos(Beta)            +  sin(Gamma)*sin(Beta)*sin(Alpha);
    Rotation[1][0] = sin(Gamma)*cos(Alpha);
    Rotation[2][0] = sin(Gamma)*cos(Beta)*sin(Alpha) - cos(Gamma)*sin(Beta);

    Rotation[0][1] = cos(Gamma)*sin(Beta)*sin(Alpha) -  sin(Gamma)*cos(Beta);
    Rotation[1][1] = cos(Gamma)*cos(Alpha);
    Rotation[2][1] = sin(Gamma)*sin(Beta)+cos(Gamma)*cos(Beta)*sin(Alpha);

    Rotation[0][2] = cos(Alpha)*sin(Beta);
    Rotation[1][2] = -sin(Alpha);
    Rotation[2][2] = cos(Beta)*cos(Alpha);

    cout << "RotMatrix:"<<endl<<Rotation.toString()<<endl<<Rotation.getDeterminant()<<endl;
    Matrix sym(3,3,0.0);

    sym[0][0] = cos(Theta);
    sym[1][0] = sin(Theta);

    sym[0][1] = -sin(Theta);
    sym[1][1] = cos(Theta);

    sym[2][2] = 1.0;

    lastFitHelix.deletePointers();

    CartesianPoint trans(Tx,0,Tz);
    CartesianPoint tmp(0,0,0);

    int numResidues = 0;
    if (ch != NULL){
	    numResidues = ch->size();
    }

    if (caOnly != NULL){
	    numResidues = caOnly->size();
    }
    for(uint i=0; i<numResidues;i++){


	    lastFitHelix.push_back(new Atom("CA"));
	    lastFitHelix.back()->setResidueNumber(i);

	    // ***********************
	    // * Equations for helix *
	    // *********************

	    CartesianPoint ideal;
	    ideal.setX(2.203*cos(98.75*M_PI/180*(double)i));
	    ideal.setY(2.203*sin(98.75*M_PI/180*(double)i));
	    ideal.setZ(1.53*(double)i);
	    
	    
	    
	    //tmp = ideal;

	    //tmp[0]  = Rotation[0][0]*ideal.getX()  + Rotation[0][1]*ideal.getY() + Rotation[0][2]*ideal.getZ() + Tx;
	    //tmp[1]  = Rotation[1][0]*ideal.getX()  + Rotation[1][1]*ideal.getY() + Rotation[1][2]*ideal.getZ() + 0.0;
	    //tmp[2]  = Rotation[2][0]*ideal.getX()  + Rotation[2][1]*ideal.getY() + Rotation[2][2]*ideal.getZ() + Tz;


	    //tmp  = CartesianGeometry::instance()->matrixTimesCartesianPoint(ideal,Rotation) + trans;
	    tmp = ideal*Rotation + trans;
	    tmp *= CartesianGeometry::instance()->getZRotationMatrix(Theta*M_PI/180);
	    //tmp[0]  = sym[0][0]*tmp.getX()  + sym[0][1]*tmp.getY() + sym[0][2]*tmp.getZ();
	    //tmp[1]  = sym[1][0]*tmp.getX()  + sym[1][1]*tmp.getY() + sym[1][2]*tmp.getZ();
	    //tmp[2]  = sym[2][0]*tmp.getX()  + sym[2][1]*tmp.getY() + sym[2][2]*tmp.getZ();


	    
	    (*lastFitHelix.back()).getCoor().setX(tmp.getX());
	    (*lastFitHelix.back()).getCoor().setY(tmp.getY());
	    (*lastFitHelix.back()).getCoor().setZ(tmp.getZ());
						  
						  

   }


  return lastFitHelix;

}

// **********************************************************************
// * This function will use the Engel+DeGrado method to generate a helix *
// * and then the helix will be positioned somewhere in space using the  *
// * Euler Angles (Alpha, Beta, Gamma) and a translation.                *
// **********************************************************************
AtomVector& HelixFit::fittedHelixC2(vector<double> &_parameters){
   
    double Temp[3] = {0.00,0.00,0.00};
    double Rotation[3][3];
    double SSQ=0.00;
    double Alpha, Beta, Gamma;
	double Tx, Ty,Tz;
	double HX, HY, HZ;
	double Vx, Vz, Theta, Phi;
    int i,j,N;

    Alpha = parameters[0]; 
    Beta =  parameters[1];  
    Gamma = parameters[2]; 
    Tx = parameters[3]; 
    Ty = parameters[4]; 
    Tz = parameters[5]; 
    Vx = parameters[6]; 
    Vz = parameters[7]; 
    Theta = parameters[8]; 
    Phi = parameters[9]; 
    N = parameters.size();

    //You should have a variable to tell you how many atoms you have.
	
    Rotation[0][0] = cos(Beta)*cos(Gamma);
    Rotation[1][0] = cos(Beta)*sin(Gamma);
    Rotation[2][0] = -sin(Beta);

    Rotation[0][1] = sin(Alpha)*sin(Beta)*cos(Gamma) - cos(Alpha)*sin(Gamma);
    Rotation[1][1] = cos(Alpha)*cos(Gamma) + sin(Alpha)*sin(Beta)*sin(Gamma);
    Rotation[2][1] = sin(Alpha)*cos(Beta);

    Rotation[0][2] = sin(Alpha)*sin(Gamma)+cos(Alpha)*sin(Beta)*cos(Gamma);
    Rotation[1][2] = cos(Alpha)*sin(Beta)*sin(Gamma)-sin(Alpha)*cos(Gamma);
    Rotation[2][2] = cos(Alpha)*cos(Beta);
  
    lastFitHelix.deletePointers();

    for(i=0; i<N;i++){

	    lastFitHelix.push_back(new Atom("CA"));
	    lastFitHelix.back()->setResidueNumber(i);

	    ///***********************
	    // * Equations for helix *
	    // *********************

	    HX= 2.203*cos(98.75*M_PI/180*(double)i);
	    HY= 2.203*sin(98.75*M_PI/180*(double)i);
	    HZ= 1.53*(double)i;

	    (*lastFitHelix.back()).getCoor().setX((Rotation[0][0]*(Vx + HX*cos(Phi) - HY*sin(Phi)) +	
						   Rotation[0][1]*(HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi)) + 
						   Rotation[0][2]*(HZ*cos(Theta)+Vz*cos(Theta)+HY*cos(Phi)*sin(Theta)+HX*sin(Theta)*sin(Phi)) + Tx));

	    (*lastFitHelix.back()).getCoor().setY((Rotation[1][0]*(Vx + HX*cos(Phi) - HY*sin(Phi)) + 
						   Rotation[1][1]*(HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi)) + 
						   Rotation[1][2]*(HZ*cos(Theta)+Vz*cos(Theta)+HY*cos(Phi)*sin(Theta)+HX*sin(Theta)*sin(Phi)) + Ty));

	    (*lastFitHelix.back()).getCoor().setZ((Rotation[2][0]*(Vx + HX*cos(Phi) - HY*sin(Phi)) + 
						   Rotation[2][1]*(HY*cos(Theta)*cos(Phi)-Vz*sin(Theta)-HZ*sin(Theta)+HX*cos(Theta)*sin(Phi)) + 
						   Rotation[2][2]*(HZ*cos(Theta)+Vz*cos(Theta)+HY*cos(Phi)*sin(Theta)+HX*sin(Theta)*sin(Phi)) + Tz));

   }


  return lastFitHelix;

}
*/
