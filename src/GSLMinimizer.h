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


#include "Minimizer.h"

#ifndef __GSL__
#error message("GSLMinimizer can't compile unless GSL libraries are present")
#endif

#include <gsl/gsl_multimin.h>


template <class T> class GSLMinimizer : public Minimizer<T> {

	public:

		GSLMinimizer();  
		~GSLMinimizer();

		// Minimize
		void Minimize();	

	private:



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

		// GSL Constants for Minimization Algorithms
		const gsl_multimin_fminimizer_type    *R;
		const gsl_multimin_fdfminimizer_type  *F;

		// Double data in GSL
		gsl_vector *gslData;

		// Step size vector in GSL
		gsl_vector *ss;


	
};
#endif
