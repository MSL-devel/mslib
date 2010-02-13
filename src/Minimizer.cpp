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

#include "Minimizer.h"

template <class T> Minimizer<T>::Minimizer(){
	ptTObj = NULL;
}

template <class T> Minimizer<T>::Minimizer(T* _ptTObj){
	setObjectInstance(_ptTObj);
}

template <class T> Minimizer<T>::~Minimizer(){
	for (uint i = 0; i < springControlledAtoms.size();i++){
		delete(springControlledAtoms[i]);
	}

	springControlledAtoms.clear();
}


template <class T> void Minimizer<T>::Minimize(){

}


template <class T> void Minimizer<T>::addData(AtomPointerVector  &_av){


	for (AtomPointerVector::iterator it = _av.begin(); it != _av.end();it++){
		data.push_back((**it).getCoor().getXptr());
		data.push_back((**it).getCoor().getYptr());
		data.push_back((**it).getCoor().getZptr());
	}

}


template <class T> void Minimizer<T>::printData(){
	for (uint i = 0; i < data.size();i++){
		fprintf(stdout, "%8.3f, ",*data[i]);
	}
	fprintf(stdout, "\n");
}


template <class T> void Minimizer<T>::freezeAtoms(AtomPointerVector &_av, double _springConstant){


	for (AtomPointerVector::iterator it = _av.begin(); it != _av.end();it++){	

		Atom *a = new Atom(**it);
		springControlledAtoms.push_back(a);
		//springEnergy.addCharmmBondTerm(*a,**it,0.0, _springConstant);
		a = NULL;
	}

	
}

//  INLINES
template<class T> void Minimizer<T>::addData(double &_data) { data.push_back(&_data); }
template<class T> void Minimizer<T>::setStepSize(double _stepsize) { stepsize = _stepsize; }
template<class T> double Minimizer<T>::getStepSize()               { return stepsize; }
	       
template<class T> double Minimizer<T>::getTolerance()          { return tolerance; }
template<class T> void Minimizer<T>::setTolerance(double _tol)     { tolerance = _tol; }
	       
template<class T> int Minimizer<T>::getMaxIterations()              { return maxIterations; }
template<class T> void Minimizer<T>::setMaxIterations(int _maxIter)  { maxIterations = _maxIter; }
	       
template<class T> int Minimizer<T>::getMinimizeAlgorithm()              { return minimizeAlgorithm; }
template<class T> void Minimizer<T>::setMinimizeAlgorithm(int _algo)     { minimizeAlgorithm = _algo; }
	       
template<class T> vector<double *>& Minimizer<T>::getData()           { return data;  }
	       
template<class T> void Minimizer<T>::setFunction(double (T::*_func)())   { func = _func; }
//templa<class T> double (*Minimizer<T>::getFunction())()             { return func;}
	       
template<class T> void Minimizer<T>::setDerivative(double (T::*_df)())  { df = _df; }
//template<class T> double (*Minimizer<T>::getDerivative())()          { return df;}

template<class T> double Minimizer<T>::CallFunc() { return (*ptTObj.*func)(); }

// Template Stuff
class logTest;
template class Minimizer<logTest>;

class EnergySet;
template class Minimizer<EnergySet>;
