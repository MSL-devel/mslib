/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#ifndef CHARMMENERGY_H
#define CHARMMENERGY_H

#include <iostream>
#include <math.h>


namespace MSL { 
class CharmmEnergy {
	public:
		static CharmmEnergy * instance();

		double LJ(double _d, double _rmin, double _Emin) const;
		double switchingFunction(double _d, double _rOn, double _rOff) const;
		double spring(double _d, double _Kd, double _d0) const;
		double coulombEner(double _d, double _q1, double _q2, double _diel, double _rescalingFactor) const; 

		// merely a division function, for precomputed [rescal*kq*q1*q2/diel] / distance or distance^2 in case of
		// r-dielectric
		double coulombEnerPrecomputed(double _d, double _q1_q2_kq_diel_rescal) const; 
//		double coulombEnerRDielPrecomputed(double _d, double _q1_q2_kq_diel_rescal) const; 
		double dihedralEner(double _chiRadians, double _Kchi, double _n, double _deltaRadians) const; 	

		double EEF1Ener(double _d, double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j) const;
		// For ureyBradley and Angle -- Call spring with angle in Radians and appropriate prameters
		
		// parameter settings
		
		/*
		// rescaling of the 1-4 electostatic interactions.  should be 1 for charmm 22
		// and 0.6 for charmm 19
		void setElec14factor(double _e14);
		double getElec14factor() const;

		// the dielectric constant
		void setDielectricConstant(double _diel);
		double getDielectricConstant() const;

		// use a distance dependent dielectric
		void setUseRdielectric(bool _flag);
		bool getUseRdielectric() const;
		*/

                static const double Kq;
		static const double eef1_constant;
	protected:
		// disallow instantiation
		CharmmEnergy();
		// disallow copy
		CharmmEnergy(const CharmmEnergy & _instance);
		void operator =(const CharmmEnergy & _instance);

//		double elec14factor;
//		double dielectricConstant;
//		bool useRdielectric;



};
inline double CharmmEnergy::spring(double _d, double _Kd, double _d0) const { return _Kd * pow((_d - _d0), 2.0);}
//inline void CharmmEnergy::setElec14factor(double _e14) {elec14factor = _e14;}
//inline double CharmmEnergy::getElec14factor() const {return elec14factor;}
//inline void CharmmEnergy::setDielectricConstant(double _diel) {dielectricConstant = _diel;}
//inline double CharmmEnergy::getDielectricConstant() const {return dielectricConstant;}
//inline void CharmmEnergy::setUseRdielectric(bool _flag) {useRdielectric = _flag;}
//inline bool CharmmEnergy::getUseRdielectric() const {return useRdielectric;}
inline double CharmmEnergy::coulombEnerPrecomputed(double _d, double _q1_q2_kq_diel_rescal) const {
	//std::cout << "UUU " << _q1_q2_kq_diel_rescal << " / " <<  _d << std::endl;
	return _q1_q2_kq_diel_rescal / _d;
}
/*
inline double CharmmEnergy::coulombEnerRDielPrecomputed(double _d, double _q1_q2_kq_diel_rescal, double _dummy, double _dummy2, double _dummy3) const {
	return _q1_q2_kq_diel_rescal / (_d * _d);
}
*/




}

#endif
