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

#include "CoiledCoils.h"

using namespace MSL;
using namespace std;


/*
  Should go in CoiledCoils.h?

  Conversion from old notation:

  rSuper         = r0             Super Helical Radius
  Pitch          = p0             Super Helical Pitch
  helicalPhase   = ??


  theta          ~ w0             Angular Frequency , radians per residue. (102 = 2*pi / 3.5 )

  risePerResidue = risePerRes                                              (1.5115)
  rAlpha         = _r1            Helical Radius 


  Conversion from Gevorg/Crick to Dan/North notation:
  Gevorg/Crick:
  R0 = super helical radius (Angstroms)
  R1 = alpha helical radius (Angstroms)
  w0 = super helical frequency (radians/residue)
  w1 = alpha helical frequency (radians/residue)
  alpha = helix crossing angle (radians)
  phi1 = alpha helical phase (radians)

  Dan/North:
  r0 = super helical radius (Angstroms)
  r1 = alpha helical radius (Angstroms)
  theta = alpha helical frequency (radians/residue)
  phase = alpha helical phase (radians)
  h = rise/residue in alpha helix (Angstroms)
  pitch = super helical pitch, distance between repeats in superhelix (Angstroms)

  Conversion:
  r0 = R0
  r1 = R1
  theta = w1
  phase = phi1

  h = -(w0 * R0) / sin(alpha)
  pitch = 2 * pi * R0 / tan(alpha)

 */

CoiledCoils::CoiledCoils(){
	sys = new System();
}

CoiledCoils::~CoiledCoils(){
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
	if (sys != NULL){
		delete(sys);
	}
}

/*
  CRICKS PARAMETERIZATION

  Note: Only CA atoms.
  Ref for Cricks parameterization?

 */
void CoiledCoils::cricksCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1){
	atoms.clear();


	double w0 = 2*M_PI*1.479 / _p0;

	for (unsigned int i=1;i<=_nRes; i++) {

		Atom *CA = new Atom();
		CA->setResidueName("ALA");
		CA->setName("CA");
		CA->setElement("C");
		CA->setResidueNumber(i+1);
		CA->setChainId("A");

		double alpha = atan(2*M_PI*_r0 / _p0);
		double x     =   _r0 * cos(w0 * i)  + _r1 * cos(w0*i)  * cos(_w1*i)  + _r1 * cos(alpha)*sin(w0*i)*sin(_w1*i);
		double y     = -(_r0 * sin(w0 * i)  + _r1 * sin(w0*i)  * cos(_w1*i)) + _r1 * cos(alpha)*cos(w0*i)*sin(_w1*i);
		double z     =   _p0*w0*i/(2*M_PI)  + _r1 * sin(alpha) * sin(_w1*i);

		CA->setCoor(x,y,z);

		atoms.push_back(CA);
		CA = NULL;
	}


	

}

/*

  OFFERS PARAMETERIZATION

  Note: Only CA atoms.
  Ref for Offers parameterization?

 */
void CoiledCoils::offersCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1){


	atoms.clear();

	double alpha = atan(2*M_PI*_r0 / _p0);
	double w0    = 2*M_PI / (sqrt(pow(_p0,2)+ pow((2*M_PI*_r0),2)));


	for (unsigned int i=1;i<=_nRes; i++) {

		double t = (i - 1)*_risePerRes;
		double theta = _w1 + 4*M_PI * (i-1) / 7;

		Atom *CA = new Atom();
		CA->setResidueName("ALA");
		CA->setName("CA");
		CA->setElement("C");
		CA->setResidueNumber(i+1);
		CA->setChainId("A");

		double x = _r0   * cos(w0 * t)  + _r1 * cos(w0*t)  * cos(theta)  + _r1 * cos(alpha)*sin(w0*t)*sin(theta);
		double y = -(_r0 * sin(w0 * t ) + _r1 * sin(w0*t)  * cos(theta)) + _r1 * cos(alpha)*cos(w0*t)*sin(theta);
		double z = _p0*w0*t/(2*M_PI)    + _r1 * sin(alpha) * sin(theta);

		CA->setCoor(x,y,z);

		atoms.push_back(CA);
		CA  = NULL;
	}
	
}

/*
  Ben North CoiledCoil Formulation
  Ref.
  NOTE: ALL Atoms, CB built directly from backbone atoms of coil.

*/

void CoiledCoils::northCoiledCoils(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _theta, double _helicalPhase){
	atoms.clear();

	Atom *N, *CA, *CB, *C, *O;
	N = CA = CB = C = O = NULL;

	for (uint i = 0 ; i < _nRes;i++){

		// calc some constants
		double walpha = -1.0*(_theta*2*M_PI/360)*sqrt(pow(_r0,2)+pow((_p0/(2*M_PI)),2))/_risePerRes;
		double t1 = _theta*M_PI*(1.0-double(_nRes))/(360.0*walpha);
		double t2 = (_theta*2.0*M_PI/360.0)/walpha;
		double t = t1 + (double(i))*t2;


		double Nt     = t + (-1.015 / walpha);
		double Nphase = _helicalPhase + 0.555*180/M_PI;
		double Nalpha = _r1 + -0.75;

		N = new Atom();
		N->setResidueName("ALA");
		N->setName("N");
		N->setElement("N");
		N->setResidueNumber(i+1);
		N->setChainId("A");

		N->setCoor(generateNorthCoor(Nt, Nalpha, walpha, Nphase,_r0,_p0));


		CA = new Atom();
		CA->setResidueName("ALA");
		CA->setName("CA");
		CA->setElement("C");
		CA->setResidueNumber(i+1);
		CA->setChainId("A");
		CA->setCoor(generateNorthCoor(t, _r1, walpha, _helicalPhase,_r0,_p0));

		// C
		double Ct     = t + (1.19 / walpha);
		double Cphase = _helicalPhase + -0.71*180/M_PI;
		double Calpha = _r1 + -0.54;

		C = new Atom();
		C->setResidueName("ALA");
		C->setName("C");
		C->setElement("C");
		C->setResidueNumber(i+1);
		C->setChainId("A");

		C->setCoor(generateNorthCoor(Ct, Calpha, walpha, Cphase,_r0,_p0));

		// O
		double Ot = t + (2.512 / walpha);
		double Ophase = _helicalPhase + -2.092*180/M_PI;
		double Oalpha = _r1 + -0.18;

		O = new Atom();
		O->setResidueName("ALA");
		O->setName("O");
		O->setElement("O");
		O->setResidueNumber(i+1);
		O->setChainId("A");
		O->setCoor(generateNorthCoor(Ot, Oalpha, walpha, Ophase,_r0,_p0));


		// CB Build from N,CA,C atoms
		CB = new Atom(*CA);
		CB->setResidueName("ALA");
		CB->setName("CB");
		CB->setElement("C");
		CB->setChainId("A");
		CB->setResidueNumber(i+1);
		CB->setCoor(CartesianGeometry::instance()->build(CA->getCoor(), N->getCoor(), C->getCoor(), 1.521, 110.5, -122.5));			


		atoms.push_back(N);
		atoms.push_back(CA);
		atoms.push_back(CB);
		atoms.push_back(C);
		atoms.push_back(O);

		N = CA = CB = C = O = NULL;
	}
}


CartesianPoint CoiledCoils::generateNorthCoor(double _t, double _r1, double _wAlpha, double _helicalPhase, double _r0, double _p0){

	//double z = _p0*(-_t/(2*M_PI))+ (_r1*2*M_PI*_r0/sqrt(pow((2*M_PI*_r0),2) + pow(_p0,2)))*sin(_wAlpha*_t+_helicalPhase*2*M_PI/360) + zTrans; 
	double z = _p0*(-_t/(2*M_PI))+ (_r1*2*M_PI*_r0/sqrt(pow((2*M_PI*_r0),2) + pow(_p0,2)))*sin(_wAlpha*_t+_helicalPhase*2*M_PI/360);

	double variable_r0 = _r0;

	double x = variable_r0*cos(_t) + _r1*cos(_t)*cos(_wAlpha*_t+_helicalPhase*2*M_PI/360)- _r1*(_p0/sqrt(pow((2*M_PI*variable_r0),2)+pow(_p0,2)))* sin(_t)*sin(_wAlpha*_t+_helicalPhase*2*M_PI/360);
	double y = variable_r0*sin(_t) + _r1*sin(_t)*cos(_wAlpha*_t+_helicalPhase*2*M_PI/360)+ _r1*(_p0/sqrt(pow((2*M_PI*variable_r0),2)+pow(_p0,2)))* cos(_t)*sin(_wAlpha*_t+_helicalPhase*2*M_PI/360);

	CartesianPoint p(x,y,z);
	return p;
	
}


void CoiledCoils::applyCoiledCoil(AtomPointerVector &_av, double _p0 ){

	atoms.clear();
	for (uint i = 0; i < _av.size(); i++){

		double theta = -2 * M_PI * _av(i).getZ() / _p0;

		double x = _av(i).getX() * cos(theta) - _av(i).getY() * sin(theta);
		double y = _av(i).getX() * sin(theta) + _av(i).getY() * cos(theta);
		
		Atom *a = new Atom(_av(i));
		a->setCoor(x,y,_av(i).getZ());
		atoms.push_back(a);
		
	}

}
void CoiledCoils::sotoCoiledCoils(double _r0, double _risePerRes, double _r1, int _nRes, double _resPerTurn, double _alpha, double _helicalPhase){
	double OmegaSuper;
	double OmegaAlpha;
	double Temp[3];
	double x,y,z;
	double domega,dr,dz;
	double dt,Pitch=0.00;
	int t,i;

	bool POSITIVE=false;
	/*************************************************************
	 * Be certain to decide whether you are going 
	 * to be generating a left handed (negative) or right handed coil    
	 *************************************************************/  
	if (POSITIVE){
		_alpha = _alpha*M_PI/180; Pitch = (2*M_PI*_r0)/(tan(_alpha));
		OmegaSuper = (2*M_PI)/(Pitch/_risePerRes);
		OmegaAlpha = (2*M_PI)/(_resPerTurn);
		_helicalPhase *=M_PI/180;
	}else{
		_alpha = -1.0*_alpha*M_PI/180;
		Pitch = (2*M_PI*_r0)/(tan(_alpha));
		OmegaSuper = (2*M_PI)/(((2*M_PI*_r0)/(tan(_alpha)))/_risePerRes);
		OmegaAlpha = (2*M_PI)/(_resPerTurn);
		_helicalPhase *=M_PI/180;
	}


	atoms.clear();
	Atom *N, *CA, *CB, *C, *O;
	N = CA = CB = C = O = NULL;
	for(t=-_nRes/2,i=0;t<_nRes/2+(_nRes%2);t++,i=i+12){

		/***************************
		 * For the nitrogen atom   *
		 * domega = -26.89*PI/180  *
		 * dr = -0.77              *
		 * dz = -0.90              *
		 ***************************/
		domega = -26.89*M_PI/180;
		dr     =  0.77;
		dz     = -0.90;

		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha-domega));
		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha-domega)-_helicalPhase) + (dz)*sin(_alpha);
		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha-domega)-_helicalPhase) - (dz)*cos(_alpha);

		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);

		N = new Atom();
		N->setResidueName("ALA");
		N->setName("N");
		N->setElement("N");
		N->setResidueNumber(i+1);
		N->setChainId("A");
		N->setCoor(x,y,z);

		/**********
		 *   CA   *
		 **********/
		x =         _r0*cos(t*OmegaSuper) + _r1*cos(t*OmegaSuper)*cos(t*OmegaAlpha-_helicalPhase) - _r1*sin(t*OmegaSuper)*cos(_alpha)*sin(t*OmegaAlpha-_helicalPhase);
		y =         _r0*sin(OmegaSuper*t) + _r1*cos(t*OmegaAlpha-_helicalPhase)*sin(OmegaSuper*t) +  _r1*sin(t*OmegaAlpha-_helicalPhase)*cos(t*OmegaSuper)*cos(_alpha);
		z =    -1.0*_r1*sin(_alpha)*sin(t*OmegaAlpha-_helicalPhase) +  (_r0*OmegaSuper*t)/tan(_alpha); 

		CA = new Atom();
		CA->setResidueName("ALA");
		CA->setName("CA");
		CA->setElement("C");
		CA->setResidueNumber(i+1);
		CA->setChainId("A");
		CA->setCoor(x,y,z);


		/*********************************
		 * For the carboxyl carbon       *
		 * delta omega = 27.42*PI/180    *
		 * delta r = -0.55               *
		 * delta z = 1.05                *
		 *********************************/
		domega = 27.42*M_PI/180;
		dr     = 0.55;
		dz     = 1.05;

		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha+domega));
		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) + (dz)*sin(_alpha);
		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) - (dz)*cos(_alpha);
	  
		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);

		C = new Atom();
		C->setResidueName("ALA");
		C->setName("C");
		C->setElement("C");
		C->setResidueNumber(i+1);
		C->setChainId("A");

		C->setCoor(x,y,z);

		/**********************************
		 * For the oxygen atom            *
		 * delta omega = 23.85*PI/180     *
		 * delta r = -0.19                *
		 * delta z = 2.23                 *
		 **********************************/
		domega = 23.85*M_PI/180;
		dr     = 0.19;
		dz     = 2.23;

		     
		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha+domega));
		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) + (dz)*sin(_alpha);
		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) - (dz)*cos(_alpha);

		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);


		O = new Atom();
		O->setResidueName("ALA");
		O->setName("O");
		O->setElement("O");
		O->setResidueNumber(i+1);
		O->setChainId("A");
		O->setCoor(x,y,z);


		// CB Build from N,CA,C atoms
		CB = new Atom(*CA);
		CB->setResidueName("ALA");
		CB->setName("CB");
		CB->setElement("C");
		CB->setChainId("A");
		CB->setResidueNumber(i+1);
		CB->setCoor(CartesianGeometry::instance()->build(CA->getCoor(), N->getCoor(), C->getCoor(), 1.521, 110.5, -122.5));			


		atoms.push_back(N);
		atoms.push_back(CA);
		atoms.push_back(CB);
		atoms.push_back(C);
		atoms.push_back(O);

		N = CA = CB = C = O = NULL;
  
	}

}
void CoiledCoils::useBBQTable(string _bbqTable){
	bbqTable.openReader(_bbqTable);
}
Chain & CoiledCoils::getChain(){


	if (sys != NULL){
		delete(sys);
		sys = new System();
	}
	sys->addAtoms(atoms);
	bbqTable.fillInMissingBBAtoms((*sys)("A"));
	return (*sys)("A");
}

System & CoiledCoils::getSystem(){

	if (sys != NULL){
		delete(sys);
		sys = new System();
	}
	sys->addAtoms(atoms);
	bbqTable.fillInMissingBBAtoms((*sys)("A"));
	return (*sys);
}
