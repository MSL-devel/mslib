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

/********************************************
 *
 *  SYSTEM TO COILED COIL BUNDLE 
 *
 ********************************************/
 //read in parameters which are the atomId of each starting residue of the chain and then the length of each chain
bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions) {
	CartesianPoint cp;
	string atomId;
	int bundleSize = _N;
	if (_symmetry == "D" || _symmetry == "d") {
		bundleSize = 2*_N;
	}
	vector<int> startPositionIndices; 

	//Test to see if the number of parameters matches the bundle size
	//if (_symmetry == "C" || _symmetry == "c") bundleSize = _N; 
	//else 
	if (_startingPositions.size() != bundleSize) {
		cerr << "No. of starting positions does not match no. of bundles" << endl;
		return false;
	}

	//Test to see if all starting positions are contained within the system
	for (int i = 0; i < _startingPositions.size(); i++){
		if (_sys.positionExists(_startingPositions[i])) {
			startPositionIndices.push_back(_sys.getPositionIndex(&_sys.getLastFoundPosition()));
		} else {
			cerr << "Not all starting positions are contained within the bundle" << endl;
			return false;
		}
	}

	//Test to see if starting positions overlap with one another 
	for (int i = 0; i < startPositionIndices.size()-1; i++) {
		if (abs(startPositionIndices[i] - startPositionIndices[i+1]) > _nRes) {
			cerr << "Starting positions overlap with one another" << endl;
			return false;
		}
	}
	if ((startPositionIndices[startPositionIndices.size()-1] + _nRes) > _sys.positionSize()) {
		cerr << "Starting position overlaps with end of the sequence" << endl;
		return false;
	}

	//Create the coiled coil bundle in the atomPointerVector atoms
	getCoiledCoilBundle(_r0, _risePerRes, _pitch, _r1, _w1, _phi1, _dZ, _nRes, _symmetry, _N);

	//Divide out the APV into the chainVector
	vector< vector< map<string, Atom*> > > chainVector = atoms.subdivideByChainAndPosition();

	for (int k = 0; k < bundleSize; k++) {
		int resCount = 0;
		for (int j = startPositionIndices[k]; j < (startPositionIndices[k]+_nRes); j++) {
			for (int i = 0; i < _sys.getPosition(j).atomSize(); i++) {
				string atomName = _sys.getPosition(j).getAtom(i).getName();
				if (atomName == "CA" || atomName == "C" || atomName == "O" || atomName == "N") {
					_sys.getPosition(j).getAtom(i).setCoor(chainVector[k][resCount][atomName]->getCoor());
				}
			}
			resCount++;
		}
	}

	return true;
}

/********************************************
 *
 * GET C2 OR D2 COILED COIL BUNDLE 
 *
 ********************************************/
AtomPointerVector& CoiledCoils::getCoiledCoilBundle(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  string _symmetry, int _N){

	getCoiledCoil(_r0, _risePerRes, _pitch, _r1, _w1, _phi1, _dZ, _nRes);

	Transforms tr;
	Symmetry symm;
	//double angle = 0;

	/***************************
	 * C2 Symmetry
	 **************************/
	if (_symmetry == "c" || _symmetry == "C"){
		symm.applyCN(atoms, _N, CartesianPoint(0.0,0.0,1.0), true);
	}
	/***************************
	 * D2 Symmetry
	 **************************/
	else if (_symmetry == "d" || _symmetry == "D"){
		if (_N%2 == 0) {
			tr.Zrotate(atoms, (45.0 / (_N/2)));
		}
		symm.applyDN(atoms, _N, CartesianPoint(0.0,0.0,1.0), CartesianPoint(0.0,1.0,0.0), true);
	}
	else {
		cerr << "Must specify 'C' or 'D' symmetry" << endl;
	}

	return atoms; 
}

 
AtomPointerVector& CoiledCoils::getCoiledCoilBundleCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes, std::string _symmetry, int _N) {
	double risePerRes = (_w0 * _r0 ) / sin(_a);
	double pitch = 2 * M_PI * _r0 / tan (-_a);

	getCoiledCoilBundle(_r0, risePerRes, pitch, _r1, _w1, _phi1, _dZ, _nRes,  _symmetry, _N);

	return atoms;
}

/********************************************
 *
 *  DEGREE CALL WITH NORTH PARAM
 *
 ********************************************/
AtomPointerVector& CoiledCoils::getCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes){
	double radW1 = _w1 * (M_PI/180);
	double radPhi1 = _phi1 * (M_PI/180);

	radCoiledCoil(_r0, _risePerRes, _pitch, _r1, radW1, radPhi1, _dZ, _nRes);

	return atoms;
}

/********************************************
 *
 *  DEGREE CALL WITH CRICKS PARAM
 *
 ********************************************/
AtomPointerVector& CoiledCoils::getCoiledCoilCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes){
	double radW0 = _w0 * (M_PI/180);
	double radAlpha = _a * (M_PI/180);
	double radW1 = _w1 * (M_PI/180);
	double radPhi1 = _phi1 * (M_PI/180);

	double risePerRes = (radW0 * _r0 ) / sin(radAlpha);
	double pitch = 2 * M_PI * _r0 / tan (-radAlpha);

	radCoiledCoil(_r0, risePerRes, pitch, _r1, radW1, radPhi1, _dZ, _nRes);

	return atoms;
}




/********************************************
 *
 *  Ben North CoiledCoil Formulation
 *
 ********************************************/
void CoiledCoils::radCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes){

	atoms.clear();

	/***************************
	 * Calculated L, 
	 * wa, t 
	 ***************************/
	double L = sqrt(pow(_r0,2)+pow((_pitch /(2*M_PI)),2)); // the length of the spiral
	double wa = -1.0*_w1*L / _risePerRes;
	double t = _risePerRes / (-L);
	double zCenter =  (-1.0)*(_pitch*(-(t*_nRes)/(2*M_PI)) +   _r1*(2*M_PI*_r0/sqrt(pow((2*M_PI*_r0),2) + pow(_pitch,2)))*sin((t*_nRes)*wa + _phi1)) / 2;

	for (uint i = 0; i < _nRes; i++){

		/***************************
		 * Nitrogen placement
		 ***************************/
		Atom *N = new Atom();
		N->setResidueName("ALA");
		N->setName("N");
		N->setElement("N");
		N->setResidueNumber(i+1);
		N->setChainId("A");

		double Nphase = _phi1+ 0.555;
		double nR1 = _r1 + -0.75;
		double corr = -1.015;


		generateNorthCoor(N, i, nR1, wa, Nphase, _r0, _pitch, _risePerRes, L, corr, zCenter);


		/***************************
		 * C-alpha placement
		 ***************************/
		Atom *CA = new Atom();
		CA->setResidueName("ALA");
		CA->setName("CA");
		CA->setElement("C");
		CA->setResidueNumber(i+1);
		CA->setChainId("A");

		corr = 0;


		generateNorthCoor(CA, i, _r1, wa, _phi1, _r0, _pitch, _risePerRes, L, corr, zCenter);

		/***************************
		 * Carbonyl Carbon placement
		 ***************************/
		Atom *C = new Atom();
		C->setResidueName("ALA");
		C->setName("C");
		C->setElement("C");
		C->setResidueNumber(i+1);
		C->setChainId("A");

		double Cphase = _phi1+ -0.71;
		double cR1 = _r1 + -0.54;
		corr = 1.19;


		generateNorthCoor(C, i, cR1, wa, Cphase, _r0, _pitch, _risePerRes, L, corr, zCenter);

		/***************************
		 * Carbonyl Oxygen placement
		 ***************************/
		Atom *O = new Atom();
		O->setResidueName("ALA");
		O->setName("O");
		O->setElement("O");
		O->setResidueNumber(i+1);
		O->setChainId("A");

		double Ophase = _phi1+ -2.092;
		double oR1 = _r1 + -0.18;
		corr = 2.512;


		generateNorthCoor(O, i, oR1, wa, Ophase, _r0, _pitch, _risePerRes, L, corr, zCenter);


		/***************************
		 * Add Atoms to atoms 
		 ***************************/
		atoms.push_back(N);
		atoms.push_back(CA);
		atoms.push_back(C);
		atoms.push_back(O);

		N = CA = C = O = NULL;
	}

	/***************************
	 * Rotate to X-axis 
	 **************************/
	Transforms tr;
	double angle = 0;
	double x = _r0 * cos( ((double)(_nRes)/2.0) *(t));
	double y = _r0 * sin( ((double)(_nRes)/2.0) *(t));
	double quad = 1;

	if (y > 0) {quad *= -1;}

	angle = CartesianGeometry::angle(CartesianPoint(x,y,0.0) , CartesianPoint(0.0,0.0,0.0), CartesianPoint(1.0,0.0,0.0));
	Matrix rotMat = CartesianGeometry::getRotationMatrix(quad*angle, CartesianPoint(0.0,0.0,1.0));

	tr.rotate(atoms, rotMat);
}


void CoiledCoils::generateNorthCoor(Atom * _pAtom, unsigned int i, double _r1, double _wa, double _phi1, double _r0, double _pitch, double _risePerRes, double _L, double _corr, double _zCenter){

	double t = (_risePerRes/(-_L));
	double _c2 = _corr / _wa;

	double x = _r0*cos((t*i+_c2))           +   _r1*cos((t*i+_c2))*cos((t*i+_c2)*_wa + _phi1)      -   _r1*(_pitch/sqrt(pow((2*M_PI*_r0),2)+pow(_pitch,2)))*sin((t*i+_c2))*sin((t*i+_c2)*_wa + _phi1);
	double y = _r0*sin((t*i+_c2))           +   _r1*sin((t*i+_c2))*cos((t*i+_c2)*_wa + _phi1)      +   _r1*(_pitch/sqrt(pow((2*M_PI*_r0),2)+pow(_pitch,2)))*cos((t*i+_c2))*sin((t*i+_c2)*_wa + _phi1);
	double z = _pitch*(-(t*i+_c2)/(2*M_PI)) +   _r1*(2*M_PI*_r0/sqrt(pow((2*M_PI*_r0),2) + pow(_pitch,2)))*sin((t*i+_c2)*_wa + _phi1)    +   _zCenter;

	_pAtom->setCoor(x,y,z);	
}

///********************************************
// *
// *  CRICKS/GEVORG PARAMETERIZATION
// *
// ********************************************/
//void CoiledCoils::gevorgCoiledCoil(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes){
// 	atoms.clear();
//
//	/***************************
//	 * Calc North's parameters 
//	 ***************************/
//	double pitch = (-1.0)*((2*M_PI*_r0)/tan(_a));
//	double rpr = (((_w0)*_r0)/sin(_a)); //Rise per Residue
//
//	double L = sqrt(pow(_r0,2)+pow((pitch /(2*M_PI)),2)); // the length of the spiral
//	double wa = -1.0*(_w1)*L / rpr;
//	double zCenter = 0;
//
//	
//	for (unsigned int i=0; i<_nRes; i++){
//		
//		/***************************
//		 * Nitrogen placement
//		 ***************************/
//		Atom *N = new Atom();
//		N->setResidueName("GLY");
//		N->setName("N");
//		N->setElement("N");
//		N->setResidueNumber(i+1);
//		N->setChainId("A");
//
//		double corrR1 = _r1 - 0.75;
//		double corrP1 = _phi1 + 0.555;
//		double corr = -1.015;
//		
//		generateGevorgCoor(N, i, _r0, _w0, _a, corrR1, _w1, corrP1, _dZ, corr, wa, zCenter);
//
//		/***************************
//		 * C-alpha placement
//		 ***************************/
//
//		Atom *CA = new Atom();
//		CA->setResidueName("GLY");
//		CA->setName("CA");
//		CA->setElement("C");
//		CA->setResidueNumber(i+1);
//		CA->setChainId("A");
//
//		corr = 0.0;
//
//		generateGevorgCoor(CA, i, _r0, _w0, _a, _r1, _w1, _phi1, _dZ, corr, wa, zCenter); 
//
//		/***************************
//		 * Carbonyl Carbon placement
//		 ***************************/
//		Atom *C = new Atom();
//		C->setResidueName("GLY");
//		C->setName("C");
//		C->setElement("C");
//		C->setResidueNumber(i+1);
//		C->setChainId("A");
//
//		corrR1 = _r1 - 0.54;
//		corrP1 = _phi1 - 0.71;
//		corr = 1.19;
//		
//		generateGevorgCoor(C, i, _r0, _w0, _a, corrR1, _w1, corrP1, _dZ, corr, wa, zCenter); 
//
//		/***************************
//		 * Carbonyl Oxygen placement
//		 ***************************/
//		Atom *O = new Atom();
//		O->setResidueName("GLY");
//		O->setName("O");
//		O->setElement("O");
//		O->setResidueNumber(i+1);
//		O->setChainId("A");
//
//		corrR1 = _r1 - 0.18;
//		corrP1 = _phi1 - 2.092;
//		corr = 2.512;
//
//		generateGevorgCoor(O, i, _r0, _w0, _a, corrR1, _w1, corrP1, _dZ, corr, wa, zCenter); 
//
//		/***************************
//		 * Add Atoms to atoms 
//		 ***************************/
//		atoms.push_back(N);
//		atoms.push_back(CA);
//		atoms.push_back(C);
//		atoms.push_back(O);
//
//		CA = N = C = O = NULL;
//	}
//}

//void CoiledCoils::generateGevorgCoor(Atom * _pAtom, unsigned int i, double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, double _corr, double _wa, double _zCenter) {
//
//	//double phi0prime = _dZ * tan(_a) / _r0;
//	double c1 = _corr;
//	double c2 = c1 / _wa;
//
//	double x = _r0*cos(_w0*i + c2)       +   _r1*cos(_w0*i + c2)*cos(_w1*i + c1 + _phi1)      -   _r1*cos(_a)*sin(_w0*i + c2)*sin(_w1*i+ c1 + _phi1);
//	double y = _r0*sin(_w0*i + c2)       +   _r1*sin(_w0*i + c2)*cos(_w1*i + c1 + _phi1)      +   _r1*cos(_a)*cos(_w0*i + c2)*sin(_w1*i+ c1 + _phi1);
//	double z = _r0*(_w0*i + c2)/tan(_a)  -   _r1*sin(_a)*sin(_w1*i + c1 + _phi1)    +   _zCenter;
//
//	_pAtom->setCoor(x,y,z);
//}

//void CoiledCoils::applyCoiledCoil(AtomPointerVector &_av, double _p0 ){
//
//	atoms.clear();
//	for (uint i = 0; i < _av.size(); i++){
//
//		double theta = -2 * M_PI * _av(i).getZ() / _p0;
//
//		double x = _av(i).getX() * cos(theta) - _av(i).getY() * sin(theta);
//		double y = _av(i).getX() * sin(theta) + _av(i).getY() * cos(theta);
//		
//		Atom *a = new Atom(_av(i));
//		a->setCoor(x,y,_av(i).getZ());
//		atoms.push_back(a);
//		
//	}
//
//}
//
//void CoiledCoils::useBBQTable(string _bbqTable){
//	bbqTable.openReader(_bbqTable);
//}
//
//Chain & CoiledCoils::getChain(){
//
//
//	if (sys != NULL){
//		delete(sys);
//		sys = new System();
//	}
//	sys->addAtoms(atoms);
//	bbqTable.fillInMissingBBAtoms((*sys)("A"));
//	return (*sys)("A");
//}
//
//System & CoiledCoils::getSystem(){
//
//	if (sys != NULL){
//		delete(sys);
//		sys = new System();
//	}
//	sys->addAtoms(atoms);
//	bbqTable.fillInMissingBBAtoms((*sys)("A"));
//	return (*sys);
//}

//void CoiledCoils::offersCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1){
//
//
//	atoms.clear();
//
//	double alpha = atan(2*M_PI*_r0 / _p0);
//	double w0    = 2*M_PI / (sqrt(pow(_p0,2)+ pow((2*M_PI*_r0),2)));
//
//
//	for (unsigned int i=1;i<=_nRes; i++) {
//
//		double t = (i - 1)*_risePerRes;
//		double theta = _w1 + 4*M_PI * (i-1) / 7;
//
//		Atom *CA = new Atom();
//		CA->setResidueName("ALA");
//		CA->setName("CA");
//		CA->setElement("C");
//		CA->setResidueNumber(i+1);
//
//
//		double x = _r0   * cos(w0 * t)  + _r1 * cos(w0*t)  * cos(theta)  + _r1 * cos(alpha)*sin(w0*t)*sin(theta);
//		double y = -(_r0 * sin(w0 * t ) + _r1 * sin(w0*t)  * cos(theta)) + _r1 * cos(alpha)*cos(w0*t)*sin(theta);
//		double z = _p0*w0*t/(2*M_PI)    + _r1 * sin(alpha) * sin(theta);
//
//		CA->setCoor(x,y,z);
//
//		atoms.push_back(CA);
//		CA  = NULL;
//	}
//	
//}

//void CoiledCoils::sotoCoiledCoils(double _r0, double _risePerRes, double _r1, int _nRes, double _resPerTurn, double _alpha, double _helicalPhase){
//	double OmegaSuper;
//	double OmegaAlpha;
//	double Temp[3];
//	double x,y,z;
//	double domega,dr,dz;
//	double dt,Pitch=0.00;
//	int t,i;
//
//	bool POSITIVE=false;
//	/*************************************************************
//	 * Be certain to decide whether you are going 
//	 * to be generating a left handed (negative) or right handed coil    
//	 *************************************************************/  
//	if (POSITIVE){
//		_alpha = _alpha*M_PI/180; 
//		Pitch = (2*M_PI*_r0)/(tan(_alpha));
//		OmegaSuper = (2*M_PI)/(Pitch/_risePerRes);
//		OmegaAlpha = (2*M_PI)/(_resPerTurn);
//		_helicalPhase *=M_PI/180;
//	}else{
//		_alpha = -1.0*_alpha*M_PI/180;
//		Pitch = (2*M_PI*_r0)/(tan(_alpha));
//		OmegaSuper = (2*M_PI)/(((2*M_PI*_r0)/(tan(_alpha)))/_risePerRes);
//		OmegaAlpha = (2*M_PI)/(_resPerTurn);
//		_helicalPhase *=M_PI/180;
//	}
//
//
//	atoms.clear();
//	Atom *N, *CA, *CB, *C, *O;
//	N = CA = CB = C = O = NULL;
//	for(t=-_nRes/2,i=0;t<_nRes/2+(_nRes%2);t++,i=i+12){
//
//		/***************************
//		 * For the nitrogen atom   *
//		 * domega = -26.89*PI/180  *
//		 * dr = -0.77              *
//		 * dz = -0.90              *
//		 ***************************/
//		domega = -26.89*M_PI/180;
//		dr     =  0.77;
//		dz     = -0.90;
//
//		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha-domega));
//		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha-domega)-_helicalPhase) + (dz)*sin(_alpha);
//		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha-domega)-_helicalPhase) - (dz)*cos(_alpha);
//
//		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
//		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
//		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);
//
//		N = new Atom();
//		N->setResidueName("ALA");
//		N->setName("N");
//		N->setElement("N");
//		N->setResidueNumber(i+1);
//		N->setChainId("A");
//		N->setCoor(x,y,z);
//
//		/**********
//		 *   CA   *
//		 **********/
//		x =         _r0*cos(t*OmegaSuper) + _r1*cos(t*OmegaSuper)*cos(t*OmegaAlpha-_helicalPhase) - _r1*sin(t*OmegaSuper)*cos(_alpha)*sin(t*OmegaAlpha-_helicalPhase);
//		y =         _r0*sin(OmegaSuper*t) + _r1*cos(t*OmegaAlpha-_helicalPhase)*sin(OmegaSuper*t) +  _r1*sin(t*OmegaAlpha-_helicalPhase)*cos(t*OmegaSuper)*cos(_alpha);
//		z =    -1.0*_r1*sin(_alpha)*sin(t*OmegaAlpha-_helicalPhase) +  (_r0*OmegaSuper*t)/tan(_alpha); 
//
//		CA = new Atom();
//		CA->setResidueName("ALA");
//		CA->setName("CA");
//		CA->setElement("C");
//		CA->setResidueNumber(i+1);
//		CA->setChainId("A");
//		CA->setCoor(x,y,z);
//
//
//		/*********************************
//		 * For the carboxyl carbon       *
//		 * delta omega = 27.42*PI/180    *
//		 * delta r = -0.55               *
//		 * delta z = 1.05                *
//		 *********************************/
//		domega = 27.42*M_PI/180;
//		dr     = 0.55;
//		dz     = 1.05;
//
//		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha+domega));
//		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) + (dz)*sin(_alpha);
//		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) - (dz)*cos(_alpha);
//	  
//		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
//		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
//		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);
//
//		C = new Atom();
//		C->setResidueName("ALA");
//		C->setName("C");
//		C->setElement("C");
//		C->setResidueNumber(i+1);
//		C->setChainId("A");
//
//		C->setCoor(x,y,z);
//
//		/**********************************
//		 * For the oxygen atom            *
//		 * delta omega = 23.85*PI/180     *
//		 * delta r = -0.19                *
//		 * delta z = 2.23                 *
//		 **********************************/
//		domega = 23.85*M_PI/180;
//		dr     = 0.19;
//		dz     = 2.23;
//
//		     
//		Temp[0] =      (_r1-dr)*cos(_helicalPhase -     (t*OmegaAlpha+domega));
//		Temp[1] =      (_r1-dr)*cos(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) + (dz)*sin(_alpha);
//		Temp[2] =      (_r1-dr)*sin(_alpha)*sin((t*OmegaAlpha+domega)-_helicalPhase) - (dz)*cos(_alpha);
//
//		x = cos(OmegaSuper*t)*Temp[0] - sin(OmegaSuper*t)*Temp[1] + _r0*cos(OmegaSuper*t);
//		y = sin(OmegaSuper*t)*Temp[0] + cos(OmegaSuper*t)*Temp[1] + _r0*sin(OmegaSuper*t);
//		z = -1.0*Temp[2] + (_r0*OmegaSuper*t)/tan(_alpha);
//
//
//		O = new Atom();
//		O->setResidueName("ALA");
//		O->setName("O");
//		O->setElement("O");
//		O->setResidueNumber(i+1);
//		O->setChainId("A");
//		O->setCoor(x,y,z);
//
//
//		// CB Build from N,CA,C atoms
//		CB = new Atom(*CA);
//		CB->setResidueName("ALA");
//		CB->setName("CB");
//		CB->setElement("C");
//		CB->setChainId("A");
//		CB->setResidueNumber(i+1);
//		CB->setCoor(CartesianGeometry::build(CA->getCoor(), N->getCoor(), C->getCoor(), 1.521, 110.5, -122.5));			
//
//
//		atoms.push_back(N);
//		atoms.push_back(CA);
//		atoms.push_back(CB);
//		atoms.push_back(C);
//		atoms.push_back(O);
//
//		N = CA = CB = C = O = NULL;
//  
//	}
//
//}
