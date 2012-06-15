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

#include "CoiledCoils.h"

using namespace MSL;
using namespace std;


CoiledCoils::CoiledCoils(){
	//sys = new System();
	CAname = "CA";
	Nname = "N";
	Cname = "C";
	Oname = "O";

}

CoiledCoils::~CoiledCoils(){
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
//	if (sys != NULL){
//		delete(sys);
//	}
}

bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions) {
	cerr << "WARNING deprecated function bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions), use setSystemToCoiledCoil instead" << endl;
	return setSystemToCoiledCoil(_r0, _risePerRes, _pitch, _r1, _w1, _phi1, _dZ, _nRes, _symmetry, _N, _sys, _startingPositions);
}

/********************************************
 *
 *  SYSTEM TO COILED COIL BUNDLE 
 *
 ********************************************/
bool CoiledCoils::setSystemToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions) {
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
		cerr << "WARNING No. of starting positions does not match no. of bundles in bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions)" << endl;
		return false;
	}

	//Test to see if all starting positions are contained within the system
	for (int i = 0; i < _startingPositions.size(); i++){
		if (_sys.positionExists(_startingPositions[i])) {
			startPositionIndices.push_back(_sys.getPositionIndex(&_sys.getLastFoundPosition()));
		} else {
			cerr << "WARNING Not all starting positions are contained within the bundle in bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions)" << endl;
			return false;
		}
	}

	//Test to see if starting positions overlap with one another 
	for (int i = 0; i < startPositionIndices.size()-1; i++) {
		if (abs(startPositionIndices[i] - startPositionIndices[i+1]) > _nRes) {
			cerr << "WARNING Starting positions overlap with one another in bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions)" << endl;
			return false;
		}
	}
	if ((startPositionIndices[startPositionIndices.size()-1] + _nRes) > _sys.positionSize()) {
		cerr << "WARNING Starting position overlaps with end of the sequence in bool CoiledCoils::primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes, string _symmetry, int _N, System& _sys, vector<string> _startingPositions)" << endl;
		return false;
	}

	//Create the coiled coil bundle in the atomPointerVector atoms
	getCoiledCoilBundle(_r0, _risePerRes, _pitch, _r1, _w1, _phi1, _dZ, _nRes, _symmetry, _N);

	//Divide out the APV into the chainVector
	vector< vector< map<string, Atom*> > > chainVector = atoms.subdivideByChainAndPosition();

	// apply the coordinates of the backbone atoms to the system
	for (int k = 0; k < bundleSize; k++) {
		int resCount = 0;
		for (int j = startPositionIndices[k]; j < (startPositionIndices[k]+_nRes); j++) {
			Residue * pRes = &(_sys.getIdentity(j));

			if (pRes->atomExists(Nname) && chainVector[k][resCount].find(Nname) != chainVector[k][resCount].end()) {
				pRes->getLastFoundAtom().setCoor(chainVector[k][resCount][Nname]->getCoor());
			}

			if (pRes->atomExists(CAname) && chainVector[k][resCount].find(CAname) != chainVector[k][resCount].end()) {
				pRes->getLastFoundAtom().setCoor(chainVector[k][resCount][CAname]->getCoor());
			}

			if (pRes->atomExists(Cname) && chainVector[k][resCount].find(Cname) != chainVector[k][resCount].end()) {
				pRes->getLastFoundAtom().setCoor(chainVector[k][resCount][Cname]->getCoor());
			}

			if (pRes->atomExists(Oname) && chainVector[k][resCount].find(Oname) != chainVector[k][resCount].end()) {
				pRes->getLastFoundAtom().setCoor(chainVector[k][resCount][Oname]->getCoor());
			}

			/*
			for (int i = 0; i < _sys.getPosition(j).atomSize(); i++) {
				string atomName = _sys.getPosition(j).getAtom(i).getName();
				if (atomName == "CA" || atomName == "C" || atomName == "O" || atomName == "N") {
					_sys.getPosition(j).getAtom(i).setCoor(chainVector[k][resCount][atomName]->getCoor());
				}
			}
			*/
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
		//N->setName("N");
		N->setName(Nname);
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
		//CA->setName("CA");
		CA->setName(CAname);
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
		//C->setName("C");
		C->setName(Cname);
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
		//O->setName("O");
		O->setName(Oname);
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

