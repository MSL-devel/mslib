#include "AtomBondBuilder.h"

using namespace MSL;
using namespace std;


AtomBondBuilder::AtomBondBuilder() {
	setup();
}

AtomBondBuilder::AtomBondBuilder(const AtomBondBuilder & _abb) {
	setup();
	atomRadii = _abb.atomRadii;
}

AtomBondBuilder::~AtomBondBuilder() {
}

void AtomBondBuilder::setup() {
	useDefaultRadii = true;
	atomRadii["C"] = 2.00;
	atomRadii["H"] = 1.09;
	atomRadii["N"] = 1.85;
	atomRadii["O"] = 1.735;
	atomRadii["CA"] = 1.71;
	atomRadii["FE"] = 0.65;
	atomRadii["S"] = 2.06;
	atomRadii["ZN"] = 1.09;
	atomRadii["HE"] = 1.48;
	atomRadii["NE"] = 1.53;
	atomRadii["NA"] = 1.36375;
	atomRadii["K"] = 1.76375;
	atomRadii["CL"] = 2.27;
	atomRadii["MG"] = 1.185;
	atomRadii["CS"] = 2.1;
	atomRadii["CU"] = 1.38;
	atomRadii["P"] = 2.15;
	atomRadii["SE"] = 2.20;
	factor = 0.47;
}

void AtomBondBuilder::buildConnections(AtomPointerVector & _atoms) {

	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		for (AtomPointerVector::iterator l=k+1; l!=_atoms.end(); l++) {

			string atomk = (*k)->getElement();
			string atoml = (*l)->getElement();

			float kRad = 0;
			float lRad = 0;

			if (atomRadii.find(atomk) != atomRadii.end()) {
				kRad = atomRadii[atomk];
			}

			if (atomRadii.find(atomk) != atomRadii.end()) {
				lRad = atomRadii[atoml];
			}

			float dist = (*k)->distance(**l);
			float vdwrSum = kRad + lRad; 

			if (dist < (vdwrSum*factor)) {
		//		cout << (**k) << " to " << (**l) << " at: " << dist << endl;
				(*k)->setBoundTo(*l);
			}


		}
	}

	/*
	// Future testing: use short PDB file and CHARMM and compare the number of bonds present	
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		cout << (**k) << endl;
		vector<Atom*> bondMap = (*k)->getBonds();
		for (vector<Atom*>::iterator l=bondMap.begin(); l!=bondMap.end(); l++) {
			cout << "   " << **l << endl;
		}
		cout << endl;
	}
	*/


}
