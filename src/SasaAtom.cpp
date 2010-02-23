#include<SasaAtom.h>

using namespace MSL;


SasaAtom::SasaAtom() {
	/*
	atom  = NULL;
	ss = NULL; 
	sasa = INVALID_SASA;
	probeRadius = 0.0;
	totalRadius = 0.0;
	radius = 0.0;
	numberOfPoints = 0;
	*/
	setup(NULL, 0.0, 0.0, 0);
}

SasaAtom::SasaAtom(Atom  * _atom, double _probeRadius, unsigned int _numberOfPoints) {
	/*
	atom  = _atom;
	ss = new SurfaceSphere();
	sasa = INVALID_SASA;
	probeRadius = _probeRadius;
	radius = _atom->getRadius();
	totalRadius = radius + probeRadius;
	numberOfPoints = 0;
	*/
	setup(&_atom->getCoor(), _atom->getRadius(), _probeRadius, _numberOfPoints);
}

SasaAtom::SasaAtom(CartesianPoint & _center, double _radius, double _probeRadius, unsigned int _numberOfPoints) {
	setup(&_center, _radius, _probeRadius, _numberOfPoints);
}


SasaAtom::SasaAtom(const SasaAtom &  _sasaAtom) {
	//atom  = _sasaAtom.atom;
	pCenter  = _sasaAtom.pCenter;
	if(_sasaAtom.ss != NULL) {
		ss = new SurfaceSphere(*(_sasaAtom.ss));
	} else {
		ss = NULL;
	}
	sasa = _sasaAtom.sasa;
	probeRadius = _sasaAtom.probeRadius;
	radius = _sasaAtom.radius;
	totalRadius = _sasaAtom.totalRadius;
	numberOfPoints = _sasaAtom.numberOfPoints;
}


void SasaAtom::removeOccluded(Atom & _atom, bool _purge) {
	removeOccluded(_atom.getCoor(), _atom.getRadius(), _purge);
}

void SasaAtom::removeOccluded(SasaAtom & _sasaAtom, bool _purge) {
	removeOccluded(*_sasaAtom.getCenter(), _sasaAtom.getRadius(), _purge);
}

void SasaAtom::removeOccluded(CartesianPoint & _center, double _radius, bool _purge) {
	// remove all points that are occluded by the other atom	
	double effectiveRadius = _radius + probeRadius;
	unsigned int counter = 0;
	for (int j = 0; j < ss->size(); j++) {
		if (!ss->getOccluded(j) && (*ss)[j].distance(_center) <= effectiveRadius) {
			ss->setOccluded(j);
			counter++;
		}
	}
	if (_purge || counter > 100 || counter == ss->size()) {
		// purge when enough points have been flagged
		// this could be optimized more rationally with
		// a large number of test cases
		ss->purgeOccluded();
	}
}


void SasaAtom::setup(CartesianPoint * _center, double _radius, double _probeRadius, unsigned int _numberOfPoints) {
	pCenter = _center;
	probeRadius = _probeRadius;
	radius = _radius;
	totalRadius = radius + probeRadius;
	numberOfPoints = _numberOfPoints;
	sasa = 0.0;
	ss = new SurfaceSphere();
	if (numberOfPoints > 0) {
		ss->buildOcclusionPoints(*pCenter, totalRadius, numberOfPoints);
	}
}
