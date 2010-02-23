#include<SurfaceSphere.h>

using namespace MSL;
using namespace std;


		
SurfaceSphere::SurfaceSphere() {
	originalNumberOfPoints = 0;
}

SurfaceSphere::SurfaceSphere(SurfaceSphere & _ss ) {
	deletePointers();
	for (int i = 0; i < _ss.points.size(); i++) {
		points.push_back(new CartesianPoint(*(_ss.points[i])));
		occluded.push_back(_ss.occluded[i]);
	}
	originalNumberOfPoints = _ss.originalNumberOfPoints;

}
 
void SurfaceSphere::operator=(const SurfaceSphere & _ss ) {
	deletePointers();
	for (int i = 0; i < _ss.points.size(); i++) {
		points.push_back(new CartesianPoint(*(_ss.points[i])));
		occluded.push_back(_ss.occluded[i]);
	}
	originalNumberOfPoints = _ss.originalNumberOfPoints;

}
SurfaceSphere::~SurfaceSphere () { 
	deletePointers();
}

void SurfaceSphere::deleteAllPoints() {
	deletePointers();
	originalNumberOfPoints = 0;
}

void SurfaceSphere::deletePointers() {
	for(int i = 0; i < points.size(); i++) {
		delete points[i];	
	}
	points.clear();
	occluded.clear();
}


 	
bool SurfaceSphere::buildOcclusionPoints(CartesianPoint & _pt, double _rad, unsigned int _n) {
	/**************************************************************
	 *   http://cgafaq.info/wiki/Evenly_distributed_points_on_sphere
	 *
	 *   [...] A related method chooses successive longitudes according to the 
	 *   "most irrational number"  (known as the golden section) so that no two 
	 *   nodes in nearby bands come too near each other in longitude. This method 
	 *   obtains a slightly better packing than the above:
	 *
	 *        dlong := pi*(3-sqrt(5))
	 *        long := 0
	 *        dz := 2.0/N
	 *        z := 1 - dz/2
	 *        for k := 0 .. N-1
	 *            r := sqrt(1-z*z)
	 *            pt[k] := (cos(long)*r, sin(long)*r, z)
	 *            z := z - dz
	 *            long := long + dlong
	 **************************************************************/
	originalNumberOfPoints = _n;
	deletePointers();
	if (_n == 0) {
		return true;
	}
	double dlongitude = M_PI*(3-sqrt(5));
	double longitude = 0.0;
	double dz = 2.0/(double)_n;
	double z = 1.0 - dz/2.0;
	for (int i=0; i < _n; i++) {
		double r = sqrt(1-z*z);
		CartesianPoint p;
		p.setX(cos(longitude)*r);
		p.setY(sin(longitude)*r);
		p.setZ(z);
		p = p*_rad + _pt;
		points.push_back(new CartesianPoint(p));
		occluded.push_back(false);
		z = z-dz;
		longitude = longitude + dlongitude;
	}
	return true;
}

void SurfaceSphere::purgeOccluded() {
	unsigned int counter = 0;
	vector<CartesianPoint*> tmp;
	tmp.reserve(points.size());
	for (vector<CartesianPoint*>::iterator k=points.begin(); k!=points.end(); k++) {
		if (occluded[counter]) {
			delete *k;
			*k = NULL;
		} else {
			tmp.push_back(*k);
		}
		counter++;
	}
	points = tmp;
	occluded = vector<bool>(points.size(), false);
}

