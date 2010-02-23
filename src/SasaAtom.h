#ifndef SASAATOM_H
#define SASAATOM_H

#include<Atom.h>
#include<SurfaceSphere.h>

#define INVALID_SASA -1.0

/******************************************************************
 *  A container for calculating the SASA of an atom
 ******************************************************************/

namespace MSL { 
class SasaAtom {
	public:
		SasaAtom();
		SasaAtom(Atom  * _atom, double _probeRadius=0.0, unsigned int _numberOfPoints=0);
		SasaAtom(CartesianPoint & _center, double _radius, double _probeRadius=0.0, unsigned int _numberOfPoints=0);
		SasaAtom(const SasaAtom & _sasaAtom);
		~SasaAtom();
	

		SurfaceSphere* & getSurfaceSphere();
		
		void deleteOcclusionPoints();
		double calcSasa();
		double getSasa() const;
		//void setSasa(double _s);
		CartesianPoint * getCenter() const;

		bool buildOcclusionPoints(unsigned int _n);

		unsigned int getSurfaceSphereCurrentSize() const; // the current number of points
		unsigned int getSurfaceSphereOiginalSize() const; // the size of the original sphere
		double getSurfaceFraction() const; // return the ration between the current size and the original number of points

		double getProbeRadius() const;
		void setProbeRadius(double _probeRadius);
		double getRadius() const;
		void setRadius(double _radius);

		void removeOccluded(Atom & _atom, bool _purge=false);
		void removeOccluded(CartesianPoint & _center, double _radius, bool _purge=false);
		void removeOccluded(SasaAtom & _sasaAtom, bool _purge=false);

		double distance(const SasaAtom & _sasaAtom) const;

	private:
		void setup(CartesianPoint * _center, double _radius, double _probeRadius, unsigned int _numberOfPoints);

		//Atom *atom;
		CartesianPoint * pCenter;
		SurfaceSphere *ss;
		double sasa;
		double probeRadius;
		double radius;
		double totalRadius;
		unsigned int numberOfPoints;
};

inline SasaAtom::~SasaAtom() { delete ss;};
inline SurfaceSphere* & SasaAtom::getSurfaceSphere() {
	return ss;
}
inline double SasaAtom::calcSasa() {
	if (numberOfPoints == 0) {
		sasa = 0.0;
		return 0.0;
	}
	sasa = ss->getSurfaceFraction() * 4.0 * M_PI * pow(totalRadius,2.0);
	return sasa;
}

inline double SasaAtom::getSasa() const { return sasa;}
inline CartesianPoint * SasaAtom::getCenter() const { return pCenter;}
inline void SasaAtom::deleteOcclusionPoints() { ss->deleteAllPoints();}
inline bool SasaAtom::buildOcclusionPoints(unsigned int _n) { numberOfPoints = _n; return ss->buildOcclusionPoints(*pCenter, totalRadius, numberOfPoints);}
inline double SasaAtom::getProbeRadius() const { return probeRadius;}
inline void SasaAtom::setProbeRadius(double _probeRadius) { probeRadius = _probeRadius; totalRadius = probeRadius + radius;}
inline double SasaAtom::getRadius() const { return radius;}
inline void SasaAtom::setRadius(double _radius) { radius = _radius; totalRadius = probeRadius + radius;}
inline unsigned int SasaAtom::getSurfaceSphereCurrentSize() const {return ss->size();}
inline unsigned int SasaAtom::getSurfaceSphereOiginalSize() const {return ss->originalSize();}
inline double SasaAtom::getSurfaceFraction() const {return ss->getSurfaceFraction();}
inline double SasaAtom::distance(const SasaAtom & _sasaAtom) const {return pCenter->distance(*_sasaAtom.pCenter);}

}

#endif
