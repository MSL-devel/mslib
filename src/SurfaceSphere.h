#ifndef SURFACE_SPHERE_H
#define SURFACE_SPHERE_H

#include <math.h>
#include "Real.h"
#include "CartesianPoint.h"


namespace MSL { 
class SurfaceSphere {

	public:
		SurfaceSphere();
		SurfaceSphere (SurfaceSphere & _ss ); 
		~SurfaceSphere();
		
		void operator=(const SurfaceSphere &_ss);
		unsigned int size() const; // the current number of points
		unsigned int originalSize() const; // the size of the original sphere
		double getSurfaceFraction() const; // return the ration between the current size and the original number of points

		CartesianPoint & operator[](size_t _n);
		bool buildOcclusionPoints(CartesianPoint & _pt, double _rad, unsigned int _n);

		void setOccluded(unsigned int _n);
		bool getOccluded(unsigned int _n) const;
		void purgeOccluded();

		void deleteAllPoints();
	
		int getNumberofOccludedPoints();

		std::vector<CartesianPoint*> & getOcclusionPoints() ;
		

	private:

		void deletePointers();
		std::vector<CartesianPoint*> points;
		std::vector<bool> occluded;
		unsigned int originalNumberOfPoints;

};

inline void SurfaceSphere::setOccluded(unsigned int _n) {
	occluded[_n] = true;
}
inline bool SurfaceSphere::getOccluded(unsigned int _n) const {return occluded[_n];}
		
inline CartesianPoint & SurfaceSphere::operator[](size_t _n) {
	return *points[_n];
}
inline unsigned int SurfaceSphere::size() const {return points.size();}
inline unsigned int SurfaceSphere::originalSize() const {return originalNumberOfPoints;}
inline double SurfaceSphere::getSurfaceFraction() const {return (double)points.size()/(double)originalNumberOfPoints;}

inline int SurfaceSphere::getNumberofOccludedPoints() {
	int ret = 0;
	for (unsigned int i = 0; i < occluded.size(); i++) {
		if(occluded[i]) {
			ret++;	
		}
	}
	return ret;
}

	
inline std::vector<CartesianPoint*> & SurfaceSphere::getOcclusionPoints () {
	return points;
}
	
}

#endif
