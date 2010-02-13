#ifndef SASABUILDER_H
#define SASABUILDER_H

#include <sstream>

#include "AtomPointerVector.h"
#include "SurfaceSphere.h"
#include "SasaAtom.h"
#include "SurfaceSphere.h"
#include "Atom3DGrid.h"

//#include "math.h"
//#include "Real.h"

#define DEFAULT_OCCLUSION_POINTS (2000)
#define DEFAULT_MIN_ATOMSIZE_FOR_ATOMGRID (2000)

#define DEFAULT_PROBE_RADIUS 1.4
#define INVALID_SASA -1.0

using namespace std;
/*********************************************************
  TO BE IMPLEMENTED:
  - external file reader for radii
  - return the sasa points that have not been occluded
**********************************************************/

class SasaCalculator {
	public:
		SasaCalculator();
		SasaCalculator(SasaCalculator &_sasaCalculator);
		SasaCalculator(AtomPointerVector& _atoms, double _probeRadius=DEFAULT_PROBE_RADIUS, const int _no_OcclusionPoints = DEFAULT_OCCLUSION_POINTS);
		~SasaCalculator() { deletePointers();};

		void addAtoms(AtomPointerVector& _atoms);
		vector<SasaAtom*> & getAtoms();
		void calcSasa();
		double getResidueSasa(string _chainId,int _resNumber);
		string getSasaTable(bool _byAtom=true); // if _byAtom == false print residue sasa
		void printSasaTable(bool _byAtom=true);
//		vector<CartesianPoint>  getOcclusionPoints();
//		void printOcclusionPoints();
//		void printCubes();
		void deletePointers();
		void setProbeRadius (double _probeRadius);
		double getProbeRadius() const;

		void setUseDefaultRadii(bool _flag); // if true, zero radii are replaced by element default radii
		bool getUseDefaultRadii() const;

		map <string,double> getRadiiMap() const;
		void setRadiiMap(const map <string,double> & _radiiMap);

		void setTempFactorWithSasa(bool _flag); // if true saves the sasa also on the B factor
		bool getTempFactorWithSasa() const;

//		bool readRadiiMap(string _mapFile); TO BE IMPLEMENTED!!!

	private:
		void setup(int _noOcclusionPoints, double _probeRadius);
		//vector<vector<SasaAtom*> > findNeighbors();
		vector<AtomPointerVector> findNeighbors();

		map <string,double> atomRadii;

		//map <string ,map <int, double> > residueSasaMap;
		AtomPointerVector atoms;
		vector<SasaAtom*> sasaAtoms;
		int noOcclusionPoints;
		bool useDefaultRadii;
		double probeRadius;

		bool setTFactor;
};


inline void SasaCalculator::setUseDefaultRadii(bool _flag) { useDefaultRadii = _flag;}
inline bool SasaCalculator::getUseDefaultRadii() const {return useDefaultRadii;}
inline void SasaCalculator::setProbeRadius (double _probeRadius) {probeRadius = _probeRadius;}
inline double SasaCalculator::getProbeRadius() const {return probeRadius;}
inline map <string,double> SasaCalculator::getRadiiMap() const {return atomRadii;}
inline void SasaCalculator::setRadiiMap(const map <string,double> & _radiiMap) {atomRadii = _radiiMap;}
inline void SasaCalculator::setTempFactorWithSasa(bool _flag) {setTFactor = _flag;}
inline bool SasaCalculator::getTempFactorWithSasa() const {return setTFactor;}

//inline bool SasaCalculator::readRadiiMap(string _mapFile) {
//	ifstream file_fs;
//	file_fs.open(_mapFile.c_str());
//	if (file_fs.fail()) {
//		return false;
//	}
//	
//}


#endif
