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
#ifndef SASACALCULATOR_H
#define SASACALCULATOR_H

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

/*********************************************************
  TO BE IMPLEMENTED:
  - external file reader for radii
  - return the sasa points that have not been occluded
  - support for atom, residue, chain SASA
**********************************************************/

namespace MSL { 
class SasaCalculator {
	public:
		SasaCalculator();
		SasaCalculator(SasaCalculator &_sasaCalculator);
		SasaCalculator(AtomPointerVector& _atoms, double _probeRadius=DEFAULT_PROBE_RADIUS, const int _no_OcclusionPoints = DEFAULT_OCCLUSION_POINTS);
		~SasaCalculator() { deletePointers();};

		void addAtoms(AtomPointerVector& _atoms);
		std::vector<SasaAtom*> & getAtomPointers();
		void calcSasa();
	//	double getAtomSasa(std::string _atomId); // use "A 7 CA" or "A,7,CA"
		double getResidueSasa(std::string _positionId); // use "A 7" or "A,7"
		std::string getSasaTable(bool _byAtom=true); // if _byAtom == false print residue sasa
		void printSasaTable(bool _byAtom=true);
		std::string getResidueSasaTable();
		void printResidueSasaTable();
		double getTotalSasa();
//		vector<CartesianPoint>  getOcclusionPoints();
//		void printOcclusionPoints();
//		void printCubes();
		void deletePointers();
		void setProbeRadius (double _probeRadius);
		double getProbeRadius() const;

		void setUseDefaultRadii(bool _flag); // if true, zero radii are replaced by element default radii
		bool getUseDefaultRadii() const;

		std::map <std::string,double> getRadiiMap() const;
		void setRadiiMap(const std::map <std::string,double> & _radiiMap);

		void setTempFactorWithSasa(bool _flag); // if true saves the sasa also on the B factor
		bool getTempFactorWithSasa() const;

//		bool readRadiiMap(std::string _mapFile); TO BE IMPLEMENTED!!!

	private:
		void setup(int _noOcclusionPoints, double _probeRadius);
		//std::vector<std::vector<SasaAtom*> > findNeighbors();
		std::vector<AtomPointerVector> findNeighbors();

		std::map <std::string,double> atomRadii;
	//	std::map <std::string, std::map<std::string, std::map<std::string, double> > > AtomSasa;
	//	std::map <std::string, std::map<std::string, double> > ResidueSasa;
	//	std::map <std::string, double> ChainSasa;

		//std::map <std::string ,std::map <int, double> > residueSasaMap;
		AtomPointerVector atoms;
		std::vector<SasaAtom*> sasaAtoms;
		int noOcclusionPoints;
		bool useDefaultRadii;
		double probeRadius;

		bool setTFactor;
};


inline void SasaCalculator::setUseDefaultRadii(bool _flag) { useDefaultRadii = _flag;}
inline bool SasaCalculator::getUseDefaultRadii() const {return useDefaultRadii;}
inline void SasaCalculator::setProbeRadius (double _probeRadius) {probeRadius = _probeRadius;}
inline double SasaCalculator::getProbeRadius() const {return probeRadius;}
inline std::map <std::string,double> SasaCalculator::getRadiiMap() const {return atomRadii;}
inline void SasaCalculator::setRadiiMap(const std::map <std::string,double> & _radiiMap) {atomRadii = _radiiMap;}
inline void SasaCalculator::setTempFactorWithSasa(bool _flag) {setTFactor = _flag;}
inline bool SasaCalculator::getTempFactorWithSasa() const {return setTFactor;}
inline std::string SasaCalculator::getResidueSasaTable() {return getSasaTable(false);}
inline void SasaCalculator::printResidueSasaTable() {printSasaTable(false);}
//inline double getAtomSasa(std::string _chain_resnumr_name) {
//	if (atoms->exist(_chain_resnumr_name)) {
//		return atoms->getLastFoundAtom
//inline double getResidueSasa(std::string _positionId); // use "A 7" or "A,7"

//inline bool SasaCalculator::readRadiiMap(std::string _mapFile) {
//	std::ifstream file_fs;
//	file_fs.open(_mapFile.c_str());
//	if (file_fs.fail()) {
//		return false;
//	}
//	
//}


}

#endif
