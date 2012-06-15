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

#ifndef PRINICPLECOMPONENTANALYSIS_H
#define PRINICPLECOMPONENTANALYSIS_H

#include "AtomPointerVector.h"
#include "Line.h"
#include "CartesianPoint.h"
#include "Transforms.h"




// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>
#endif

namespace MSL { 
class PrincipleComponentAnalysis {
	public:
	       PrincipleComponentAnalysis();	
	       PrincipleComponentAnalysis(const PrincipleComponentAnalysis & _pca);
	       ~PrincipleComponentAnalysis();	

	       void operator=(const PrincipleComponentAnalysis & _pca); // assignment

	       void computePrincipleComponents(AtomPointerVector &_av);

	       void printPrincpleComponents();

	       // Set/Get Functions
	       std::vector<Line> getLines();
	       std::vector<std::vector<double> > getEigenVectors() { return eigenvectors; }
	       std::vector<double> getEigenValues() { return eigenvalues; }  
	       void setPrintImaginary(bool _flag) { printImag = _flag; }

	private:	

	       void copy(const PrincipleComponentAnalysis & _pca);

	       Matrix compMat;
	       std::vector<std::vector<double> > eigenvectors;
	       std::vector<double> eigenvalues;
	       bool printImag;

	       CartesianPoint toGeoCenter;
	       CartesianPoint toOrigin;   

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(std::string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(std::string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & compMat;
			ar & eigenvectors;
			ar & eigenvalues;
			ar & printImag;
			ar & toGeoCenter;
			ar & toOrigin;
		}
#endif

};
}

#endif
