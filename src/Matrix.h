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


#ifndef MATRIX_H
#define MATRIX_H

// MSL Includes
#include <string>
#include <vector>
#include <iostream>
#include <complex>
#include <fstream>
#include <cstdlib>
#include <sys/types.h>

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#endif

#ifdef __GSL__
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#endif


// Forward Declarations
namespace MSL { 


// Namespaces


class Matrix {

	public:
		Matrix();
		Matrix(unsigned int _rows, unsigned int _cols);
		Matrix(unsigned int _rows, unsigned int _cols, double _val); // fill with this value
		Matrix(std::vector<std::vector<double> > _matrixValues); // fill with these values
		Matrix(const Matrix  & _m); // copy constructor

		~Matrix();

		std::vector<double> & operator[](size_t _n);
		double getElement(unsigned int _row, unsigned int _col) const;
		unsigned int getRows() const {return rows;};
		unsigned int getCols() const {return cols;};
		std::vector<double> getRow(unsigned int _row) const;

		Matrix operator*(const Matrix & _m) const;
		void operator*=(const Matrix & _m);
		friend std::ostream & operator<<(std::ostream &_os, const Matrix & _mtrx) {_os << _mtrx.toString(); return _os;};


		double getDeterminant() const;
		Matrix getMinor(unsigned int _row, unsigned int _col) const;
		Matrix getSubMatrix(unsigned int _rowStart, unsigned int _rowEnd, unsigned int _colStart, unsigned int _colEnd) const;
		Matrix getTranspose() const;
		std::vector<std::vector<double> > getEigenvectorsGSL(); // returns 4x4
		std::vector<std::vector<double> > getEigenvectorsEigenValuesGSL(); // returns 4x5

		void addRow(std::vector<double> _vals);
		void addCol(std::vector<double> _vals);

		void initialize();
		void initialize(unsigned int _rows, unsigned int _cols);
		void initialize(unsigned int _rows, unsigned int _cols, double _val);


		std::string toString() const;


	protected:
	
		std::vector<std::vector<double> > matrix;
		unsigned int rows;
		unsigned int cols;
		std::string archiveType;

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(std::string filename) const{

			if (archiveType == "binary"){
				std::ofstream fout(filename.c_str(),std::ios::binary);
				boost::archive::binary_oarchive oa(fout);
				oa << (*this);
			} else if (archiveType == "xml"){
				std::ofstream fout(filename.c_str());
				boost::archive::xml_oarchive oa(fout);
				oa << boost::serialization::make_nvp("Matrix",*this);
			} else {
				std::ofstream fout(filename.c_str());
				boost::archive::text_oarchive oa(fout);
				oa << (*this);
			}

		}

		void load_checkpoint(std::string filename){

			if (archiveType == "binary"){
				std::ifstream fin(filename.c_str(), std::ios::binary);
				boost::archive::binary_iarchive ia(fin);
				ia >> (*this);
			} else if (archiveType == "xml"){
				std::ifstream fin(filename.c_str());
				boost::archive::xml_iarchive ia(fin);
				ia >> boost::serialization::make_nvp("Matrix",*this);
			} else {
				std::ifstream fin(filename.c_str());
				boost::archive::text_iarchive ia(fin);
				ia >> (*this);
			}
		}


	private:
		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			using boost::serialization::make_nvp;

			ar & make_nvp("matrix",matrix);
			ar & make_nvp("rows",rows);
			ar & make_nvp("cols",cols);
		}
#else
	public:
		void save_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF SAVE_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
		void load_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF LOAD_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}		
#endif

};

}

#endif
