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


#include "Matrix.h"
#include "CartesianGeometry.h"

#ifdef __GSL__
  #include <gsl/gsl_eigen.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_math.h>
  #include <gsl/gsl_complex_math.h>

using namespace MSL;
using namespace std;

#endif

using namespace MSL;
using namespace std;

Matrix::Matrix() {
	initialize(0, 0, 0.0);
}

Matrix::Matrix(unsigned int _rows, unsigned int _cols) {
	initialize(_rows, _cols, 0.0);
}

Matrix::Matrix(unsigned int _rows, unsigned int _cols, double _val) {
	initialize(_rows, _cols, _val);
}

Matrix::Matrix(vector<vector<double> > _matrixValues) {
	// make sure all the rows have the same number of columns,
	// or else pad with zeros
	rows = _matrixValues.size();
	cols = 0;
	for (int i=0; i<_matrixValues.size(); i++) {
		if (_matrixValues[i].size() > cols) {
			cols = _matrixValues[i].size();
		}
	}

	for (int i=0; i<_matrixValues.size(); i++) {
		while (_matrixValues[i].size() < cols) {
			_matrixValues[i].push_back(0.0);
		}
	}

	matrix = _matrixValues;
}

Matrix::Matrix(const Matrix  & _m) {
	matrix = _m.matrix;
	rows = _m.rows;
	cols = _m.cols;
}

Matrix::~Matrix() {
}

vector<double> & Matrix::operator[](size_t n) {
	if (n > rows){
		cerr << "ERROR 1942 Matix::operator[] n is larger than number of rows. n = "<<n<<"; rows = "<<rows<<endl;
		exit(1942);
	}
	return matrix[n];
}

void Matrix::initialize() {
	initialize(rows, cols, 0.0);
}

void Matrix::initialize(unsigned int _rows, unsigned int _cols) {
	initialize(_rows, _cols, 0.0);
}

void Matrix::initialize(unsigned int _rows, unsigned int _cols, double _val) {
	matrix.clear();
	rows = _rows;
	cols = _cols;
	for (unsigned int i=0; i<rows; i++) {
		matrix.push_back(vector<double>(cols, _val));
	}
	archiveType = "binary";
}

double Matrix::getElement(unsigned int _row, unsigned int _col) const {
	return matrix[_row][_col];
}

vector<double> Matrix::getRow(unsigned int _row) const {
	return matrix[_row];
}

Matrix Matrix::operator*(const Matrix & _m) const {
	if (cols != _m.rows) {
		cerr << "ERROR 9232: incorrect size for multiply matrices ("<< rows << ", " << cols << ") ("<< _m.rows << ", " << _m.cols << ") in Matrix Matrix::operator*(Matrix & _m) const" << endl;
		exit(9232);
	}
	Matrix out(rows, _m.cols, 0.0);
	for (unsigned int i=0; i<rows; i++) {
		for (unsigned int j=0; j<_m.cols; j++) {
			for (unsigned int k=0; k<_m.rows; k++) {
				out[i][j] += matrix[i][k] * _m.matrix[k][j];
			}
		}
	}

	return out;
}

void Matrix::operator*=(const Matrix & _m) {
	Matrix product = *this * _m;
	*this = product;
}




vector<vector<double> > Matrix::getEigenvectorsEigenValuesGSL(){

	vector<vector<double> > PrincipleComponents(4,vector<double>(5,0.0));

#ifdef __GSL__

  // Make 4x4 out of a 3x3 if needed
  if (getRows() == 3){
      vector<double> valsRow(3,0);
      addRow(valsRow);

      vector<double> valsCol(4,0);
      addCol(valsCol);
      
  }

  if (getRows() != 4){
     cerr << "ERROR 2345 matrix::getEigenvectorsGSL() matrix is not 4x4\n";
     exit(2345);
  }


  //we first need to copy the contents of this matrix into an array...
  double *mtx = new double[16];
    
  uint k = 0;
  for (uint i = 0; i < rows; i++){
      for (uint j = 0; j < cols; j++){
        mtx[k++] = matrix[i][j];
      }

  }
  //now use the GNU Scientific Library to solve the eigenvalue
  //problem for this matrix

  gsl_matrix_view m = gsl_matrix_view_array(mtx,4,4);
    
  //allocate space for the resulting eigenvectors and eigenvalues
  gsl_vector_complex *eig_val = gsl_vector_complex_alloc(4);
  gsl_matrix_complex *eig_vec = gsl_matrix_complex_alloc(4,4); 
  gsl_matrix_complex_set_zero(eig_vec);


  //now allocate some workspace for the calculation...
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc (4);


  //perform the calculation
  gsl_eigen_nonsymmv(&m.matrix, eig_val, eig_vec, w);

    
  //free the space used by the calculation
  gsl_eigen_nonsymmv_free(w);

  //now sort the eigenvectors from the largest eigenvalue to smallest
  gsl_eigen_nonsymmv_sort(eig_val, eig_vec, GSL_EIGEN_SORT_ABS_DESC);


  // now put results in PrincipleComponents
  for (uint i = 0; i < 4; i++){
    
      gsl_complex eval_i = gsl_vector_complex_get (eig_val, i);
      gsl_vector_complex_view eig_vec_i = gsl_matrix_complex_column (eig_vec, i);
     
      for (uint j = 0; j < 4; ++j) {
          gsl_complex z = gsl_vector_complex_get(&eig_vec_i.vector,j);
	  complex<double> foo(GSL_REAL(z),GSL_IMAG(z));

	  // FOR NOW THROWING IMAG part away! std::complex couldn't work in boost::serialization projec
	  //PrincipleComponents[i][j] = foo;
          PrincipleComponents[i][j] = GSL_REAL(z);
      }

      PrincipleComponents[i][4] = GSL_REAL(eval_i);
  }
  

  //free up the memory used by the GSL data...
  gsl_vector_complex_free(eig_val);
  gsl_matrix_complex_free(eig_vec);


  //free up the memory held by the copy of this matrix
  delete[] mtx;

#endif


 return PrincipleComponents;

}

vector<vector<double> > Matrix::getEigenvectorsGSL(){

	vector<vector<double> > PrincipleComponents(4, vector<double>(4, 0.0));

#ifdef __GSL__

  // Make 4x4 out of a 3x3 if needed
  if (getRows() == 3){
      vector<double> valsRow(3,0);
      addRow(valsRow);

      vector<double> valsCol(4,0);
      addCol(valsCol);
      
  }

  if (getRows() != 4){
     cerr << "ERROR 2345 matrix::getEigenvectorsGSL() matrix is not 4x4\n";
     exit(2345);
  }


  //we first need to copy the contents of this matrix into an array...
  double *mtx = new double[16];
    
  uint k = 0;
  for (uint i = 0; i < rows; i++){
      for (uint j = 0; j < cols; j++){
        mtx[k++] = matrix[i][j];
      }

  }
  //now use the GNU Scientific Library to solve the eigenvalue
  //problem for this matrix

  gsl_matrix_view m = gsl_matrix_view_array(mtx,4,4);
    
  //allocate space for the resulting eigenvectors and eigenvalues
  gsl_vector_complex *eig_val = gsl_vector_complex_alloc(4);
  gsl_matrix_complex *eig_vec = gsl_matrix_complex_alloc(4,4); 
  gsl_matrix_complex_set_zero(eig_vec);


  //now allocate some workspace for the calculation...
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc (4);


  //perform the calculation
  gsl_eigen_nonsymmv(&m.matrix, eig_val, eig_vec, w);

    
  //free the space used by the calculation
  gsl_eigen_nonsymmv_free(w);

  //now sort the eigenvectors from the largest eigenvalue to smallest
  gsl_eigen_nonsymmv_sort(eig_val, eig_vec, GSL_EIGEN_SORT_ABS_DESC);


  // now put results in PrincipleComponents
  for (uint i = 0; i < 4; i++){
    
      gsl_vector_complex_view eig_vec_i = gsl_matrix_complex_column (eig_vec, i);
     
      for (uint j = 0; j < 4; ++j) {
          gsl_complex z = gsl_vector_complex_get(&eig_vec_i.vector,j);
	  complex<double> foo(GSL_REAL(z),GSL_IMAG(z));

	  // FOR NOW THROWING IMAG part away! std::complex couldn't work in boost::serialization projec
	  //PrincipleComponents[i][j] = foo;
          PrincipleComponents[i][j] = GSL_REAL(z);
      }


  }


  //free up the memory used by the GSL data...
  gsl_vector_complex_free(eig_val);
  gsl_matrix_complex_free(eig_vec);


  //free up the memory held by the copy of this matrix
  delete[] mtx;

#endif


 return PrincipleComponents;

}

Matrix Matrix::getTranspose() const {
	Matrix m(cols, rows, 0.0);
	
	for (unsigned int i=0; i < rows; i++) {
		for (unsigned int j=0; j < cols; j++) {
			m[j][i] = matrix[i][j];
		}
	}
	return m;
}

Matrix Matrix::getSubMatrix(unsigned int _rowStart, unsigned int _rowEnd, unsigned int _colStart, unsigned int _colEnd) const {
	if (_rowStart > rows || _colStart > cols || _rowEnd > rows || _colEnd > cols) {
		cerr << "ERROR 9189: out of range indeces (" << _rowStart << ", " << _rowEnd << ", " << _colStart << ", " << _colEnd << ") for matrix ("<< rows << ", " << cols << ") in Matrix Matrix::getSubMatrix(unsigned int _rowStart, unsigned int _rowEnd, unsigned int _colStart, unsigned int _colEnd) const" << endl;
		exit(9189);
	}
	if (_rowStart > _rowEnd) {
		unsigned int tmp = _rowStart;
		_rowStart = _rowEnd;
		_rowEnd = tmp;
	}
	if (_colStart > _colEnd) {
		unsigned int tmp = _colStart;
		_colStart = _colEnd;
		_colEnd = tmp;
	}

	Matrix m(_rowEnd - _rowStart + 1, _colEnd - _colStart + 1, 0.0);
	
	for (unsigned int i=_rowStart, ii=0; i<=_rowEnd; i++, ii++) {
		for (unsigned int j=_colStart, jj=0; j<=_colEnd; j++, jj++) {
			m[ii][jj] = matrix[i][j];
		}
	}
	return m;
}

double Matrix::getDeterminant() const {

	// Check for non-square matrix
	if (rows != cols) {
		cerr << "ERROR 9203: cannot get determinant for non-square matrix ("<< rows << ", " << cols << ") in double Matrix::getDeterminant() const" << endl;
		exit(9203);
	}
	
	// recursive algorithm until size 2x2 is reached
	double det = 0.0;
	if (rows == 2) {
		det =  matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	} else {
		for (unsigned int i=0; i<cols; i++) {
			det += pow(-1.0, (double)i) * matrix[0][i] * getMinor(0, i).getDeterminant();
		}
	}
	return det;

}

Matrix Matrix::getMinor(unsigned int _row, unsigned int _col) const {
	if (matrix.size() == 0) {
		cerr << "ERROR 9187: cannot get minor for zero dimension matrix ("<< rows << ", " << cols << ") in Matrix Matrix::getMinor(unsigned int _row, unsigned int _col)" << endl;
		exit(9187);
	}
	if (_row >= rows || _col >= cols) {
		cerr << "ERROR 9187: out of range indeces ("<< _row << ", " << _col << ") for matrix ("<< rows << ", " << cols << ") in Matrix Matrix::getMinor(unsigned int _row, unsigned int _col)" << endl;
		exit(9187);
	}

	Matrix m(rows - 1 , cols - 1, 0.0);
	
	for (unsigned int i=0; i<rows-1; i++) {
		int incrementedI = i;
		if (i>= _row) {
			incrementedI = i+1;
		}
		for (unsigned int j=0; j<cols-1; j++) {
			int incrementedJ = j;
			if (j>= _col) {
				incrementedJ = j+1;
			}
			m[i][j] = matrix[incrementedI][incrementedJ];
		}
	}
	return m;
}

void Matrix::addRow(vector<double> _val) {
	// case emtpy matrix
	if (rows == 0 && cols == 0) {
		rows = 1;
		cols = _val.size();
	} else {

		if (_val.size() != cols) {
			cerr << "ERROR 9181: vector size (" << _val.size() << ") out of range in void Matrix::addRow(vector<double> _val)" << endl;
			exit(9181);
		}
		rows++;
	}
	matrix.push_back(_val);
}

void Matrix::addCol(vector<double> _val) {
	// case emtpy matrix
	if (rows == 0 && cols == 0) {
		rows = _val.size();
		cols = 1;
		matrix = vector<vector<double> >(_val.size(), vector<double>());
	} else {
		if (_val.size() != rows) {
			cerr << "ERROR 9183: vector size (" << _val.size() << ") out of range in void Matrix::addCol(vector<double> _val)" << endl;
			exit(9183);
		}
		cols++;
	}
	for (unsigned int i=0; i<rows; i++) {
		matrix[i].push_back(_val[i]);
	}
}

string Matrix::toString() const {

	stringstream s;
	s << "["<<endl;
	for (uint i = 0; i < matrix.size(); i++){
		for (uint j = 0; j < matrix[i].size();j++){
			char tmp[12];
			sprintf(tmp, "%8.3f ",matrix[i][j]);
			s << tmp;
		}
		s <<endl;
	}
	s << "]"<<endl;

	return s.str();
}

