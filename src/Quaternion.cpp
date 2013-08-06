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

#include "Quaternion.h"

using namespace MSL;
using namespace std;




Quaternion::Quaternion() {
	a = 0.0;
	v = vector<double>(3, 0.0);
}

Quaternion::Quaternion(double _a, double _v0, double _v1, double _v2) {
	a = _a;
	v = vector<double>(3, 0.0);
	v[0] = _v0;
	v[1] = _v1;
	v[2] = _v2;
}

Quaternion::Quaternion(double _a, vector<double> _v) {
	a = _a;
	if (_v.size() != 3) {
		cerr << "ERROR 1740: vector argument must have size 3 in Quaternion::Quaternion(double _a, vector<double> _v)" << endl;
		exit(1740);
	}
	v = _v;

}

Quaternion::~Quaternion() {
}

void Quaternion::operator+=(const Quaternion & _q) {
	a += _q.a;
	v[0] += _q.v[0];
	v[1] += _q.v[1];
	v[2] += _q.v[2];
}

void Quaternion::operator-=(const Quaternion & _q) {
	a -= _q.a;
	v[0] -= _q.v[0];
	v[1] -= _q.v[1];
	v[2] -= _q.v[2];
}

void Quaternion::operator*=(double _b) {
	a *= _b;
	v[0] *= _b;
	v[1] *= _b;
	v[2] *= _b;

}

void Quaternion::operator/=(double _b) {
	a /= _b;
	v[0] /= _b;
	v[1] /= _b;
	v[2] /= _b;
}

void Quaternion::operator*=(const Quaternion & _q) {
	/**********************************************
	 * if quaternions are (a, A) and (b, B):
	 * 
	 *   ab = (a0*b0 - A.B, bA + aB + AxB)
	 **********************************************/
	a = a * _q.a - dotVec(v, _q.v);
	
	vector<double> aB = multiplyVec(_q.v, a);
	vector<double> bA = multiplyVec(v, _q.a);
	vector<double> AxB = crossVec(v, _q.v);
	v[0] = aB[0] + bA[0] + AxB[0];
	v[1] = aB[1] + bA[1] + AxB[1];
	v[2] = aB[2] + bA[2] + AxB[2];
}

Quaternion Quaternion::operator+(const Quaternion & _q) {
	Quaternion out = *this;
	out += _q;
	return out;
}

Quaternion Quaternion::operator-(const Quaternion & _q) {
	Quaternion out = *this;
	out -= _q;
	return out;
}
	
Quaternion Quaternion::operator*(const Quaternion & _q) {
	Quaternion out = *this;
	out *= _q;
	return out;
}

Quaternion Quaternion::operator*(double _b) {
	Quaternion out = *this;
	out *= _b;
	return out;
}

Quaternion Quaternion::operator/(double _b) {
	Quaternion out = *this;
	out /= _b;
	return out;
}

void Quaternion::setConiugate() {
	v[0] *= -1;
	v[1] *= -1;
	v[2] *= -1;
}

Quaternion Quaternion::getConiugate() const {
	Quaternion out = *this;
	out.setConiugate();
	return out;
}

double Quaternion::norm() const {
	return a*a + dotVec(v, v);
}

void Quaternion::setInverse() {
	Quaternion q = getConiugate();
	q /= norm();
	a = q.a;
	v = q.v;
}

Quaternion Quaternion::getInverse() const {
	Quaternion out = *this;
	out.setInverse();
	return out;
}
	
vector<double>Quaternion::multiplyVec(vector<double> _v, double _b) const {
	vector<double> out = _v;

	out[0] *= _b;
	out[1] *= _b;
	out[2] *= _b;

	return out;
}

double Quaternion::dotVec(vector<double> _v1, vector<double> _v2) const {
	return _v1[0] * _v2[0] + _v1[1] * _v2[1] + _v1[2] * _v2[2];
}

vector<double>Quaternion::crossVec(vector<double> _v1, vector<double> _v2) const {
	vector<double> out(3, 0.0);

	out[0] = _v1[1] * _v2[2] - _v1[2] * _v2[1];
	out[1] = _v1[2] * _v2[0] - _v1[0] * _v2[2];
	out[2] = _v1[0] * _v2[1] - _v1[1] * _v2[0];

	return out;
}



bool Quaternion::convertToRotationMatrix(Matrix &matrix){

  double n = norm();

  double s = (n > 0.0) ? 2.0 / n : 0.0;

  double xx,yy,zz,ww,xs, ys, zs, wxs, wys, wzs, xys, xzs, yzs;

  xx  = v[0]*v[0]; yy = v[1]*v[1]; zz = v[2]*v[2]; ww = a*a;
  xs  = v[0]*s;    ys = v[1]*s;    zs = v[2]*s;
  xys = v[0]*ys;  xzs = v[0]*zs; 
  yzs = v[1]*zs; 
  wxs = a*xs;     wys = a*ys;     wzs = a*zs; 

  // Resize to insure proper dimensions 
  matrix[0][0] = xx + yy -zz - ww;
  matrix[1][0] = yzs + wxs;
  matrix[2][0] = wys - xzs ;

  matrix[0][1] = yzs - wxs;
  matrix[1][1] = xx + zz - yy - ww;
  matrix[2][1] = xys + wzs;

  matrix[0][2] = xzs + wys;
  matrix[1][2] = wzs - xys;
  matrix[2][2] = xx + ww - yy - zz;

  return true;
}

bool Quaternion::convertToQuaternion(Matrix &matrix){

  int i,j,k;
  i = 0; j = 1; k = 2;  

  double tr = matrix[i][i] + matrix[j][j] + matrix[k][k];
  double r,s;
  r = s = 0.0;
  
  if (tr < 0.0){

    if (matrix[j][j] > matrix[i][i]) { i = 1; j = 2; k = 0; }
    if (matrix[2][2] > matrix[i][i]) { i = 2; j = 0; k = 1; }

    r = sqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + 1.0);
    s = 0.5 / r;

    v[0] = 0.5 * r;
    v[1] = (matrix[i][j] + matrix[j][i])*s;
    v[2] = (matrix[k][i] + matrix[i][k])*s;
    a    = (matrix[k][j] - matrix[j][k])*s;

  }else {

    r = sqrt(tr + 1.0);
    s = 0.5 / r;

    v[0] = (matrix[k][j] - matrix[j][k])*s;
    v[1] = (matrix[i][k] - matrix[k][i])*s;
    v[2] = (matrix[j][i] - matrix[i][j])*s;
    a    = 0.5 * r;
    
  }
  

  return true;
}

/* Make a quaternion for rotation purposes from an axis of rotation and an angle */
bool Quaternion::makeQuaternion(CartesianPoint &axis, const double theta){

  v[0] = axis[0] * sin(theta/2);
  v[1] = axis[1] * sin(theta/2);
  v[2] = axis[2] * sin(theta/2);
  a    = cos(theta/2);
  
  return true;
}

#ifdef __GSL__
bool Quaternion::makeQuaternion(AtomPointerVector &_align, AtomPointerVector &_ref){


  if (_ref.size() != _align.size()) return false;

  /*
  GC1[0] = GC1[1] = GC1[2] = GC2[0] = GC2[1] = GC2[2] = 0.0; 
  for(int i=0;i<_ref.size();i++){
    GC1[0] += (*_ref[i])[0]  ; GC1[1] += (*_ref[i])[1]   ;GC1[2] += (*_ref[i])[2];
    GC2[0] += (*_align[i])[0]; GC2[1] += (*_align[i])[1] ;GC2[2] += (*_align[i])[2];
  }

  for(int i=0;i<3;i++){
    GC1[i]=GC1[i]/(double)_ref.size();
    GC2[i]=GC2[i]/(double)_align.size();
  }
  */

  //_ref.updateGeometricCenter();
  //_align.updateGeometricCenter();
  CartesianPoint GC1 = _ref.getGeometricCenter();
  CartesianPoint GC2 = _align.getGeometricCenter();

 // Calculate 4x4 quaternion matrix..
 vector<vector<double > > qmat;
 qmat.resize(4);
 qmat[0].resize(4); 
 qmat[1].resize(4); 
 qmat[2].resize(4); 
 qmat[3].resize(4); 


 for (uint i = 0; i < _ref.size(); i++){
 
   Quaternion m;
   m.setX(((*_ref[i])[0]-GC1[0]) - ((*_align[i])[0]-GC2[0]));
   m.setY(((*_ref[i])[1]-GC1[1]) - ((*_align[i])[1]-GC2[1]));
   m.setZ(((*_ref[i])[2]-GC1[2]) - ((*_align[i])[2]-GC2[2]));
   m.setW(0);

   Quaternion p;
   p.setX(((*_ref[i])[0]-GC1[0]) + ((*_align[i])[0]-GC2[0]));
   p.setY(((*_ref[i])[1]-GC1[1]) + ((*_align[i])[1]-GC2[1]));
   p.setZ(((*_ref[i])[2]-GC1[2]) + ((*_align[i])[2]-GC2[2]));
   p.setW(0);

 
   qmat[0][0] += m.norm(); 
   qmat[0][1] += p.getY() * m.getZ() - m.getY()*p.getZ();
   qmat[0][2] += m.getX() * p.getZ() - p.getX()*m.getZ();
   qmat[0][3] += p.getX() * m.getY() - m.getX()*p.getY();
 
   qmat[1][0] +=  p.getY()*m.getZ()-m.getY()*p.getZ();
   qmat[1][1] +=  p.getY()*p.getY()+p.getZ()*p.getZ()+m.getX()*m.getX();
   qmat[1][2] +=  m.getX()*m.getY()-p.getX()*p.getY();
   qmat[1][3] +=  m.getX()*m.getZ()-p.getX()*p.getZ();                        
 
   qmat[2][0] +=  m.getX()*p.getZ()-p.getX()*m.getZ();
   qmat[2][1] +=  m.getX()*m.getY()-p.getX()*p.getY();
   qmat[2][2] +=  p.getX()*p.getX()+p.getZ()*p.getZ()+m.getY()*m.getY();
   qmat[2][3] +=  m.getY()*m.getZ()-p.getY()*p.getZ();
 
   qmat[3][0] +=  p.getX()*m.getY()-m.getX()*p.getY();     
   qmat[3][1] +=  m.getX()*m.getZ()-p.getX()*p.getZ();
   qmat[3][2] +=  m.getY()*m.getZ()-p.getY()*p.getZ();
   qmat[3][3] +=  p.getX()*p.getX()+p.getY()*p.getY()+m.getZ()*m.getZ();
 
  }

  //diagonalise this matrix
  getPrincipalAxes(qmat);


  // The rotation matrix that maximises the RMSD is the smallest (4th) eigenvector 
  // So store this quaternion
  setX(qmat[0][3]);
  setY(qmat[1][3]);
  setZ(qmat[2][3]);
  setW(qmat[3][3]);
 
  return true;
}
#endif

bool Quaternion::rotatePoint(CartesianPoint &pt){
  return rotatePoint(pt,pt);
}

bool Quaternion::rotatePoint(CartesianPoint &pt,CartesianPoint &pt_rot){

  Quaternion v;
  v.setX(pt[0]);
  v.setY(pt[1]);
  v.setZ(pt[2]);
  v.setW(0);

  Quaternion conjugate = getConiugate();
  Quaternion result = (*this) * v * conjugate;

  // Is this '= result.getX()' or '= pt[0] + result.getX()' ?
  pt_rot[0] = result.getX();
  pt_rot[1] = result.getY();
  pt_rot[2] = result.getZ();
  
  return true;
}


string Quaternion::toString() {

  // Output Quaternion values
  char Q[100];
  sprintf(Q,"Quaternion:\n\t% 10.6f  % 10.6f  % 10.6f  % 10.6f \n",v[0],v[1],v[2],a);
  string Qstr(Q);

  // Create the corresponding rotation matrix
  Matrix m(3,3);
  convertToRotationMatrix(m);

  char RM[1000];
  sprintf(RM,"Rotation Matrix:\n\t% 10.6f  % 10.6f  % 10.6f\n\t% 10.6f  % 10.6f  % 10.6f\n\t% 10.6f  % 10.6f  % 10.6f\n",m[0][0],m[0][1],m[0][2],
	  m[1][0],m[1][1],m[1][2],
	  m[2][0],m[2][1],m[2][2]);

  string RMstr(RM);


  string result = Qstr + RMstr;

  

  return result;
}

#ifdef __GSL__
void Quaternion::getPrincipalAxes(vector<vector<double> > &mat) {



  //assume that this matrix is symmetrical - we will thus
  //only look at the values in the upper-right diagonal
    
  //we first need to copy the contents of this matrix into an array...
  double *sym_mtx = new double[16];
    
  int x,y,z,a;
  x = 0; y = 1; z = 2; a = 3;

  sym_mtx[0]  = mat[x][x]; 
  sym_mtx[1]  = mat[x][y];
  sym_mtx[2]  = mat[x][z];
  sym_mtx[3]  = mat[x][a];
    
  sym_mtx[4]  = mat[x][y];
  sym_mtx[5]  = mat[y][y];
  sym_mtx[6]  = mat[y][z];
  sym_mtx[7]  = mat[y][a];
    
  sym_mtx[8]  = mat[x][z];
  sym_mtx[9]  = mat[y][z];
  sym_mtx[10] = mat[z][z];
  sym_mtx[11] = mat[z][a];
    
  sym_mtx[12] = mat[x][a];
  sym_mtx[13] = mat[y][a];
  sym_mtx[14] = mat[z][a];
  sym_mtx[15] = mat[a][a];    
    



  //now use the GNU Scientific Library to solve the eigenvalue
  //problem for this matrix

  gsl_matrix_view m = gsl_matrix_view_array(sym_mtx,4,4);
    
  //allocate space for the resulting eigenvectors and eigenvalues
  gsl_vector *eig_val = gsl_vector_alloc(4);
  gsl_matrix *eig_vec = gsl_matrix_alloc(4,4); gsl_matrix_set_zero(eig_vec);
    
  //now allocate some workspace for the calculation...
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4);
  
  //perform the calculation
  gsl_eigen_symmv(&m.matrix, eig_val, eig_vec, w);
    
  //free the space used by the calculation
  gsl_eigen_symmv_free(w);

  //now sort the eigenvectors from the smallest eigenvalue to 
  //the largest
  gsl_eigen_symmv_sort (eig_val, eig_vec, GSL_EIGEN_SORT_ABS_ASC);
    
  //now copy the results back into a new Matrix
  mat[x][x] = gsl_matrix_get(eig_vec,x,a);
  mat[x][y] = gsl_matrix_get(eig_vec,x,z);
  mat[x][z] = gsl_matrix_get(eig_vec,x,y);
  mat[x][a] = gsl_matrix_get(eig_vec,x,x);

  mat[y][x] = gsl_matrix_get(eig_vec,y,a);
  mat[y][y] = gsl_matrix_get(eig_vec,y,z);
  mat[y][z] = gsl_matrix_get(eig_vec,y,y);
  mat[y][a] = gsl_matrix_get(eig_vec,y,x);

  mat[z][x] = gsl_matrix_get(eig_vec,z,a);
  mat[z][y] = gsl_matrix_get(eig_vec,z,z);
  mat[z][z] = gsl_matrix_get(eig_vec,z,y);
  mat[z][a] = gsl_matrix_get(eig_vec,z,x);

  mat[a][x] = gsl_matrix_get(eig_vec,a,a);
  mat[a][y] = gsl_matrix_get(eig_vec,a,z);
  mat[a][z] = gsl_matrix_get(eig_vec,a,y);
  mat[a][a] = gsl_matrix_get(eig_vec,a,x);

  //free up the memory used by the GSL data...
  gsl_vector_free(eig_val);
  gsl_matrix_free(eig_vec);



  //free up the memory held by the copy of this matrix
  delete[] sym_mtx;

//    Matrix m(mat);
//    vector<vector<complex<double > > > cm = m.getEigenvectors();
//
//    for (uint i = 0; i < 4; i++){
//      for (uint j = 0; j < 4; j++){
//	cout << cm[i][j] << endl;
//	mat[i][j] = cm[i][j].real();
//      }
//    }


}    
#endif
