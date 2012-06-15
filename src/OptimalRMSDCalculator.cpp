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

#include <iostream>
#include <math.h>
#include "OptimalRMSDCalculator.h"
#include "AtomPointerVector.h"

using namespace std;
using namespace MSL;

vector<double> OptimalRMSDCalculator::lastTranslation() {
	vector<double> trans(3, 0.0);
	for (int i = 0; i < 3; i++) trans[i] = t[i];
	return trans;
}

vector<vector<double> > OptimalRMSDCalculator::lastRotation() {
	vector<vector<double> > rot(3);
	for (int i = 0; i < 3; i++) {
		rot[i].resize(3, 0);
		for (int j = 0; j < 3; j++ ) {
			rot[i][j] = u[i][j];
		}
	}
	return rot;
}

double OptimalRMSDCalculator::bestRMSD(AtomPointerVector &_align, AtomPointerVector &_ref, bool* _suc) {
    rmsd = 999999.0;
    if (Kabsch(_align, _ref, 0)) { if (_suc != NULL) *_suc = true; }
    else { if (_suc != NULL) *_suc = false; }
    return rmsd;
}

bool OptimalRMSDCalculator::align(AtomPointerVector &_align, AtomPointerVector &_ref, AtomPointerVector& _moveable) {
    rmsd = 999999.0;
	bool suc = Kabsch(_align, _ref, 1);

	if (suc) {
		double x[3],x1[3];
		for(int k=0; k<_moveable.size(); k++) {
			x[0]=_moveable[k]->getX();
			x[1]=_moveable[k]->getY();
			x[2]=_moveable[k]->getZ();
			x1[0] = t[0]+u[0][0]*x[0]+u[0][1]*x[1]+u[0][2]*x[2];
			x1[1] = t[1]+u[1][0]*x[0]+u[1][1]*x[1]+u[1][2]*x[2];
			x1[2] = t[2]+u[2][0]*x[0]+u[2][1]*x[1]+u[2][2]*x[2];
			_moveable[k]->setCoor(x1[0],x1[1],x1[2]);
		}
	}
	return suc;
}


/**************************************************************************
  Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
  x    - x(i,m) are coordinates of atom m in set x            (input)
  y    - y(i,m) are coordinates of atom m in set y            (input)
  n    - n is number of atom pairs                            (input)
  mode  - 0:calculate rmsd only                               (input)
          1:calculate rmsd,u,t                                (takes longer)
  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
  u    - u(i,j) is   rotation  matrix for best superposition  (output)
  t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
bool OptimalRMSDCalculator::Kabsch(AtomPointerVector &_align, AtomPointerVector &_ref, int mode) {
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;  
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3=1.73205080756888, tol=0.01;
	int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
	int ip2312[]={1, 2, 0, 1};
	
	int a_failed=0, b_failed=0;
	double epsilon=0.00000001;
	
	int n=_ref.size();
	if(n != _align.size()) {
		cout << "Two proteins have different length!" << endl;
		return false;
	}
	
	//initializtation
	rmsd=0;
	rms1=0;
	e0=0;
	for (i=0; i<3; i++) {
		xc[i]=0.0;
		yc[i]=0.0;
		t[i]=0.0;
		for (j=0; j<3; j++) {
			u[i][j]=0.0;
			r[i][j]=0.0;
			a[i][j]=0.0;
			if (i==j) {
				u[i][j]=1.0;
				a[i][j]=1.0;
			}
		}
	}
	
	if (n<1) {
		cout << "Protein length is zero!" << endl;
		return false;
	} 

	//compute centers for vector sets x, y
	for(i=0; i<n; i++){
		xc[0] += _align[i]->getX();
		xc[1] += _align[i]->getY();
		xc[2] += _align[i]->getZ();
		
		yc[0] += _ref[i]->getX();
		yc[1] += _ref[i]->getY();
		yc[2] += _ref[i]->getZ();
	}
	for(i=0; i<3; i++){
		xc[i] = xc[i]/(double)n;
		yc[i] = yc[i]/(double)n;        
	}
	
	//compute e0 and matrix r
	for (m=0; m<n; m++) {
		e0 += (_align[m]->getX()-xc[0])*(_align[m]->getX()-xc[0]) \
		  +(_ref[m]->getX()-yc[0])*(_ref[m]->getX()-yc[0]);
		e0 += (_align[m]->getY()-xc[1])*(_align[m]->getY()-xc[1]) \
		  +(_ref[m]->getY()-yc[1])*(_ref[m]->getY()-yc[1]);
		e0 += (_align[m]->getZ()-xc[2])*(_align[m]->getZ()-xc[2]) \
		  +(_ref[m]->getZ()-yc[2])*(_ref[m]->getZ()-yc[2]);
		r[0][0] += (_ref[m]->getX() - yc[0])*(_align[m]->getX() - xc[0]);
		r[0][1] += (_ref[m]->getX() - yc[0])*(_align[m]->getY() - xc[1]);
		r[0][2] += (_ref[m]->getX() - yc[0])*(_align[m]->getZ() - xc[2]);
		r[1][0] += (_ref[m]->getY() - yc[1])*(_align[m]->getX() - xc[0]);
		r[1][1] += (_ref[m]->getY() - yc[1])*(_align[m]->getY() - xc[1]);
		r[1][2] += (_ref[m]->getY() - yc[1])*(_align[m]->getZ() - xc[2]);
		r[2][0] += (_ref[m]->getZ() - yc[2])*(_align[m]->getX() - xc[0]);
		r[2][1] += (_ref[m]->getZ() - yc[2])*(_align[m]->getY() - xc[1]);
		r[2][2] += (_ref[m]->getZ() - yc[2])*(_align[m]->getZ() - xc[2]);
	}
	//compute determinat of matrix r
	det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )		\
	- r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )		\
	+ r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] ); 
	sigma=det;
	
	//compute tras(r)*r
	m = 0;
	for (j=0; j<3; j++) {
		for (i=0; i<=j; i++) {            
			rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j];
			m++;
		}
	}
	
	double spur=(rr[0]+rr[2]+rr[5]) / 3.0;
	double cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5])	\
		  - rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
	det = det*det; 
	
	for (i=0; i<3; i++){
		e[i]=spur;
	}
	
	if (spur>0) {
		d = spur*spur;
		h = d - cof;
		g = (spur*cof - det)/2.0 - spur*h;

		if (h>0) {
			sqrth = sqrt(h);
			d = h*h*h - g*g;
			if(d<0.0) d=0.0;
			d = atan2( sqrt(d), -g ) / 3.0;			
			cth = sqrth * cos(d);
			sth = sqrth*sqrt3*sin(d);
			e[0]= (spur + cth) + cth;
			e[1]= (spur - cth) + sth;            
			e[2]= (spur - cth) - sth;

			if (mode!=0) {//compute a                
				for (l=0; l<3; l=l+2) {
					d = e[l];  
					ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
					ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
					ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
					ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
					ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];                
					ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1]; 

					if (fabs(ss[0])<=epsilon) ss[0]=0.0;
					if (fabs(ss[1])<=epsilon) ss[1]=0.0;
					if (fabs(ss[2])<=epsilon) ss[2]=0.0;
					if (fabs(ss[3])<=epsilon) ss[3]=0.0;
					if (fabs(ss[4])<=epsilon) ss[4]=0.0;
					if (fabs(ss[5])<=epsilon) ss[5]=0.0;

					if (fabs(ss[0]) >= fabs(ss[2])) {
						j=0;                    
						if( fabs(ss[0]) < fabs(ss[5])){
							j = 2;
						}
					} else if ( fabs(ss[2]) >= fabs(ss[5]) ){
						j = 1;
					} else {
						j = 2;
					}

					d = 0.0;
					j = 3 * j;
					for (i=0; i<3; i++) {
						k=ip[i+j];
						a[i][l] = ss[k];
						d = d + ss[k]*ss[k];						
					} 


					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > epsilon) d = 1.0 / sqrt(d);
					else d=0.0;
					for (i=0; i<3; i++) {
						a[i][l] = a[i][l] * d;
					}               
				}//for l

				d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2])) {
					m1=2;
					m=0;
				} else {
					m1=0;
					m=2;                
				}
				p=0;
				for(i=0; i<3; i++){
					a[i][m1] = a[i][m1] - d*a[i][m];
					p = p + a[i][m1]*a[i][m1];
				}
				if (p <= tol) {
					p = 1.0;
					for (i=0; i<3; i++) {
						if (p < fabs(a[i][m])){
							continue;
						}
						p = fabs( a[i][m] );
						j = i;                    
					}
					k = ip2312[j];
					l = ip2312[j+1];
					p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] ); 
					if (p > tol) {
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m]/p;
						a[l][m1] =  a[k][m]/p;                                                       
					} else {//goto 40
						a_failed=1;
					}     
				} else {//if p<=tol
					p = 1.0 / sqrt(p);
					for(i=0; i<3; i++){
						a[i][m1] = a[i][m1]*p;
					}                                  
				}//else p<=tol  
				if (a_failed!=1) {
					a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
					a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
					a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];       
				}                                   
			}//if(mode!=0)       
		}//h>0

		//compute b anyway
		if (mode!=0 && a_failed!=1) {//a is computed correctly
			//compute b
			for (l=0; l<2; l++) {
				d=0.0;
				for(i=0; i<3; i++){
					b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
					d = d + b[i][l]*b[i][l];
				}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > epsilon) d = 1.0 / sqrt(d);
				else d=0.0;
				for (i=0; i<3; i++) {
					b[i][l] = b[i][l]*d;
				}                
			}            
			d = b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
			p=0.0;

			for (i=0; i<3; i++) {
				b[i][1] = b[i][1] - d*b[i][0];
				p += b[i][1]*b[i][1];
			}

			if (p <= tol) {
				p = 1.0;
				for (i=0; i<3; i++) {
					if (p<fabs(b[i][0])) {
						continue;
					}
					p = fabs( b[i][0] );
					j=i;
				}
				k = ip2312[j];
				l = ip2312[j+1];
				p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] ); 
				if (p > tol) {
					b[j][1] = 0.0;
					b[k][1] = -b[l][0]/p;
					b[l][1] =  b[k][0]/p;        
				} else {
					//goto 40
					b_failed=1;
				}                
			} else {//if( p <= tol )
				p = 1.0 / sqrt(p);
				for(i=0; i<3; i++){
					b[i][1]=b[i][1]*p;
				}
			}            
			if (b_failed!=1){
				b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
				b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
				b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0]; 
				//compute u
				for (i=0; i<3; i++){
					for(j=0; j<3; j++){
						u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1]	\
								+ b[i][2]*a[j][2];
					}
				}
			}

			//compute t
			for(i=0; i<3; i++){
				t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1])	\
						- u[i][2]*xc[2];
			}            
		}//if(mode!=0 && a_failed!=1)
	} else {//spur>0, just compute t and errors
		//compute t
		for (i=0; i<3; i++) {
			t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
		}
	} //else spur>0 

	//compute rmsd
	for(i=0; i<3; i++){
		if( e[i] < 0 ) e[i] = 0;
		e[i] = sqrt( e[i] );           
	}            
	d = e[2];
	if( sigma < 0.0 ){
		d = - d;
	}
	d = (d + e[1]) + e[0];
	rms1 = (e0 - d) - d; 
	if( rms1 < 0.0 ) rms1 = 0.0;  

	rmsd=sqrt(rms1/(double)n);

	return true;
}
