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

#ifndef __BOOST__
#error message("Boost libraries not defined, can't use testBoost without them.")
#endif

#include <fstream>
#include <string>
#include <iostream>
#include <vector>

//#include "BoostGPS.cpp"
#include "CartesianPoint.h"
#include "Atom.h"
#include "AtomPointerVector.h"

using namespace std;
using namespace MSL;




#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


int main() {


    CartesianPoint p1(1.1,2.2,3.3);
    p1.save_checkpoint("point.ckpt");

    
    CartesianPoint p2;
    p2.load_checkpoint("point.ckpt");

    cout << "CartesianPoints\n";
    cout << "P1: "<<p1.toString()<<endl;
    cout << "P2: "<<p2.toString()<<endl;


    Matrix m1(3,3,5.5);
    m1.save_checkpoint("matrix.ckpt");
    
    Matrix m2;
    m2.load_checkpoint("matrix.ckpt");


    cout << "Matrix\n";
    cout << "M1: "<<endl;
    for (uint i = 0 ; i < 3; i++){
	    for(uint j = 0; j < 3; j++){
		    fprintf(stdout, "%8.3f ",m1[i][j]);
	    }
	    fprintf(stdout,"\n");
    }
    cout << "M2: "<<endl;
    for (uint i = 0 ; i < 3; i++){
	    for(uint j = 0; j < 3; j++){
		    fprintf(stdout, "%8.3f ",m1[i][j]);
	    }
	    fprintf(stdout,"\n");
    }
    


    Atom a;
    a.setResidueName("ALA");
    a.setCoor(1.0,1.0,1.0);
    a.setName("CA");

    a.save_checkpoint("atom.ckpt");


    Atom b;
    b.load_checkpoint("atom.ckpt");

    fprintf(stdout, "A: %s\n",a.toString().c_str());
    fprintf(stdout, "B: %s\n",b.toString().c_str());


    cout << "AtomPointerVector1: "<<endl;
    AtomPointerVector av;
    av.push_back(&a);
    for (uint i = 0;i < av.size();i++){
	    cout << av(i)<<endl;
    }
    av.save_checkpoint("atomvector.ckpt");
    
    AtomPointerVector av2;
    av2.load_checkpoint("atomvector.ckpt");

    cout << "AtomPointerVector2: "<<endl;
    for (uint i = 0;i < av2.size();i++){
	    cout << av2(i)<<endl;
    }


    // testing
    AtomPointerVector av3;
    av3.load_checkpoint("/Users/dwkulp/software/mslib/out.db");
    cout << "AV3.size(): "<<av3.size()<<endl;

    return 0;
}
