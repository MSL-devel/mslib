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


#include <RInside.h>                    // for the embedded R via RInside

using namespace std;

int main(int argc, char *argv[]) {

    // create an embedded R instance 
    RInside R;

    /* 
       rinside_sample0 from RInside package
	 Shows
	     how you can call R commands using a parseEvalQ(string);
     */
    cout << "\n\n ************   TEST 1 : Hello World ***************** \n\n\n";
    R["txt"] = "Hello, world!\n";	// assign a char* (string) to 'txt'
    R.parseEvalQ("cat(txt)");           // eval the init string, ignoring any returns




    /* 
       rinside_sample1 from RInside package
        Shows 
              how you can get data into R and 
              how you can get data back from R
    */
    cout << "\n\n ************   TEST 2 : Passing Data ***************** \n\n\n";
    SEXP ans;                                 // SEXP data type for getting data back from R

    // Create a NumericMatrix  (4 by 4)
    Rcpp::NumericMatrix M(4,4);
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            M(i,j) = i*10+j; 
        }
    }

    R["M"] = M;                                 // assign C++ matrix M to R's 'M' var

    std::string evalstr = "\
        cat('Running ls()\n'); print(ls());                    \
        cat('Showing M\n'); print(M);                          \
        cat('Showing colSums()\n'); Z <- colSums(M); print(Z); \
        Z";                     // returns Z

    ans = R.parseEval(evalstr);                 // eval the init string -- Z is now in ans
                                                
    Rcpp::NumericVector v(ans);                 // convert SEXP ans to a vector of doubles
    for (int i=0; i< v.size(); i++) {           // show the result
      std::cout << "In C++ element " << i << " is " << v[i] << " corresponds to sum of column "<<i+1<<" of matrix M"<<std::endl;
    }


    cout << "\n\n ************* Done *************"<<endl;

    exit(0);
}

