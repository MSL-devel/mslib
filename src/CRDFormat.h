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

#ifndef CRDFORMAT_H
#define CRDFORMAT_H
/*
  CRD Formating information goes here.  As well as a simple structure for containing a single line
  of a CRD.

  For each line type define start, stop, length, justification of fixed-width fields.
  The format should be:
  S_NAME_OF_FIELD  = NUM,  for starting positions
  E_NAME_OF_FIELD  = NUM,  for ending positions
  L_NAME_OF_FIELD  = NUM,  for length of field
  J_NAME_OF_FIELD  = NUM,  for justification ; 0 = right, 1 = left jusftified
  
  E.G.

  enum NAME {
   S_ATOM_NAME = 12,  
   E_ATOM_NAME = 15,  
   L_ATOM_NAME = 4,  
   J_ATOM_NAME = 1,  
  };


Example CRD lines (Put standard + any strange ones here)

          1         2         3         4         5         6         7                     
01234567890123456789012345678901234567890123456789012345678901234567890123456789
<---><--->x<-->x<--><--------><--------><-------->x<-->x<--><-------->
    1    1 ALA  N      2.14273   1.32845   0.00000 A    1     -0.30000
    2    1 ALA  HT1    3.17887   1.23889   0.00000 A    1      0.33000
    3    1 ALA  HT2    1.84011   1.84730   0.84901 A    1      0.33000


format from io.doc:

	title
	NATOM (I5)
	ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
	I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5
		 1         2         3         4         5         6         7
	----+----|----+----|----+----|----+----|----+----|----+----|----+----|
	<---><--->x<-->x<--><--------><--------><-------->x<-->x<--><-------->
	    1    1 ALA  N      2.14273   1.32845   0.00000 A    1     -0.30000
	    2    1 ALA  HT1    3.17887   1.23889   0.00000 A    1      0.33000
	    3    1 ALA  HT2    1.84011   1.84730   0.84901 A    1      0.33000
	*/





//MSL Includes
#include "MslTools.h"
#include "Atom.h"

//STL Includes
#include <string>
#include <iostream>
#include <cstring>


namespace MSL { 
class CRDFormat {
	public:

		enum AtomField { 

			// Start Record Position
			S_ATOM_NO    =  0,
			S_ABS_RES    =  5,
			S_RES_NAME   = 11,
			S_ATOM_NAME  = 16,
			S_X          = 20,
			S_Y          = 30,
			S_Z          = 40,
			S_CHAIN_ID   = 51,
			S_RES_NUM    = 56,
			S_CHARGE     = 60,

			// End Record Position
			E_ATOM_NO    =  4,
			E_ABS_RES    =  9,
			E_RES_NAME   = 14,
			E_ATOM_NAME  = 19,
			E_X          = 29,
			E_Y          = 39,
			E_Z          = 49,
			E_CHAIN_ID   = 54,
			E_RES_NUM    = 59,
			E_CHARGE     = 69,


			// Length of Record
			L_ATOM_NO    =  5,
			L_ABS_RES    =  5,
			L_RES_NAME   =  4,
			L_ATOM_NAME  =  4,
			L_X          = 10,
			L_Y          = 10,
			L_Z          = 10,
			L_CHAIN_ID   =  4,
			L_RES_NUM    =  4,
			L_CHARGE     = 10,
			

			// Justification within record 0 = right-justified, 1 = left-justified.
			J_ATOM_NO    =  0,
			J_ABS_RES    =  0,
			J_RES_NAME   =  1,
			J_ATOM_NAME  =  1,
			J_X          =  0,
			J_Y          =  0,
			J_Z          =  0,
			J_CHAIN_ID   =  1,
			J_RES_NUM    =  1,
			J_CHARGE     =  0,

			
		};

		struct AtomData {

			int D_ATOM_NO;
			int D_ABS_RES;  
			std::string D_RES_NUM;  
			double D_X;
			double D_Y;
			double D_Z;
			double D_CHARGE;   
			char D_RES_NAME[sizeof(char) * (L_RES_NAME+1)];
			char D_ATOM_NAME [sizeof(char)*(L_ATOM_NAME+1)];
			char D_CHAIN_ID[sizeof(char)*(L_CHAIN_ID +1)];

			AtomData() {
				clear();
			};

			AtomData(const AtomData &_at){
				D_ATOM_NO = _at.D_ATOM_NO;
				D_ABS_RES = _at.D_ABS_RES;  
				D_RES_NUM = _at.D_RES_NUM;  
				D_X  =  _at.D_X;
				D_Y  =  _at.D_Y;
				D_Z  =  _at.D_Z;
				D_CHARGE = _at.D_CHARGE;   
				strcpy(D_RES_NAME,_at.D_RES_NAME);
				strcpy(D_ATOM_NAME,_at.D_ATOM_NAME);
				strcpy(D_CHAIN_ID,_at.D_CHAIN_ID);
			}
			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_RES_NAME   ,"                           ",L_RES_NAME);     
				strncpy(D_ATOM_NAME     ,"                           ",L_ATOM_NAME);      
				strncpy(D_CHAIN_ID      ,"                           ",L_CHAIN_ID);

				// Terminate c-style std::strings properly..
				D_ATOM_NAME[L_ATOM_NAME]           = '\0';
				D_RES_NAME[L_RES_NAME]             = '\0';      
				D_CHAIN_ID[L_CHAIN_ID]             = '\0';      

				// Initialize numeric variables
				D_ATOM_NO   = 1;
				D_ABS_RES   = 1;        
				D_RES_NUM   = "1";        
				D_X         = 0.0;       
				D_Y         = 0.0;       
				D_Z         = 0.0;       
				D_CHARGE    = 0.0;     
			};
		};
		static AtomData parseAtomLine(const std::string &_crdAtomLine);
		static AtomData createAtomData(const Atom &_at, unsigned int _atomNum=1, unsigned int _absres=1);
		static AtomData createAtomData(unsigned int _atomnum, unsigned int _absres, std::string _resname, std::string _atomname, Real &_x, Real &_y, Real &_z, std::string _chainid, int _resnum, double _charge);
		static std::string createAtomLine(const AtomData &ad, unsigned int _atomNum=1, unsigned int _absres=1);
	

	protected:
	private:


};
}

#endif

