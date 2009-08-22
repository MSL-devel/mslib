/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#ifndef PDBFORMAT_H
#define PDBFORMAT_H
/*
  PDB Formating information goes here.  As well as a simple structure for containing a single line
  of a PDB.

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


Example PDB lines (Put standard + any strange ones here)

          1         2         3         4         5         6         7                     
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM     15  CD  LYS A   1      17.517  50.534  29.567  1.00 29.97           C
HETATM 1101  O   HOH   117      30.011  51.798   1.270  1.00 43.26           O




 */

//MSL Includes
#include "MslTools.h"
#include "Atom.h"

//STL Includes
#include <string>
#include <iostream>
#include <cstring>

using namespace std;

class PDBFormat {
	public:

		enum AtomField { 

			// Start Record Position
			S_RECORD_NAME    =  0,
			S_SERIAL         =  6,
			S_ATOM_NAME      = 12,
			S_ALT_LOC        = 16,
			S_RES_NAME       = 17,
			S_CHAIN_ID       = 21,
			S_RES_SEQ        = 22,
			S_I_CODE         = 26,
			S_X              = 30,
			S_Y              = 38,
			S_Z              = 46,
			S_OCCUP          = 54,
			S_TEMP_FACT      = 60,
			S_SEG_ID         = 72,
			S_ELEMENT_SYMBOL = 76,
			S_CHARGE         = 78,

			// End Record Position
			E_RECORD_NAME    =  5,
			E_SERIAL         = 10,
			E_ATOM_NAME      = 15,
			E_ALT_LOC        = 16,
			E_RES_NAME       = 19,
			E_CHAIN_ID       = 21,
			E_RES_SEQ        = 25,
			E_I_CODE         = 26,
			E_X              = 37,
			E_Y              = 45,
			E_Z              = 53,
			E_OCCUP          = 59,
			E_TEMP_FACT      = 65,
			E_SEG_ID         = 75,
			E_ELEMENT_SYMBOL = 77,
			E_CHARGE         = 79,


			// Length of Record
			L_RECORD_NAME    =  6,
			L_SERIAL         =  5,
			L_ATOM_NAME      =  4,
			L_ALT_LOC        =  1,
			L_RES_NAME       =  3,
			L_CHAIN_ID       =  1,
			L_RES_SEQ        =  4,
			L_I_CODE         =  1,
			L_X              =  8,
			L_Y              =  8,
			L_Z              =  8,
			L_OCCUP          =  6,
			L_TEMP_FACT      =  6,
			L_SEG_ID         =  4,
			L_ELEMENT_SYMBOL =  2,
			L_CHARGE         =  2,


			// Justification within record 0 = right-justified, 1 = left-justified.
			// TODO: Make this correct!
			J_RECORD_NAME    = 1,
			J_SERIAL         = 1,
			J_ATOM_NAME      = 1,
			J_ALT_LOC        = 1,
			J_RES_NAME       = 1,
			J_CHAIN_ID       = 1,
			J_RES_SEQ        = 1,
			J_I_CODE         = 1,
			J_X              = 1,
			J_Y              = 1,
			J_Z              = 1,
			J_OCCUP          = 1,
			J_TEMP_FACT      = 1,
			J_SEG_ID         = 1,
			J_ELEMENT_SYMBOL = 1,
			J_CHARGE         = 1

			
		};

		struct AtomData {
			char   D_RECORD_NAME[sizeof(char)*(L_RECORD_NAME+1)];
			int    D_SERIAL;
			char   D_ATOM_NAME[sizeof(char)*(L_ATOM_NAME+1)];
			char   D_ALT_LOC[sizeof(char)*(L_ALT_LOC+1)];
			char   D_RES_NAME[sizeof(char)*(L_RES_NAME+1)];
			char   D_CHAIN_ID[sizeof(char)*(L_CHAIN_ID+1)];
			int    D_RES_SEQ;
			char   D_I_CODE[sizeof(char)*(L_I_CODE+1)];
			double D_X;
			double D_Y;
			double D_Z;
			double D_OCCUP;
			double D_TEMP_FACT;
			char   D_SEG_ID[sizeof(char)*(L_SEG_ID+1)];
			char   D_ELEMENT_SYMBOL[sizeof(char)*(L_ELEMENT_SYMBOL+1)];
			double D_CHARGE;

			AtomData() {
				clear();
			};

			AtomData(const AtomData &_at){

				strcpy(D_RECORD_NAME,_at.D_RECORD_NAME);
				strcpy(D_ATOM_NAME,_at.D_ATOM_NAME);
				strcpy(D_ALT_LOC,_at.D_ALT_LOC);
				strcpy(D_RES_NAME,_at.D_RES_NAME);
				strcpy(D_CHAIN_ID,_at.D_CHAIN_ID);
				strcpy(D_I_CODE,_at.D_I_CODE);
				strcpy(D_SEG_ID,_at.D_SEG_ID);
				strcpy(D_ELEMENT_SYMBOL,_at.D_ELEMENT_SYMBOL);

				D_SERIAL = _at.D_SERIAL ;       
				D_RES_SEQ = _at.D_RES_SEQ ;      
				D_X = _at.D_X ;	       
				D_Y = _at.D_Y ;	       
				D_Z = _at.D_Z ;	       
				D_OCCUP = _at.D_OCCUP ;	       
				D_TEMP_FACT = _at.D_TEMP_FACT ;    


			}
			void clear() {

				// Create a blank string of size L_FIELD_NAME
				strncpy(D_RECORD_NAME   ,"                           ",L_RECORD_NAME);     
				strncpy(D_ATOM_NAME     ,"                           ",L_ATOM_NAME);      
				strncpy(D_ALT_LOC       ,"                           ",L_ALT_LOC);
				strncpy(D_RES_NAME      ,"                           ",L_RES_NAME);
				strncpy(D_CHAIN_ID      ,"                           ",L_CHAIN_ID);
				strncpy(D_I_CODE        ,"                           ",L_I_CODE);
				strncpy(D_SEG_ID        ,"                           ",L_SEG_ID);
				strncpy(D_ELEMENT_SYMBOL,"                           ",L_ELEMENT_SYMBOL); 	

				// Terminate c-style strings properly..
				D_RECORD_NAME[L_RECORD_NAME]       = '\0';
				D_ATOM_NAME[L_ATOM_NAME]           = '\0';
				D_ALT_LOC[L_ALT_LOC]               = '\0';       
				D_RES_NAME[L_RES_NAME]             = '\0';      
				D_CHAIN_ID[L_CHAIN_ID]             = '\0';      
				D_I_CODE[L_I_CODE]                 = '\0';        
				D_SEG_ID[L_SEG_ID]                 = '\0';        
				D_ELEMENT_SYMBOL[L_ELEMENT_SYMBOL] = '\0';

				// Initialize numeric variables
				D_SERIAL    = 0;
				D_RES_SEQ   = 0;        
				D_X         = 0.0;       
				D_Y         = 0.0;       
				D_Z         = 0.0;       
				D_OCCUP     = 0.0;       
				D_TEMP_FACT = 0.0;     
				D_CHARGE    = 0.0;     
			};
		};


		static AtomData parseAtomLine(const string &_pdbAtomLine);
		static AtomData createAtomData(const Atom &_at);
		static AtomData createAtomData(string _resName, Real &_x, Real &_y, Real &_z, string _element);
		static string createAtomLine(const AtomData &ad);
		static string createTerLine(const AtomData &ad);


	protected:
	private:


};
#endif
