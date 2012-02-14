/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

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

#ifndef FORMATCONVERTER_H
#define FORMATCONVERTER_H

#include "System.h"


/*! \brief  This Object converts charmm variables from and to pdb variables
 *  (atom and residue names )
 *
 */

using namespace std;

class FormatConverter {
	public:

		FormatConverter();
		FormatConverter(System * _pSys);

		~FormatConverter();

		void setPdbFromCharmm(string charmmVersion);
		void setPdbFromCharmm(System * _pSys, string _charmmVersion);
		void setPdbFromCharmm(System * _pSys, string _charmmVersion, bool _Nterminal, bool _Cterminal); 

		void setCharmmFromPdb(string charmmVersion);
		void setCharmmFromPdb(System * _pSys, string _charmmVersion);
		void setCharmmFromPdb(System * _pSys, string _charmmVersion, bool _Nterminal, bool _Cterminal); 

	private:
		void setCharmmFromPdb(Chain & _rChain, string _charmmVersion);
		void setCharmmFromPdb(Residue & _rRes, string _charmmVersion, bool _Nterminal, bool _Cterminal);
		void setCharmmFromPdb(Atom & _rAtom, string _charmmVersion, bool _Nterminal, bool _Cterminal);

		void setPdbFromCharmm(Chain & _rChain, string _charmmVersion);
		void setPdbFromCharmm(Residue & _rRes, string _charmmVersion, bool _Nterminal, bool _Cterminal);
		void setPdbFromCharmm(Atom & _rAtom, string _charmmVersion, bool _Nterminal, bool _Cterminal);
		
		System * pSys;

};

#endif

