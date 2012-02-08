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

