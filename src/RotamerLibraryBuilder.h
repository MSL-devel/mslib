#ifndef ROTAMERLIBRARYBUILDER_H
#define ROTAMERLIBRARYBUILDER_H

#include <iostream>

#include "Residue.h"
#include "CartesianGeometry.h"
#include "RotamerLibrary.h"


/***************************************************************************
 *  This object is an interface for building a RotamerLibrary object fom other
 *  objects, like a residue (currently implemented) or a whole system (to be
 *  implemented. 
 ***************************************************************************/

namespace MSL { 
class RotamerLibraryBuilder {
	public:
		RotamerLibraryBuilder();
		RotamerLibraryBuilder(RotamerLibrary * _pRotlib);
		RotamerLibraryBuilder(const RotamerLibraryBuilder & _rotLibbuilder);
		~RotamerLibraryBuilder();

		void setRotamerLibrary(RotamerLibrary * _pRotlib);

		RotamerLibrary * getRotamerLibrary() const;

		bool addRotamer(Residue& _res, std::string _libName);

	private:

		RotamerLibrary * pRotlib;

};

inline void RotamerLibraryBuilder::setRotamerLibrary(RotamerLibrary * _pRotlib) { pRotlib = _pRotlib; }
inline RotamerLibrary * RotamerLibraryBuilder::getRotamerLibrary() const {return pRotlib;}

}

#endif
