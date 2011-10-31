#ifndef ROTAMERLIBRARYBUILDER_H
#define ROTAMERLIBRARYBUILDER_H

#include <iostream>

#include "Residue.h"
#include "CartesianGeometry.h"
#include "RotamerLibrary.h"


/***********************************************************************************************
 *  This object is an interface for building a RotamerLibrary object fom other
 *  objects, like a residue (currently implemented) or a whole system (to be
 *  implemented. 

     Usage:
 *   Create a Rotamer Library object (say pRotlib) by reading an existing RotamerLibrary file.
 *   Remove/ Retain the original rotamers in pRotlib  as desired
 *   Pass &pRotlib to the RotamerLibraryBuilder.
 *   use addRotamer(Residue& _res, std::string _libName) to add rotamers to pRotlib

 Eg: RotamerLibrary rotLib;
     rotLib.readFile("original_lib.txt"); // open a library as a template
     rotLib.removeAllConformations(); // clear the rotamers
     
     RotamerLibraryBuilder rBuild(&rotLib);
     
     if (sys.residueExists("A,73")) {
        Residue & res = sys.getLastFoundResidue();
        rBuild.addResidue(res, "A_73");
     }
     
     rotLib.writeFile("newLib.txt");

 Note: 1) Use methods in the RotamerLibrary object to read/write libraries or remove rotamers etc.,
       2) Be aware that this object modifies the RotamerLibrary object that gets passed to it.
       3) This object can be made to work on a different RotamerLibrary object using the SetRotamerLibrary method.
 ***********************************************************************************************/

namespace MSL { 
class RotamerLibraryBuilder {
	public:
		RotamerLibraryBuilder();
		RotamerLibraryBuilder(RotamerLibrary * _pRotlib); // rotamers are added to *_pRotlib
		RotamerLibraryBuilder(const RotamerLibraryBuilder & _rotLibbuild); // WARNING: Both "this" and _rotLibbuild will work on the same RotamerLibrary object
		~RotamerLibraryBuilder();

		void setRotamerLibrary(RotamerLibrary * _pRotlib);

		RotamerLibrary * getRotamerLibrary() const;

		bool addRotamer(Residue & _res, std::string _copyLibName, std::string _newLibName, unsigned int _bin = 0) ;
		bool addRotamer(Residue& _res, std::string _libName, unsigned int _bin = 0);

	private:

		RotamerLibrary * pRotlib;

};

inline void RotamerLibraryBuilder::setRotamerLibrary(RotamerLibrary * _pRotlib) { pRotlib = _pRotlib; }
inline RotamerLibrary * RotamerLibraryBuilder::getRotamerLibrary() const {return pRotlib;}

}

#endif
