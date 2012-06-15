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

#include "CRDFormat.h"

using namespace MSL;
using namespace std;

CRDFormat::AtomData CRDFormat::parseAtomLine(const string &_crdAtomLine){
	AtomData atom;

	int lineLength = _crdAtomLine.size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_ATOM_NO)    atom.D_ATOM_NO         = MslTools::toInt(MslTools::trim(_crdAtomLine.substr(S_ATOM_NO,L_ATOM_NO)), "Atom number(ATOM_NO) problem");
		if (lineLength >= E_ABS_RES)    atom.D_ABS_RES         = MslTools::toInt(MslTools::trim(_crdAtomLine.substr(S_ABS_RES,L_ABS_RES)), "Residue Number(ABS_RES) problem");
		
		if (lineLength >= E_RES_NAME)        strcpy(atom.D_RES_NAME,        MslTools::toUpper(MslTools::trim(_crdAtomLine.substr(S_RES_NAME,L_RES_NAME))).c_str());
		if (lineLength >= E_ATOM_NAME)       strcpy(atom.D_ATOM_NAME,       MslTools::toUpper(MslTools::trim(_crdAtomLine.substr(S_ATOM_NAME, L_ATOM_NAME))).c_str());


		if (lineLength >= E_X)         atom.D_X              = MslTools::toDouble(MslTools::trim(_crdAtomLine.substr(S_X,L_X)), "X-coord problem");
		if (lineLength >= E_Y)         atom.D_Y              = MslTools::toDouble(MslTools::trim(_crdAtomLine.substr(S_Y,L_Y)), "Y-coord problem");
		if (lineLength >= E_Z)         atom.D_Z              = MslTools::toDouble(MslTools::trim(_crdAtomLine.substr(S_Z,L_Z)), "Z-coord problem");


		if (lineLength >= E_CHAIN_ID)       strcpy(atom.D_CHAIN_ID,       MslTools::toUpper(MslTools::trim(_crdAtomLine.substr(S_CHAIN_ID, L_CHAIN_ID))).c_str());

		if (lineLength >= E_RES_NUM)   atom.D_RES_NUM        = MslTools::trim(_crdAtomLine.substr(S_RES_NUM,L_RES_NUM), "Residue Number problem");

		// Charge is  a problem
		if (lineLength >= E_CHARGE && MslTools::trim(_crdAtomLine.substr(S_CHARGE, L_CHARGE)) != "")    atom.D_CHARGE         = MslTools::toDouble(MslTools::trim(_crdAtomLine.substr(S_CHARGE, L_CHARGE)), "CHARGE problem");

	} catch(exception &e){

		cerr << "ERROR 34919 CRDFormat parseAtomLine "<<e.what()<<endl;
		exit(34918);

	}

	return atom;
}


CRDFormat::AtomData CRDFormat::createAtomData(unsigned int _atomnum, unsigned int _absres, string _resname, string _atomname, Real &_x, Real &_y, Real &_z, string _chainid, int _resnum, double _charge){
	Atom a(_resname, _x,_y,_z);
	a.setCharge(_charge);
	a.setChainId(_chainid);
	a.setResidueNumber(_resnum);
	a.setName(_atomname);
	return createAtomData(a, _atomnum, _absres);
		
}
CRDFormat::AtomData CRDFormat::createAtomData(const Atom &_at, unsigned int _atomNum, unsigned int _absres){
	CRDFormat::AtomData atom;


// Dont know Insertion Code
//	strncpy(atom.D_I_CODE, _at.getResidueIcode().c_str(), CRDFormat::L_I_CODE);

	atom.D_ATOM_NO  = _atomNum;
	atom.D_ABS_RES = _absres;
	atom.D_RES_NUM = MslTools::intToString(_at.getResidueNumber())+_at.getResidueIcode();

	atom.D_X = _at.getX();
	atom.D_Y = _at.getY();
	atom.D_Z = _at.getZ();

	strncpy(atom.D_RES_NAME, _at.getResidueName().c_str(), CRDFormat::L_RES_NAME);
	strncpy(atom.D_ATOM_NAME, _at.getName().c_str(), CRDFormat::L_ATOM_NAME);
	strncpy(atom.D_CHAIN_ID, _at.getChainId().c_str(), CRDFormat::L_CHAIN_ID);

	atom.D_CHARGE = _at.getCharge();

	return atom;

}

string CRDFormat::createAtomLine(const CRDFormat::AtomData &ad, unsigned int _atomNum, unsigned int _absres){

	/*
	          1         2         3         4         5         6         7                     
	01234567890123456789012345678901234567890123456789012345678901234567890123456789
	<---><--->x<-->x<--><--------><--------><-------->x<-->x<--><-------->
	    1    1 ALA  N      2.14273   1.32845   0.00000 A    1     -0.30000
	    2    1 ALA  HT1    3.17887   1.23889   0.00000 A    1      0.33000
	*/
	char c [1000];
	sprintf(c, "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4s%10.5f", _atomNum, _absres,ad.D_RES_NAME,ad.D_ATOM_NAME,ad.D_X, ad.D_Y, ad.D_Z,ad.D_CHAIN_ID,ad.D_RES_NUM.c_str(),ad.D_CHARGE);
	
	return (string)c;

}
