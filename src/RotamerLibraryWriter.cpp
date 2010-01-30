#include "RotamerLibraryWriter.h"

//string RotamerLibraryWriter::resList[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HSD","HSE","HSP","ILE","LEU","LYS","MET","PHE","SER","THR","TRP","TYR","VAL"}; // update the RotamerLibraryWriter::resList_size variable when changing this.

bool RotamerLibraryWriter::write(RotamerLibrary & _rotlib, string _charmm) {
	vector<string> libs = _rotlib.getLibraryNames();
	//cout << "UUUUUUUUUUU Libs size = " << libs.size() << endl;
	for (int i = 0; i < libs.size(); i++) {
		string 	line = "LIBRARY " + libs[i] + "\n" + _charmm + "\n";
		writeln(line);
		writeLibrary(libs[i],_rotlib); 
	}
	return true;
}

bool RotamerLibraryWriter::writeLibrary(const string &_libName, RotamerLibrary &_rotlib) {
	if(_rotlib.libraryExists(_libName) == false) {
		return false;
	}
	set<string> resList = _rotlib.getResList(_libName); 
	for(set<string>::iterator i = resList.begin(); i != resList.end(); i++) {
//		writeResidue(resList.at(i),_libName,_rotlib);
		writeResidue((*i),_libName,_rotlib);
	}
	return true;
}


bool RotamerLibraryWriter::writeResidue(const string &_res, const string &_libName, RotamerLibrary &_rotlib) {
	
	if(_rotlib.residueExists(_libName,_res) == false) {
		return false;
	}

	string Residue = "RESI ";
	Residue +=  (_res + "\n");
	//cout << "UUUU write Residue:" << Residue;
	Residue += _rotlib.getInitAtomsLine(_libName, _res) + "\n";


	vector<string> defLines =  _rotlib.getInternalCoorDefinitionLines(_libName, _res);

	for(vector<string>::iterator line = defLines.begin(); line != defLines.end(); line++) {
		Residue += ((*line) + "\n");
	}

	vector<string> coorLines =  _rotlib.getAllInternalCoorLines(_libName, _res);

	for(vector<string>::iterator line = coorLines.begin(); line != coorLines.end(); line++) {
		Residue += ((*line) + "\n");
	}

	Residue += "\n";

	
	if(Writer::write(Residue) == false) {
		return false;
	} else {
		return true;
	}	
}












