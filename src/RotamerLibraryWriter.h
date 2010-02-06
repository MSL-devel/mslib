#ifndef ROTAMERLIBRARYWRITER_H
#define ROTAMERLIBRARYWRITER_H
/*

 */

//MSL Includes
#include "Writer.h"
#include "RotamerLibrary.h"
#include "MslTools.h"

// STL Includes
#include <vector>
#include <set>
#include <map>
#include <iostream>
using namespace std;

/**
 * This class will provide an object which is able
 * to write in and interpret RotamerLibrary files.
 */
class RotamerLibraryWriter : public Writer {

	public:
		// Constructors/Destructors
		RotamerLibraryWriter();
		RotamerLibraryWriter(const string &_filename);
		virtual ~RotamerLibraryWriter();

		bool write(RotamerLibrary * _rotlib, string _charmm = "CHARMMPAR 22 27");
		bool writeResidue(const string &_res, const string &_libName, RotamerLibrary *_rotlib);
		bool writeLibrary(const string &_libname, RotamerLibrary *_rotlib);
		bool open();
		bool open(const string &_filename); // There is a default implementation
		bool open(const string &_filename, int mode); // There is a default implementation
		bool open(stringstream &_ss);
		void close();
		void writeREMARKS();
	protected:		
	private:
		void setup();
		string createInitLine( vector<string> _atoms);
		string createInternalCoorDefinition( vector<RotamerLibrary::InternalCoorDefi> icDefis);
		string createInternalCoor(vector<vector<double > > coor);
		//static const int resList_size = 22;
		//static std::string resList[17];
		//static std::string resList[resList_size];
};

//Inlines go HERE
	inline RotamerLibraryWriter::RotamerLibraryWriter() : Writer() { }
	inline RotamerLibraryWriter::RotamerLibraryWriter(const string &_filename) : Writer(_filename) {}
	inline RotamerLibraryWriter::~RotamerLibraryWriter() {}
	inline bool RotamerLibraryWriter::open() {bool success = Writer::open(); return success;}
	inline bool RotamerLibraryWriter::open(const string &_filename) {bool success = Writer::open(_filename); return success;}
	inline bool RotamerLibraryWriter::open(const string &_filename, int mode) {bool success = Writer::open(_filename, mode); return success;}
	inline bool RotamerLibraryWriter::open(stringstream &_ss) {bool success = Writer::open(_ss); return success;}
	inline void RotamerLibraryWriter::close() { Writer::close(); }

	inline void RotamerLibraryWriter::writeREMARKS() {};
	

#endif
