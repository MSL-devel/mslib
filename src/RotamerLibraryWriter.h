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

/**
 * This class will provide an object which is able
 * to write in and interpret RotamerLibrary files.
 */
namespace MSL { 
class RotamerLibraryWriter : public Writer {

	public:
		// Constructors/Destructors
		RotamerLibraryWriter();
		RotamerLibraryWriter(const std::string &_filename);
		virtual ~RotamerLibraryWriter();

		bool write(RotamerLibrary * _rotlib, std::string _charmm = "CHARMMPAR 22 27");
		bool open();
		bool open(const std::string &_filename); // There is a default implementation
		bool open(const std::string &_filename, int mode); // There is a default implementation
		bool open(std::stringstream &_ss);
		void close();
		void writeREMARKS();
	protected:		
	private:
		void setup();
		std::string createInitLine( std::vector<std::string> _atoms);
		std::string createInternalCoorDefinition( std::vector<RotamerLibrary::InternalCoorDefi> icDefis);
		bool writeResidue(const std::string &_res, const std::string &_libName, RotamerLibrary *_rotlib);
		bool writeLibrary(const std::string &_libname, RotamerLibrary *_rotlib);

		bool writeLevelInformation(RotamerLibrary *_rotlib);
		std::string createInternalCoor(std::vector<std::vector<double > > coor);
		//static const int resList_size = 22;
		//static std::string resList[17];
		//static std::string resList[resList_size];
};

//Inlines go HERE
	inline RotamerLibraryWriter::RotamerLibraryWriter() : Writer() { }
	inline RotamerLibraryWriter::RotamerLibraryWriter(const std::string &_filename) : Writer(_filename) {}
	inline RotamerLibraryWriter::~RotamerLibraryWriter() {}
	inline bool RotamerLibraryWriter::open() {bool success = Writer::open(); return success;}
	inline bool RotamerLibraryWriter::open(const std::string &_filename) {bool success = Writer::open(_filename); return success;}
	inline bool RotamerLibraryWriter::open(const std::string &_filename, int mode) {bool success = Writer::open(_filename, mode); return success;}
	inline bool RotamerLibraryWriter::open(std::stringstream &_ss) {bool success = Writer::open(_ss); return success;}
	inline void RotamerLibraryWriter::close() { Writer::close(); }

	inline void RotamerLibraryWriter::writeREMARKS() {};
	

}

#endif
