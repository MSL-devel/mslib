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


#include "File.h"

using namespace MSL;
using namespace std;



File::File() {
	string filename = "";
	init(filename, 0);
}

File::File(stringstream &_ss){
	init("", 0);
	fileName        = "string";
	stringStreamPtr = &_ss;
	fileHandler = stringstyle;
}


File::File(const string &_filename, int _mode){
	init(_filename, _mode);
}

File::~File(){
    if(is_open())
      close();
}
void File::init(const string &_filename, int _mode){
	fileName     = _filename;
	fileOpenMode = (File::fileModes)_mode;
	fileHandler = cppstyle;
	filePtr = NULL;
	stringStreamPtr = NULL;
	close();
}

bool File::open(){

	// Insure that the file is closed first
	close();
	string mode = "";
	std::ios::openmode omode = std::ios::in;
	switch (fileHandler) {
		case cstyle:
			switch (fileOpenMode) {
				case read:
					mode = "r";
					break;
				case write:
					mode = "w";
					break;
				case append:
					mode = "a";
					break;
			}	
	
			filePtr = fopen(fileName.c_str(), mode.c_str());
			break;
		case cppstyle:
			switch (fileOpenMode) {
				case read:
					omode = std::ios::in;
					break;
				case write:
					omode = std::ios::out;
					break;
				case append:
					omode = std::ios::out | std::ios::app;
					break;
			}	
			fileStream.open(fileName.c_str(), omode);

			return !fileStream.fail();					
			break;
		case stringstyle:
			setFileName("string");
			break;
	       
	}


	return true;
}

bool File::open(const string &_filename){
	setFileName(_filename);
	return open();
}

bool File::open(const string &_filename, int mode){
	setFileName(_filename);
        fileOpenMode = (fileModes)mode;
	return open();
}

bool File::open(stringstream &_ss){
	stringStreamPtr = &_ss;
	return open();
}

void File::close(){

	switch (fileHandler) {
		case cstyle:
			if (filePtr != NULL){
				fclose(filePtr);
				filePtr = NULL;
			}
			break;
		case cppstyle:

			if (fileStream.is_open()){
				fileStream.close();
			}
			fileStream.clear();
			break;
		case stringstyle:
			stringStreamPtr->clear();
			break;
	       
	}


}


bool File::doesFileExist(){
	bool fileExists = open();
	close();
	return fileExists;
}


bool File::is_open() {
    bool isOpen = false;
    switch (fileHandler) {
		case cstyle:
            isOpen = (filePtr != NULL);
			break;
		case cppstyle:
            isOpen = fileStream.is_open();
			break;
		case stringstyle:
            isOpen = (stringStreamPtr != NULL);
			break;
	}

	return isOpen;
}

bool File::endOfFileTest(){


	switch (fileHandler) {
		case cstyle:
			break;
		case cppstyle:
			if (fileStream.fail() || fileStream.eof()){
				return true;
			}
			break;
		case stringstyle:
			if (stringStreamPtr->fail()){
				return true;
			}
			break;
	}

	return false;
}


void File::setOpenMode(int _mode) { 

	fileOpenMode = (File::fileModes)_mode; 
}
File::fileModes File::getOpenMode() const{ 
	return fileOpenMode; 
}


void File::setFileHandlerType(int _fileHandler) { 

	fileHandler = (File::fileHandlerType)_fileHandler;
}
File::fileHandlerType File::getFileHandlerType() const{ 
	return fileHandler; 
}

void File::copy(const File &_anotherFile){
	
	fileName    = _anotherFile.getFileName();
	remarks     = _anotherFile.getRemarks();
	fileHandler = _anotherFile.getFileHandlerType();
	fileOpenMode= _anotherFile.getOpenMode();
}
void File::clearRemarks(){
	remarks.clear();
}
