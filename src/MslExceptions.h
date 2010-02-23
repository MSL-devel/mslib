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

#ifndef MSLEXCEPTIONS_H
#define MSLEXCEPTIONS_H
/*
Note:
Standard Library class runtime_error (defined in header file <stdexcept>). 
Class runtime_error derived class of Standard Library class exception (defined in header file <exception>)
is the C++ standard base class for representing runtime errors. 
A typical exception class that derives from the runtime_error class defines only a constructor 
that passes an error-message std::string to the base-class runtime_error constructor

every exception class that derives directly or indirectly from exception contains the virtual function what, 
which returns an exception object¡Çs error message.

Basically, it is cleaner and simpler to inherit std::runtime_error , than excpetion.
With runtime_error you only NEED to call runtime_error(std::string msg) constructor from
your constructor. (no need for fancy virtual function declarations of the 'what' function).

 */

// STL Includes
//#include <exception>
#include <stdexcept>
#include <string>


namespace MSL { 
class ConvertDoubleException : public std::runtime_error {

	public:
		ConvertDoubleException() : std::runtime_error("ConvertDoubleException") {}
		ConvertDoubleException(const std::string &_msg) : std::runtime_error("ConvertDoubleException: "+_msg) {}
};


class ConvertIntException : public std::runtime_error {

	public:
		ConvertIntException() : std::runtime_error("ConvertIntException") {}
		ConvertIntException(const std::string _msg) : std::runtime_error("ConvertIntException: "+_msg) {}

};


class MslSizeException : public std::runtime_error {
	public:
		MslSizeException() : std::runtime_error("MslSizeException") {}
		MslSizeException(const std::string _msg) : std::runtime_error("MslSizeException: "+_msg) {}
};


class MslGeneralException : public std::runtime_error {
	public:
		MslGeneralException() : std::runtime_error("MslGeneralException") {}
		MslGeneralException(const std::string _msg) : std::runtime_error("MslGeneralException: "+_msg) {}
};

}

#endif
