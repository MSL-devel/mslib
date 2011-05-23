/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#ifndef SYSENV_H
#define SYSENV_H

// STL Includes
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>


namespace MSL { 


	class SysEnv {
		public:
			SysEnv();
			~SysEnv();

			bool addEnvVariable(std::string &_var);
			bool setEnv(std::string &_var,std::string &_value);
			bool isDefined(std::string &_var);
			std::string getEnv(const std::string &_var);
			
		private:
			void setup();

			std::map<std::string,std::string> env;
			std::map<std::string,std::string>::iterator envIt;
			std::string defaultString;
			std::string undefinedString;
			
		


	};
}
#endif
