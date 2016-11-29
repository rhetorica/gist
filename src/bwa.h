/*
 	Gist bwa.h

 		header for BWA wrapper

    This file is part of Gist.

    Gist is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gist is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gist.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BWA_H_
#define BWA_H_

#include <string>
#include "Class.h"

void parseSAMs(std::string& filename);

void runBWA(vector<string>* class_names, size_t threadcount, std::string& filename);
void indexBWA(string& c);

#endif /* BWA_H_ */
