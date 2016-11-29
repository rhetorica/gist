/*
 	Gist build.h

 		headers for class construction thread management

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

#ifndef BUILD_H_
#define BUILD_H_

#include "Class.h"
#include <string>

void oneclass(vector<string>* class_filenames, Class** classes, size_t i);
void classnothreads(vector<string>* class_filenames, Class** classes);
void classthreads(vector<string>* class_filenames, Class** classes);

Class* classthread(string filename);

#endif /* BUILD_H_ */
