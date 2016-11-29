/*
 	Gist score.h

 		header for score calculation

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

#ifndef S_CLASSIFY_H_
#define S_CLASSIFY_H_

#include "Table.h"

void score_read(stringstream& lout, size_t j, size_t k);
void score(size_t start, size_t end);
void scorethreads(size_t numthreads);


#endif /* S_CLASSIFY_H_ */
