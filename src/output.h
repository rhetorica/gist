/*
 	Gist output.h

 		header for output and taxon inclusion calculations

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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "main.h"
#include <cmath>

#include <map>

extern std::map<int,number>* taxscores;

number calc_taxon_inclusion(number* return_scores, unsigned long taxid, TAXON_RANKS rank, size_t pcount, size_t j, number bloom);

void write_output();

#endif /* OUTPUT_H_ */
