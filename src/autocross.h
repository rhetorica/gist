/*
 	Gist autocross.h

 		headers for neural network

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

#ifndef AUTOCROSS_H_
#define AUTOCROSS_H_

#include "main.h"

using namespace std;

#include <string>

extern bool** readlabels;
extern number** autocross_MCW;

void autocross_load_weights(string Wfn);

void autocross_save_weights(string Wfn);

void autocross_load_labels_from_fasta();
	// assume genome filenames up to first period start the gene name up to first period
	// (e.g. NCBI accession numbers)

void autocross_train();

void autocross_prep();

void autocross_cleanup();

extern bool autocross_generate;
extern bool autocross_use;

extern int permissable_taxrank;

#endif /* AUTOCROSS_H_ */
