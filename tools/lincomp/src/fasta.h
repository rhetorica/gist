/*
 * fasta.h
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#ifndef FASTA_H_
#define FASTA_H_

#include <string>
#include "main.h"

void scrub(string& b);
int readFasta(string filename, char**& records);
int load_class_id_from_filename(const char* h);
int load_class_id_from_fasta_label(const char* h);

#endif /* FASTA_H_ */
