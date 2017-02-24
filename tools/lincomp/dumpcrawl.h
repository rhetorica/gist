/*
 * dumpcrawl.h
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#ifndef DUMPCRAWL_H_
#define DUMPCRAWL_H_

#include <string>
#include "main.h"

using namespace std;

extern size_t sz1, sz2, sz3;

string get_tax_id(const char* gi, size_t num, size_t lens, size_t lo, size_t hi);
string get_parent_node(const char* taxid, size_t num, size_t lens, size_t lo, size_t hi);
string get_name_row_up(const char* taxid, size_t lens, size_t q);
string get_name_row_down(const char* taxid, size_t lens, size_t q, size_t hi);
string get_name_row(const char* taxid, size_t num, size_t lens, size_t lo, size_t hi);
void load_db();

string gi_from_fasta(string header);

#endif /* DUMPCRAWL_H_ */
