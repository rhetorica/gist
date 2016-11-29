/*
 * SparseTable.h
 *
 *  Created on: 2015-10-09
 *      Author: samantha
 */

#ifndef SPARSETABLE_H_
#define SPARSETABLE_H_

#include "main.h"
#include <algorithm>
#include <map>

class SparseTable {
public:
	SparseTable(char* source, int prefix_size, tableType tt);
	void add_sequence(char* source);
	virtual ~SparseTable();
	number get(unsigned long int key);
	void increment(unsigned long int key);
	void flatten();

	bool smoothflag;

	//Table& operator+=(Table& rhs);
	//Table* operator*(Table& rhs);
	//Table* operator-(Table& rhs);
	//Table* operator*(number amount);
	//Table& operator/=(number amount);
	//Table* powt(number exponent);

	number logmultiprob(SparseTable* genome);

private:
	size_t bits, tt, word_length, length, base;

	std::map<unsigned long int, number>* hashmap;

	bool flattened;

	long unsigned int** flat_keys;
	number** flat_values;
	number* get_position(long unsigned int key);
	number* create_position(long unsigned int key);

	void* tree;
};

#endif /* SPARSETABLE_H_ */
