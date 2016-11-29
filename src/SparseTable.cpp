/*
 * SparseTable.cpp
 *
 *  Created on: 2015-10-09
 *      Author: samantha
 */

#include "SparseTable.h"
#include <sstream>
#include <cstring>

number SparseTable::logmultiprob(SparseTable* genome) {
	number sum = 0;

//	const number pi = atan(1)*4;

	/*
	 * the core method here is based on studying the intersection of dimensions
	 * it quickly became apparent that all that really matters is the tally of
	 * good hits - a read with 98% identity to a given segment of two chromosomes
	 * isn't inherently more likely to be from the shorter one!
	 *
	 */

	// map<long unsigned int, number> fed = map<long unsigned int, number>(*(genome->hashmap));

	number reg = 1 / (number)genome->length;

	for (std::map<long unsigned int, number>::iterator it=hashmap->begin(); it!=hashmap->end(); ++it) {
		if(genome->hashmap->count(it->first)) {
			sum += log(hashmap->at(it->first) / (number)length + reg)
				   -
				   log(it->second / (number)genome->length + reg);
		}
	}

	// perform actual calculation

	return sum;
}

number* SparseTable::get_position(unsigned long int key) {
	/*if(flattened) {
		cerr << "Implement me, Seymour.";
		exit(1);
	} else {
		size_t pos = word_length - 1;

		long unsigned int bitmask = (2 << bits) - 1;
		long unsigned int shiftedmask = bitmask << (pos * bits);
		char here = key & shiftedmask >> (pos * bits);
	}*/

	return NULL; // unneeded for map version
}

number* SparseTable::create_position(unsigned long int key) {
	/*if(flattened) {
		cerr << "Position creation is only possible in unflattened mode.";
		exit(1);
	} else {

	}*/

	hashmap->insert(pair<unsigned long int, number>(key, 0.0));
	return NULL; // unneeded for map version
}

number SparseTable::get(unsigned long int key) {
	/*number answer;

	number* pos = this->get_position(key);

	if(pos)
		answer = *pos;

	if(smoothflag)
		answer += 1 / this->length; // not laplace smoothing

	return answer;*/

	number answer = 0;

	if(hashmap->count(key) != 0)
		answer = (number)hashmap->at(key);

	if(smoothflag)
		answer++; // pseudo-laplace smoothing

	return answer / this->length;
}

void SparseTable::increment(unsigned long int key) {
	/*number* pos = this->get_position(key);
	if(pos) {
		pos = this->create_position(key);
		(*pos)++;
	}*/

	if(hashmap->count(key) == 0)
		hashmap->insert(pair<unsigned long int, number>(key, 1.0));
	else
		++(hashmap->at(key));
}

void SparseTable::add_sequence(char* source) {
	size_t srclen = strlen(source);

	size_t prefix_size = word_length - 1;

	this->length += srclen - prefix_size;

	unsigned long int mask = pow(base, prefix_size) - 1;
	unsigned long int prefix = 0;
	int skips = 0;
	for(size_t i = 0; i < srclen - 1; i++) {
		char newfix, focus;
		if(tt == NUCLEOTIDE) {
			newfix = twobits(source[i], i);
			focus = twobits(source[i+1], i+1);
		} else {
			newfix = fivebits(source[i]);
			focus = fivebits(source[i+1]);
			if(focus == -1) skips = prefix_size + 1;
		}
		prefix = ((prefix << bits) + newfix) & mask;
		if(skips != 0) {
			skips = skips - 1;
		} else if(i >= (size_t)prefix_size - 1) {
			this->increment((prefix << bits) + focus);
		}
	}
}

SparseTable::SparseTable(char* source, int prefix_size, tableType tt) {
	if(debugdump) cerr << "my address is " << this << "." << endl;

	flattened = false;

	this->tt = tt;

	this->word_length = prefix_size + 1;
	if(this->word_length == 1) {
		cerr << "Sparse Table with word length 1 is illegal." << endl;
		throw(string("I think you really should use grep for that."));
	}

	if(debugdump) cerr << "word size: " << word_length << endl;

	if(tt == NUCLEOTIDE) {
		bits = 2;
		base = 4;
	} else {
		bits = 5;
		base = 20;
	}

	size_t srclen = strlen(source);

	this->length = srclen - prefix_size;

	smoothflag = false;

	if(strcmp(source, "smooth") == 0) {
		smoothflag = true;
		this->length = 0; // todo: Table has this: pow(cols, prefix_size + 1);
	}

	if(this->length < 0) this->length = 0;

	this->hashmap = new map<unsigned long int, number>();

	this->tree = NULL; // new void*[base];
	flat_keys = NULL;
	flat_values = NULL;

	if(this->length < 1 || smoothflag) {
		return;
	}

	unsigned long int mask = pow(base, prefix_size) - 1;
	unsigned long int prefix = 0;
	int skips = 0;
	for(size_t i = 0; i < srclen - 1; i++) {
		char newfix, focus;
		if(tt == NUCLEOTIDE) {
			newfix = twobits(source[i], i);
			focus = twobits(source[i+1], i+1);
		} else {
			newfix = fivebits(source[i]);
			focus = fivebits(source[i+1]);
			if(focus == -1) skips = prefix_size + 1;
		}
		prefix = ((prefix << bits) + newfix) & mask;
		if(skips != 0) {
			skips = skips - 1;
		} else if(i >= (size_t)prefix_size - 1) {
			this->increment((prefix << bits) + focus);
		}
	}
}

SparseTable::~SparseTable() {
	this->hashmap->clear();
	delete this->hashmap;
}

