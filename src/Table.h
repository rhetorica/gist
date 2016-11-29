/*
 	Gist Table.h

 		header for read tables

 		each Table object represents either a read, a gene, or some
 		statistical generalization thereof (e.g. a mean value)

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

#ifndef TABLE_H_
#define TABLE_H_

#include "main.h"

#define Coord(row, col) ((row) * (cols) + (col))

#define adjustRow(i) do {					\
				size_t j = i, k = 0;		\
				while(j > base) {			\
					j /= base;				\
					k++;					\
				}							\
				if(j > cols - 1) {			\
					i = pow(base, k+1);		\
				}							\
				} while(0)


class Table {
public:
	Table(char*, int, tableType); // for processing sequences
	Table(string*, int, tableType); // for processing sequences
	Table(string& serialized, int, tableType); // for re-importing from a file
	virtual ~Table();
	number** mydata;
	int length;
	int prefix_length;
	number delta(Table*);
	size_t base, rows, cols, bits;
	tableType tt;
	void dump();
	void write(ofstream &);
	void read(istream &);
	void clear();
	bool smoothflag;
	void sqrt_tf();

	Table& operator+=(Table& rhs);
	Table* operator*(Table& rhs);
	Table* operator-(Table& rhs);
	Table* operator*(number amount);
	Table& operator/=(number amount);
	Table* powt(number exponent);

	void padd(Table& rhs, number exponent);

	number loggaussprob(Table& mean, Table& variance);
};

#endif /* TABLE_H_ */
