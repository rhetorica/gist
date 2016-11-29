/*
 	Gist DifferenceTable.cpp

 		header for codelta difference table

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

#ifndef DIFFERENCETABLE_H_
#define DIFFERENCETABLE_H_

#include "main.h"
#include "Table.h"

class DifferenceTable {
public:
	DifferenceTable(int plen, tableType newtt);
	DifferenceTable(string& serialized, int, tableType newtt);
	virtual ~DifferenceTable();
	void write(ofstream&);
	void read(istream&);
	number** mydata;
	int prefix_length;
	tableType tt;
	int size;
	size_t halfsize;


	// covariance-only stuff:

	// void addCovElement(Table& e, Table& mu, number weight);
	// void normalizeCov(number* masses);


	void clear();
	size_t base, bits, cols, rows;
private:
	void init();
};

#endif /* DIFFERENCETABLE_H_ */
