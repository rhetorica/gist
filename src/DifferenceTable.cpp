/*
 	Gist DifferenceTable.cpp

 		an n x n table (tightly-packed) used for codelta information

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

#include "DifferenceTable.h"
#include <sstream>

// NOTE: this cov stuff isn't used anywhere and probably wouldn't work anyway
// since the matrix is now allocated as upper diagonal only

/*
void DifferenceTable::addCovElement(Table& e, Table& mu, number weight) {
	size_t xcoord = 0;
	size_t ycoord = 0;

	foreach(i1, rows) {
		if(tt != NUCLEOTIDE) adjustRow(i1);
		if(i1 >= rows) break;
		foreach(j1, cols) {

			ycoord = Coord(i1, j1);

			foreach(i2, rows) {
				if(tt != NUCLEOTIDE) adjustRow(i2);
				if(i2 >= rows) break;
				foreach(j2, cols) {

					xcoord = Coord(i2, j2);

					// if(xcoord <= ycoord) continue;

					this->mydata[xcoord][ycoord] += (e.mydata[xcoord] - mu.mydata[xcoord]) * (e.mydata[ycoord] - mu.mydata[ycoord]) * weight;
				}
			}
		}
	}
}

void DifferenceTable::normalizeCov(number* masses) {
	size_t xcoord = 0;
	size_t ycoord = 0;

	number wsum = 0, wsumsquared = 0, wsquaredsum = 0;

	foreach(wi, Coord(rows - 1, cols - 1) ) {
		wsum += masses[wi];
	}

	wsumsquared = wsum * wsum;

	foreach(wi, Coord(rows - 1, cols - 1) ) {
		wsquaredsum += masses[wi] * masses[wi];
	}

	if(wsum != 1) {
		std::cerr << "normalizeCov: weighted sum should be 1." << std::endl;
	}

	number mass = wsum / (wsumsquared - wsquaredsum);

	foreach(i1, rows) {
		if(tt != NUCLEOTIDE) adjustRow(i1);
		if(i1 >= rows) break;
		foreach(j1, cols) {

			ycoord = Coord(i1, j1);

			foreach(i2, rows) {
				if(tt != NUCLEOTIDE) adjustRow(i2);
				if(i2 >= rows) break;
				foreach(j2, cols) {

					xcoord = Coord(i2, j2);

//					if(xcoord <= ycoord) continue;

					this->mydata[xcoord][ycoord] *= mass;
				}
			}
		}
	}
}
*/

void DifferenceTable::init() {
	if(tt == NUCLEOTIDE) {
		base = 4;
		cols = 4;
		bits = 2;
	} else {
		cols = 20;
		base = 32;
		bits = 5;
	}

	rows = (int)pow(base, prefix_length);
	size = rows * base;

#ifndef DOG_EAR
	halfsize = (int)pow(cols, prefix_length + 1);
#else
	halfsize = (int)pow(cols, prefix_length + 1) / 2;
#endif

	this->mydata = new number*[size];

	foreach(i, (size_t)size) {
		mydata[i] = new number[halfsize];

		foreach(j, (size_t)halfsize) {
			mydata[i][j] = 0;
		}
	}
}

DifferenceTable::DifferenceTable(int prefix_len, tableType newTableType) : prefix_length(prefix_len), tt(newTableType) {
	rows = 0, base = 0, cols = 0, bits = 0, size = 0, mydata = NULL, halfsize = 0;
	this->init();
}

DifferenceTable::~DifferenceTable() {
	foreach(i, (size_t)size) {
		delete[] mydata[i];
		mydata[i] = NULL;
	}

	delete[] mydata;

//	cout << "Destroying difference table." << endl;
}

void DifferenceTable::write(ofstream& f) {
	size_t xcoord = 0;
	size_t ycoord = 0;

	foreach(i1, rows) {
		if(tt != NUCLEOTIDE) adjustRow(i1);
		if(i1 >= rows) break;
		foreach(j1, cols) {

			ycoord = Coord(i1, j1);
			if(ycoord >= halfsize) goto stoppit_dtw;

			foreach(i2, rows) {
				if(tt != NUCLEOTIDE) adjustRow(i2);
				if(i2 >= rows) break;
				foreach(j2, cols) {

					xcoord = Coord(i2, j2);

					f << " " << this->mydata[xcoord][ycoord];
				}
			}
		}
	}

stoppit_dtw:
	return;

}

DifferenceTable::DifferenceTable(string& serialized, int prefix_len, tableType newtt) {
	stringstream x(serialized);
	tt = newtt;
	prefix_length = prefix_len;
	init();

	read(x);
}

void DifferenceTable::read(istream& f) {
	size_t xcoord = 0;
	size_t ycoord = 0;

	number o;

	foreach(i1, rows / 2) {
		if(tt != NUCLEOTIDE) adjustRow(i1);

		foreach(j1, cols) {

			ycoord = Coord(i1, j1);
			if(ycoord >= halfsize) goto stoppit_dtr;

			foreach(i2, rows) {
				if(tt != NUCLEOTIDE) adjustRow(i2);
				if(i2 >= rows) break;
				foreach(j2, cols) {

					xcoord = Coord(i2, j2);

					f >> o;
					this->mydata[xcoord][ycoord] = o;
				}
			}
		}
	}

stoppit_dtr:
	return;
}
