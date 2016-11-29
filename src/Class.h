/*
 	Gist Class.h

 		headers for reference genome (Class) management

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

#ifndef CLASS_H_
#define CLASS_H_

#include "main.h"
#include "Table.h"
#include "DifferenceTable.h"
#include "SparseTable.h"

using namespace std;

class Class {
public:
	Class(string);
	void saveClass(string);
	void rebuildClass(string);
	number* score(size_t row);
	number* getMMscores(Table* datum);
	virtual ~Class();
	void info();

	bool* mmnValid;
	bool* mmpValid;

	string myname;
private:
	void closestDistance(Table*, Table*, Table*, number&, number&, number&);
	void makeCodelta(tableType tty);
	void makeMM(tableType tty);
	void loadField(string, string);
	number gradeCodelta(Table* dT);

	// Table** ntransitions;
	// Table** ptransitions;

	Table*** n_tr;
	Table*** p_tr;

	Table *avgtn, *avgtp;
	Table *vartn, *vartp;

	size_t mmnK;
	size_t mmnKi;
	Table** mmnMeans;
	Table** mmnVars;

	size_t mmpK;
	size_t mmpKi;
	Table** mmpMeans;
	Table** mmpVars;

	DifferenceTable* ncodelta;
	DifferenceTable* pcodelta;

	long long int dead_recc_p;

	long long int recc;
	long long int reci_n;
	long long int reci_p;

	SparseTable* nsparse;
	SparseTable* psparse;
};

#endif /* CLASS_H_ */
