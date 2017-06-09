/*
 	Gist Table.cpp

 		core functions for read tables

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

#include "Table.h"
#include <sstream>
#include <cstring>

void Table::clear() {
	if(smoothflag) this->length = pow(cols, this->prefix_length + 1);

	if(tt == NUCLEOTIDE) {
		for(size_t i = 0; i < rows; i++) {
			for(size_t n = 0; n < cols; n++)
				if(smoothflag) {
					this->mydata[i][n] = 1 / this->length; // not laplace smoothing
					if(isnan(this->mydata[i][n])) {
						cerr << "Houston, we have a problem." << endl;
					}
				} else {
					this->mydata[i][n] = 0;
					if(isnan(this->mydata[i][n])) {
						cerr << "Canaveral, we have a problem." << endl;
					}
				}
		}
	} else {
		for(size_t i = 0; i < rows; i++) {
			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			for(size_t n = 0; n < cols; n++)
				if(smoothflag) {
					this->mydata[i][n] = 1 / this->length; // not laplace smoothing
					if(isnan(this->mydata[i][n])) {
						cerr << "Assert failed: generated infinite value after laplace smoothing during table clearing." << endl;
					}
				} else {
					this->mydata[i][n] = 0;
					if(isnan(this->mydata[i][n])) {
						cerr << "Assert failed: generated infinite value during table clearing. Your computer is broken." << endl;
					}
				}
		}
	}
}

number Table::delta(Table* d) {
	number sum = 0;
	if(tt == NUCLEOTIDE) {
		for(size_t i = 0; i < rows; i++) {
			for(size_t n = 0; n < cols; n++)
#ifdef NN_D2_STAR
				sum += (number)(mydata[i][n] * d->mydata[i][n]);
#else
				sum += pow((number)(mydata[i][n]) / (number)(this->length) - (number)(d->mydata[i][n]) / (number)(d->length), 2);
#endif
		}
	} else {
		for(size_t i = 0; i < rows; i++) {
			//cout << i;

			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			for(size_t n = 0; n < cols; n++)
#ifdef NN_D2_STAR
				sum += (number)(mydata[i][n] * d->mydata[i][n]);
#else
				sum += pow((number)(mydata[i][n]) / (number)(this->length) - (number)(d->mydata[i][n]) / (number)(d->length), 2);
#endif
		}
	}

#ifdef NN_D2_STAR
	sum /= sqrt(this->length * d->length);
#endif

	return sum;
}

string* nfrombits(unsigned long int e, size_t pl) {
	char output[pl + 1];
	output[pl] = 0;
	for(size_t i = 0; i < pl; i++) {
		char outty = e & 3;
		e = e >> 2;
		output[pl - i - 1] = nucleotides[(size_t)outty];
	}

	return new string(output);
}

string* pfrombits(unsigned long int e, size_t pl) {
	char output[pl + 1];
	output[pl] = 0;
	for(size_t i = 0; i < pl; i++) {
		char outty = e & 31;
		e = e >> 5;
		if(outty > 19) return NULL;
		output[pl - i - 1] = peptides[(size_t)outty];
	}

	return new string(output);
}

void Table::sqrt_tf() {
	if(no_sqrt) return;
	number** md = mydata;

	if(tt == NUCLEOTIDE) {
		for(size_t i = 0; i < rows; i++) {
			for(size_t n = 0; n < cols; n++) {
				md[i][n] = sqrt(md[i][n]);
			}
		}
	} else {
		for(size_t i = 0; i < rows; i++) {

			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			string* prefixmsg = pfrombits(i, this->prefix_length);
			if(prefixmsg != NULL) {
				delete prefixmsg; prefixmsg = NULL;

				for(size_t n = 0; n < cols; n++) {
					md[i][n] = sqrt(md[i][n]);
				}
			}
		}
	}
}

// Constructor for reloading table from file
Table::Table(string& serialized, int prefix_size, tableType newtt) {
	smoothflag = false;
	this->tt = newtt;

	if(tt == NUCLEOTIDE) {
		base = 4;
		cols = 4;
		bits = 2;
	} else {
		cols = 20;
		base = 32;
		bits = 5;
	}
	rows = (int)pow(base, prefix_size);
	this->prefix_length = prefix_size;

	this->length = -1; // handled below, in read()

	this->mydata = new number*[rows];

	for(size_t i = 0; i < rows; i++) {
		this->mydata[i] = new number[cols];
	}

	stringstream x(serialized);
	read(x);
}

// Read table from encoded string
void Table::read(istream& input) {
	number** md = mydata;
	number o;

	input >> o;
	this->length = o;

	if(tt == NUCLEOTIDE) {
		for(size_t i = 0; i < rows; i++) {
			for(size_t n = 0; n < cols; n++) {
				input >> o;
				//cerr << i << "/" << rows << "," << n << "/" << cols << ":" << o << endl;
				md[i][n] = o;
			}
		}
	} else {
		for(size_t i = 0; i < rows; i++) {

			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			string* prefixmsg = pfrombits(i, this->prefix_length);
			if(prefixmsg != NULL) {
				delete prefixmsg; prefixmsg = NULL;

				for(size_t n = 0; n < cols; n++) {
					input >> o;
					//cerr << i << "," << n << ":" << o << endl;
					md[i][n] = o;
				}
			}
		}
	}
}

void Table::write(ofstream& f) {
	number** md = mydata;

	f << " " << this->length;

	if(tt == NUCLEOTIDE) {
		for(size_t i = 0; i < rows; i++) {
			for(size_t n = 0; n < cols; n++) {
				f << " " << md[i][n];
			}
		}
	} else {
		for(size_t i = 0; i < rows; i++) {

			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			string* prefixmsg = pfrombits(i, this->prefix_length);
			if(prefixmsg != NULL) {
				delete prefixmsg; prefixmsg = NULL;

				for(size_t n = 0; n < cols; n++)
					f << " " << md[i][n];
			}
		}
	}
}

void Table::dump() {
	number** md = mydata;

	cout << "Length: " << length << "<br>";

	const int chroma_inflate = 2;

	if(tt == NUCLEOTIDE) {
		//cout << "prefix\tT       C       A       G<br>";

		for(size_t n = 0; n < cols; n++) {
			for(size_t i = 0; i < rows; i++) {
				string* prefixmsg = nfrombits(i, this->prefix_length);
				//cout << *prefixmsg << "\t";
				delete prefixmsg; prefixmsg = NULL;

				cout << "<div title='" << (number)(md[i][n]) << "' ";
				number cval = sqrt((number)(md[i][n])) * 255 * chroma_inflate;
				printf("style='background: rgb(%7.0f, %7.0f, %7.0f); width: 16px; height: 16px; display: inline-block'></div>", cval, cval, cval);
			}
			cout << "<br>";
		}
	} else {
		cout << "prefix\tF       L       S       Y       C       W       P       H       Q       R       I       M       T       N       K       V       A       D       E       G" << endl;

		for(size_t i = 0; i < rows; i++) {

			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}

			string* prefixmsg = pfrombits(i, this->prefix_length);
			if(prefixmsg != NULL) {
				cout << *prefixmsg << "\t";
				delete prefixmsg; prefixmsg = NULL;

				for(size_t n = 0; n < cols; n++)
					printf("<div style='background-color: rgb(%7f, %7f, %7f); width: 16px; height: 16px; display: inline-block'></div>", (number)(md[i][n]), (number)(md[i][n]), (number)(md[i][n]));

				cout << "<br>";
			}
		}
	}
}

void Table::padd(Table& rhs, number exponent) {
	if(this->tt != rhs.tt) {
		cerr << "Incompatible table types in padd() call." << endl;
		throw new string("Incompatible table types in padd() call.");
	}

	if(this->prefix_length > rhs.prefix_length) {
		cerr << "[padd] Can't extend: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Can't extend.");
	} else if(this->prefix_length < rhs.prefix_length) {
		// TODO
		cerr << "[padd] Not yet implemented: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Not yet implemented.");
	} else {
		for(size_t i = 0; i < rows; i++) {
			for(size_t j = 0; j < cols; j++) this->mydata[i][j] += pow(rhs.mydata[i][j], exponent);
		}
	}
}

Table* Table::powt(number exponent) {
	string emptystring = string("");
	Table* result = new Table(&emptystring, this->prefix_length, this->tt);

	for(size_t i = 0; i < rows; i++) {
		for(size_t j = 0; j < cols; j++) result->mydata[i][j] = pow(this->mydata[i][j], exponent);
	}

	return result;
}

Table* Table::operator*(Table& rhs) {
	if(this->tt != rhs.tt) {
		cerr << "[op*] Incompatible table types in * call." << endl;
		throw new string("Incompatible table types in * call.");
	}

	string emptystring = string("");

	Table* result = new Table(&emptystring, rhs.prefix_length, rhs.tt);

	if(this->prefix_length > rhs.prefix_length) {
		cerr << "[op*] Can't extend: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Can't extend.");
	} else if(this->prefix_length < rhs.prefix_length) {
		// TODO
		cerr << "Not yet implemented: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Not yet implemented.");
	} else {

		number** md = mydata;
		number** yd = rhs.mydata;

		if(this->length < 0 || rhs.length < 0) {
			cerr << "Multiplying tables of negative length." << endl;
			exit(3);
		}

		result->length = this->length * rhs.length;

		for(size_t i = 0; i < rows; i++) {
			for(size_t j = 0; j < cols; j++) result->mydata[i][j] = md[i][j] * yd[i][j];
		}
	}

	return result;
}

Table* Table::operator*(number amount) {

	string emptystring = string("");

	Table* result = new Table(&emptystring, this->prefix_length, this->tt);

	number** md = mydata;

	for(size_t i = 0; i < rows; i++) {
		for(size_t j = 0; j < cols; j++) {
			if(isnan(result->mydata[i][j])) {
				cerr << "OP*D: Found out that " << i << ", " << j << " is a NaN." << endl;
			}
			result->mydata[i][j] = md[i][j] * amount;
			if(isnan(result->mydata[i][j])) {
				cerr << "OP*D: Bad! Multiplication of " << i << ", " << j << " with " << amount << " caused NaN." << endl;
				exit(4);
			}
		}
	}
	result->length = this->length;

	return result;
}

Table* Table::operator-(Table& rhs) {
	if(this->tt != rhs.tt) {
		cerr << "[op-] Incompatible table types in - call." << endl;
		throw new string("Incompatible table types in - call.");
	}

	string emptystring = string("");

	Table* result = new Table(&emptystring, rhs.prefix_length, rhs.tt);

	if(this->prefix_length > rhs.prefix_length) {
		cerr << "[op-] Can't extend: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Can't extend.");
	} else if(this->prefix_length < rhs.prefix_length) {
		// TODO
		cerr << "[op-] Not yet implemented: table has prefix length of " << this->prefix_length << ", but partner is " << rhs.prefix_length << "." << endl;
		throw new string("Not yet implemented.");
	} else {

		number** md = mydata;
		number** yd = rhs.mydata;

		for(size_t i = 0; i < rows; i++) {
			for(size_t j = 0; j < cols; j++) result->mydata[i][j] = md[i][j] - yd[i][j];
		}
	}

	return result;
}

Table& Table::operator+=(Table& rhs) {
	if(this->tt != rhs.tt) {
		cerr << "[op+=] Incompatible table types in sum call." << endl;
		throw new string("Incompatible table types in sum call.");
	}

	if(this->prefix_length > rhs.prefix_length) {
		cerr << "[op+=] Can't extend: got " << this->prefix_length << " and " << rhs.prefix_length << " for prefix lengths." << endl;
		throw new string("Can't extend.");
	} else if(this->prefix_length < rhs.prefix_length) {
		// TODO
		cerr << "[op+=] Not yet implemented: got " << this->prefix_length << " and " << rhs.prefix_length << " for prefix lengths." << endl;
		throw new string("Not yet implemented.");
	} else {
		number** md = mydata;
		number** yd = rhs.mydata;

		for(size_t i = 0; i < rows; i++) {
			for(size_t j = 0; j < cols; j++) md[i][j] += yd[i][j];
		}
		if(rhs.length < 0) {
			cerr << "Adding tables with less than zero length?" << endl;
			cerr << "rhs length = " << rhs.length << endl;
			cerr << "rhs:" << endl;
			rhs.dump();
			exit(5);
		}
		this->length += rhs.length;
	}
	return *this;
}

Table& Table::operator/=(number amount) {

	if(amount == 0) {
		string oe = string("Attempt to /= ").append(this->tt == PROTEIN ? "protein table" : "nucleotide table").append(" by zero.\n");
		cerr << oe << endl;
		throw oe;
	}

	for(size_t i = 0; i < rows; i++) {
		for(size_t j = 0; j < cols; j++) {
			if(isnan(this->mydata[i][j] / amount)) {
				cerr << "IsNaN: table element " << i << ", " << j << " in " << (this->tt == PROTEIN ? "protein table" : "nucleotide table") << " after a /= by " << amount << " (was " << this->mydata[i][j] << ")" << endl;
				exit(7);
			}
			this->mydata[i][j] /= amount;
		}
	}

	this->length /= amount;

	return *this;
}

number Table::loggaussprob(Table& mean, Table& variance) {
	if(this->tt != mean.tt || this->tt != variance.tt || variance.tt != mean.tt) {
		cerr << "====================================================" << endl;
		cerr << "Incompatible table types in NBC call." << endl;
		cerr << "====================================================" << endl;
		cerr << "mean table " << &mean << " is of type " << mean.tt << endl;
		cerr << "var table " << &variance << " is of type " << variance.tt << endl;
		cerr << "data table " << this << " is of type " << this->tt << endl;
		cerr << "mean table has " << mean.rows << " rows, variance table has " << variance.rows << " rows, data table has " << this->rows << " rows." << endl;
		cerr << "====================================================" << endl;
		throw new string("Incompatible table types in NBC call.");
	}

	if(this->prefix_length > mean.prefix_length) {
		cerr << "[lgp] Prefix length mismatch: got " << this->prefix_length << " for data and " << mean.prefix_length << " for mean." << endl;
		throw new string("Can't extend.");
	} else if(this->prefix_length < mean.prefix_length) {
		// TODO
		cerr << "[lgp] Prefix length mismatch: got " << this->prefix_length << " for data and " << mean.prefix_length << " for mean." << endl;
		throw new string("Not yet implemented.");
	}

	number ans = 0;

	const number pi = atan(1)*4;

	int fvb = 0;

	for(size_t i = 0; i < rows; i++) {

		if(this->tt != NUCLEOTIDE) {
			size_t j = i, k = 0;
			while(j > base) {
				j /= base;
				k++;
			}
			if(j > cols - 1) {
				i = pow(base, k+1);
				if(i >= rows) break;
			}
		}

		for(size_t j = 0; j < cols; j++) {
			if(variance.mydata[i][j] <= 0) {
				variance.mydata[i][j] = 0.00000001;
				fvb++;
			}
			number thisparta = 2 * pi * variance.mydata[i][j];

			number thispartb = pow(thisparta, 0.5);

			if(thisparta <= 0 || thispartb <= 0) {
				cerr << "Unstable variance value  " << variance.mydata[i][j] << " at " << i << ", " << j << " in class with " << mean.length << " data samples." << endl;
				throw new string("Unstable variance value!");
			}

			//cout << variance.mydata[i][j] << ", ";
			number thispart = 1 / thispartb;

			number thatpart = 0 - ( pow(this->mydata[i][j] - mean.mydata[i][j], 2.0) / (2 * variance.mydata[i][j]) );
			//cout << thatpart << endl;

			ans += log(thispart) + thatpart;
		}
	}

	if(fvb > 0) {
		//cerr << "    Fixed " << fvb << " elements in variance by tiny-value substitution." << endl;
		//cerr << "    Result caused: " << ans << endl;
	}

	// cerr << "Calculated log gaussian scores " << ans << " for class, normal " << exp(ans) << "." << endl;

	if(isnan(ans)) {
		cerr << "Bad Gaussian score generated." << endl;
		exit(8);
	}

	return ans;
}

Table::Table(char* source, int prefix_size, tableType tt) {
	if(debugdump) cerr << "my address is " << this << "." << endl;

	this->tt = tt;

	this->prefix_length = prefix_size;
	if(this->prefix_length == 0) {
		cerr << "Table with prefix length zero is illegal." << endl;
		throw(string("I think you really should use grep for that."));
	}

	if(debugdump) cerr << "prefix size: " << prefix_size << endl;

	if(tt == NUCLEOTIDE) {
		base = 4;
		cols = 4;
		bits = 2;
	} else {
		cols = 20;
		base = 32;
		bits = 5;
	}
	rows = (int)pow(base, prefix_size);

	size_t srclen = strlen(source);

	this->length = srclen - prefix_length; // number of k-mers

	smoothflag = false;

	if(strcmp(source, "smooth") == 0) {
		smoothflag = true;
		this->length = pow(cols, prefix_size + 1);
	}

	if(this->length < 0) this->length = 0;

	this->mydata = new number*[rows];

	for(size_t i = 0; i < rows; i++) {
		this->mydata[i] = new number[cols];
		for(size_t j = 0; j < cols; j++) {
			if(smoothflag) {
				this->mydata[i][j] = 1 / this->length; // not laplace smoothing
			} else {
				this->mydata[i][j] = 0;
			}
		}
	}

	if(this->length < 1 || smoothflag) {
		return;
	}

	long unsigned int mask = pow(base, prefix_size) - 1;
	long unsigned int prefix = 0;
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
			mydata[prefix][(size_t)focus]++;
		}
	}
}

Table::Table(string* source, int prefix_size, tableType tt) {
	if(debugdump) cerr << "creating table at " << this << "." << endl;

	this->tt = tt;

	this->prefix_length = prefix_size;
	if(this->prefix_length == 0) {
		cerr << "Table with prefix length zero is illegal." << endl;
		throw(string("You should use grep for that."));
	}

	if(debugdump) cerr << "prefix size: " << prefix_size << endl;

	if(tt == NUCLEOTIDE) {
		base = 4;
		cols = 4;
		bits = 2;
	} else {
		cols = 20;
		base = 32;
		bits = 5;
	}
	rows = (int)pow(base, prefix_size);

	this->length = source->length() - prefix_length;

	smoothflag = false;

	if(*source == string("smooth")) {
		smoothflag = true;
		this->length = pow(cols, prefix_size + 1);
	}

	if(this->length < 0) this->length = 0;

	this->mydata = new number*[rows];

	for(size_t i = 0; i < rows; i++) {
		this->mydata[i] = new number[cols];
		for(size_t j = 0; j < cols; j++) {
			if(smoothflag) {
				this->mydata[i][j] = 1 / this->length; // not laplace smoothing
			} else {
				this->mydata[i][j] = 0;
			}
		}
	}

	if(this->length < 1 || smoothflag) {
		return;
	}

	long unsigned int mask = pow(base, prefix_size) - 1;
	long unsigned int prefix = 0;
	int skips = 0;
	for(size_t i = 0; i < source->length() - 1; i++) {
		char newfix, focus;
		if(tt == NUCLEOTIDE) {
			newfix = twobits((*source)[i], i);
			focus = twobits((*source)[i+1], i+1);
		} else {
			newfix = fivebits((*source)[i]);
			focus = fivebits((*source)[i+1]);
			if(focus == -1) skips = prefix_size + 1;
		}
		prefix = ((prefix << bits) + newfix) & mask;
		if(skips != 0) {
			skips = skips - 1;
		} else if(i >= (size_t)prefix_size - 1) {
			mydata[prefix][(size_t)focus]++;
		}
	}
}

Table::~Table() {
	for(size_t i = 0; i < rows; i++) {
		delete[] mydata[i];
		mydata[i] = NULL;
	}
	delete[] mydata; mydata = NULL;
}

