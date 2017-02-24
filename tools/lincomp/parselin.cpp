/*
 * parselin.cpp
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#include "parselin.h"

string tax_id_from_lin(string filename) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	// bool first_line = true;

	int last_idnum = -1;
	string first_idnum = "-1";

	if(f.is_open()) {
		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find(": ");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					string value = s.substr(offset + 2);
					string valname = value.substr(value.find(" ") + 1);
					if(valname[0] == '(' && valname[valname.length() - 1] == ')')
						valname = valname.substr(1, valname.length() - 2);
					string sval = value.substr(0, value.find(" "));
					int idnum = atoi(sval.c_str());
					
					if(nrank.count(idnum) == 0) nrank[idnum] = fieldname;
					
					if(last_idnum != -1) {
						if(parent.count(last_idnum) == 0) {
//							cerr << "Adding " << valname << " to graph as id number " << sval << " = " << idnum << endl;
							parent[last_idnum] = idnum;
						}
					} else {
						first_idnum = sval;
					}
					nodename[idnum] = valname;
					last_idnum = idnum;
				} else {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}
		
		return first_idnum;

		f.close();
	} else {
		string* fruit = new string("Couldn't load lineage file " + filename);
		throw fruit;
	}

	return "-1";
}
