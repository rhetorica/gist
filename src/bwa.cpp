/*
 	Gist bwa.cpp

 		wrapper for BWA calls

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

using namespace std;

#include "bwa.h"
#include "main.h"
#include "Class.h"
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>

#include <sys/stat.h>
#include <map>

#include "classify.h"

string bname(string& fname) {
	string cmd(fname);
	if(cmd.substr(cmd.length() - 5, 5) == string(".gist")) {
		cmd = cmd.substr(0, cmd.length() - 5);
	}

	return cmd;
}

void parseSAMs(string& filename) {

	string clock_unit = "files";
	string clock_tag = "parsing aligner results";
	time_t clock_time = time(NULL);
	size_t prog_status = 0;

	std::map<string, int> readindex;
	foreach(j, dc) {
		string readname = string(dv[j * 2]);
		string short_readname;
		size_t off = readname.find(" ");
		if(off == string::npos)
			short_readname = readname;
		else
			short_readname = readname.substr(0, off);

		readindex[short_readname] = j;
	}

	foreach(i, cc) {
		std::ostringstream ts, tc;
		ts << i;

		ifstream f;

		string sfn = string(filename).append(ts.str()).append(".sam");

		// cerr << "=== BWA results file " << sfn << " ===" << endl;

		f.open(sfn.c_str(), ios::in);

		if(f.is_open()) {
			while(f.good()) {
				string ibuff;
				getline(f, ibuff);
				if(ibuff[0] == '@') continue;

				string sname = ibuff.substr(0, ibuff.find('\t'));
				if(sname.length() == 0) continue;

				bool fe = (bool)(readindex.count(sname));

				if(fe) {
					size_t j = readindex[sname];

					int asqpos = ibuff.find("\tAS:i:") + 6;
					int asqlastpos = ibuff.find('\t', asqpos + 1);
					int asq = atoi(ibuff.substr(asqpos, asqlastpos - asqpos).c_str());

					if(asq > scoretable[NT_BWA][j][i]) scoretable[NT_BWA][j][i] = asq;

					prog_status++;
					runClock(prog_status, cc, clock_time, clock_unit, clock_tag);
				} else {
					cerr << "Unknown gene: \"" << sname << "\"" << endl;
					cerr << ".sam file does not match source FASTA." << endl;
					exit(80);
				}
			}
		} else {
			cerr << "Couldn't open SAM file " << sfn << endl;
			exit(81);
		}

		f.close();
		remove(sfn.c_str());
	}
	cerr << endl;
}

void indexBWA(string& c) {
	string fn = bname(c);
	string fn2 = fn;

	if(!exists(fn2.append(".pac"))) {
		string cmd = string("bwa index ").append(fn);
		int ioo = system(cmd.c_str());
		if(ioo != 0) {
			if(ioo == -1) {
				cerr << "Couldn't generate BWA index for class " << c << endl;
				cerr << "Command: " << cmd << endl;
				cerr << "Return code -1: the system could not create the new process." << endl;
			} else {
				cerr << "Couldn't generate BWA index for class " << c << endl;
				cerr << "BWA return code " << ioo << endl;
			}
			exit(82);
		}
	}
}

void runBWA(vector<string>* class_names, size_t threadcount, string& filename) {
	// return; // temporary!

	string clock_unit = "classes";
	string clock_tag = "running aligner";
	time_t clock_time = time(NULL);
	size_t prog_status = 0;

	foreach(i, cc) {
		string c = class_names->at(i);

		string fn = bname(c);

		indexBWA(c);

		std::ostringstream ts, tc;
		ts << i;

		string tname = string(filename).append(ts.str());

		if(threadcount == 0) threadcount = 1;

		tc << threadcount;
		string tcount = tc.str();

		string cmd = string("bwa mem -v 0 -t ").append(tcount).append(" ").append(fn).append(" ").append(filename).append(" > ").append(tname).append(".sam");
		// cerr << cmd << endl;
		int ioo = system(cmd.c_str());
		if(ioo != 0) {
			cerr << "Couldn't run BWA-MEM alignment for class " << c << endl;
			cerr << "BWA return code " << ioo << endl;
			exit(83);
		}
		prog_status++;
		runClock(prog_status, cc, clock_time, clock_unit, clock_tag);
	}

	cerr << endl;
}

