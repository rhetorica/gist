/*
 * fasta.cpp
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#include "fasta.h"

#include <fstream>

#include "progress.h"

extern bool f_mode;

void scrub(string& b) {
	for(int i = b.length() - 1; i >= 0; i--) {
		switch(b[i]) {
		case 0:
		case ' ':
		case '\r':
		case '\t':
			b.erase(i, 1);
			break;
		default:
			break;
		}
	}
}

int readFasta(string filename, char**& records) {
	ifstream f;
		size_t recordcount = 0;
		string buffer = string("");
		f.open(filename.c_str(), ios::in | ios::binary);

		if(f.is_open()) {

			char linebufpre[2048];

			// first line:

			if(f.good()) {
				f.getline(linebufpre, 2048, '\n');

				string linebuf = linebufpre;

				if(linebuf[0] == '>') {
					buffer.append(linebuf).append("\n");

					recordcount++;
				} else {

					return -2; // Not a FASTA file

					scrub(linebuf);
					buffer.append(linebuf);
				}
			}

			// subsequent lines:

			while(f.good()) {
				f.getline(linebufpre, 2048, '\n');

				string linebuf = linebufpre;

				if(linebuf[0] == '>') {
					buffer.append("\n").append(linebuf).append("\n");

					recordcount++;
				} else {
					scrub(linebuf);
					buffer.append(linebuf);
				}
			}
		} else {
			return -1;
		}

		cerr << "[load] Loading " << recordcount << " sequences from FASTA file " << filename << "." << endl;

		prepclock(recordcount, "records", "loading data");

		buffer.append(1, '\n');

		size_t bufferlength = buffer.length();
		const char* buff2 = buffer.c_str();

		records = new char*[recordcount * 2];

		int t = -1, start = 0, end = 0;

		for(size_t i = 0; i < bufferlength; i++) {
			if(buff2[i] == '\n') {
				end = i-1;

				t++;

				/*
				if(t > 0) {
					if(buff2[start] == '>') {
						start++;
						const char* endi = strchr(buff2, ' ');
						if(endi != NULL) end = endi - 1 - buff2;
					}
					// string tempa = buffer.substr(start, end - start + 1);
					// records[t - 1] = strdup(tempa.c_str());
				}
	*/

				if(buff2[start] == '>') start++;

				records[t] = new char[end - start + 2];
				strncpy(records[t], (char*)(buff2 + start), end - start + 1);
				records[t][end - start + 1] = '\0';
	/*
				cout << endl;
				cout << t << ": " << string(records[t]) << "." << endl;
				cout << "length = " << (end - start + 1) << endl;
				cout << "allocated storage = " << (end - start + 2) << endl;
				cout << "start = " << start << endl;
				cout << "end = " << end << endl;
				cout << "buff2 @ " << (long)buff2 << endl;
				cout << "i = " << i << endl;
				cout << "t = " << t << endl;
	*/
				start = i+1;

				prog_status = t / 2;
				runClock();
			}
		}

		f.close();

		cerr << endl << "[load] Loading from " << filename << " complete." << endl;

		// exit(55);

		return recordcount;
}

int load_class_id_from_fasta_label(const char* header) {
	// assume genome filenames up to first period start the gene name up to first period
	// (e.g. NCBI accession numbers)

	const char* endpos = strchr(header, '|');
	if(endpos == 0) endpos = strchr(header, ':');
	if(endpos == 0) endpos = strchr(header, '.');
	if(endpos == 0) return -1;
	size_t endset = (size_t)endpos - (size_t)header;
	char* dm = new char[endset+1];
	strncpy(dm, header, endset);
	dm[endset] = '\0';
	size_t period = (size_t)strchr(dm, '.');
	if(period != 0) dm[(size_t)period - (size_t)dm] = '\0';
	size_t colon = (size_t)strchr(dm, ':');
	if(colon != 0) dm[(size_t)colon - (size_t)dm] = '\0';
	size_t pipe = (size_t)strchr(dm, '|');
	if(pipe != 0) dm[(size_t)pipe - (size_t)dm] = '\0';
	// strcpy()
	// string dm = dv[ji]->substr(0, endset);
	// cout << "[allff] Read belongs to " << dm << endl;

	bool matched = false;

	foreach(i, class_filenames.size()) {
		size_t lastslash = class_filenames.at(i).rfind("/");

		string bn = class_filenames.at(i).substr(lastslash + 1);

		string bns = bn.substr(0, bn.find('.'));

		// cout << "BNS: " << bns << endl;

		if(strcmp(bns.c_str(), dm) == 0) {
			// cout << "Matched " << class_filenames->at(i) << " (" << taxnames[strain_taxon[i]] << ")" << endl;
			matched = true;

			return i;

			break;
		}
	}

	if(!matched && !f_mode) {
		cerr << "[allff] Couldn't find class for read!" << endl
			 << "        read name: " << header << endl
			 << "        apparent label: " << dm << endl;

	}

	delete[] dm;

	return -1;
}

int load_class_id_from_filename(const char* header) {
	// assume genome filenames from the full string
	// (e.g. NCBI accession numbers)

	size_t endset = strlen(header);
	char* dm = new char[endset+1];
	strncpy(dm, header, endset);
	dm[endset] = '\0';

	bool matched = false;

	foreach(i, class_filenames.size()) {
		size_t lastslash = class_filenames.at(i).rfind("/");

		string bn = class_filenames.at(i).substr(lastslash + 1);

		string bns = bn.substr(0, bn.find('.'));

		// cout << "BNS: " << bns << endl;

		if(strcmp(bns.c_str(), dm) == 0) {
			// cout << "Matched " << class_filenames->at(i) << " (" << taxnames[strain_taxon[i]] << ")" << endl;
			matched = true;

			return i;

			break;
		}
	}

	if(!matched) {
		cerr << "[allff] Couldn't find class for read!" << endl
			 << "        read name: " << header << endl
			 << "        apparent label: " << dm << endl;

	}

	delete[] dm;

	return -1;
}

