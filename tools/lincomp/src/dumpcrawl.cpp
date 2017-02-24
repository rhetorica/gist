/*
 * dumpcrawl.cpp
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#include "dumpcrawl.h"

FILE *fp1, *fp2, *fp3;

size_t sz1, sz2, sz3;

#define foreach(i, max) for(size_t i = 0; i < max; i++)

string gi_from_fasta(string l0) {
	int st = l0.find("|");
	int en = l0.find("|", st + 1);
/*
	if(st == -1 || en == -1) {
		int st = l0.find(":");
		int en = l0.find(":", st + 1);
	}

	if(st == -1 || en == -1) {
		int st = l0.find(".");
		int en = l0.find(".", st + 1);
	}
*/
	if(st == -1) return "";
	if(en == -1) return "";
	l0 = l0.substr(st + 1, en - st - 1);

	return l0;
}

void load_db() {
	fp1 = fopen("gi_taxid_nucl.dmp", "r");
	fseek(fp1, 0L, SEEK_END);
	sz1 = ftell(fp1);

	fp2 = fopen("nodes.dmp", "r");
	fseek(fp2, 0L, SEEK_END);
	sz2 = ftell(fp2);

	fp3 = fopen("names.dmp", "r");
	fseek(fp3, 0L, SEEK_END);
	sz3 = ftell(fp3);
}

string get_tax_id(const char* gi, size_t num, size_t lens, size_t lo, size_t hi) {
	size_t q = (lo + hi) / 2;
	// cout << q << "?" << endl;
	fseek(fp1, q, SEEK_SET);
	while(fgetc(fp1) != '\n' && q > 0) {
		q--;
		fseek(fp1, q, SEEK_SET);
	}

	// cout << q << " vs. " << lo << "!" << endl;

	char here[100];

	fgets(here, 99, fp1);

	if(strncmp(gi, here, lens) == 0) {
		if(here[lens] == '\t') {
			// exact match
//			cerr << "Win!" << endl;
			//exit(66);
			return string(here).substr(lens + 1);
		} else {
			// too high
			return get_tax_id(gi, num, lens, lo, q);
		}
	} else {
		// cout << "  " << num << " != " << here << endl;
		size_t sl = strlen(here);
		foreach(j, 99) {
			if(here[j] == '\t') here[j] = '\0';
		}
		size_t thenny = atoi(here);

		if(q == lo) {
			if(hi > lo + sl) {
				return get_tax_id(gi, num, lens, lo + sl, hi);
			} else {
//				cerr << "Lose!" << endl;
				return "";

				/*fseek(fp1, q, SEEK_SET);
				char logb[200];
				fread(logb, sizeof(char), 199, fp1);
				logb[199] = '\0';
				cerr << "Stopped at " << q << " while searching for " << num << ". Found: " << here << "." << endl;
				cerr << "Upper bound: " << hi << "." << endl;
				cerr << "Findings: " << logb << endl;
				exit(663);*/
			}
		}

		if(thenny > num) {
			return get_tax_id(gi, num, lens, lo, q);
		} else if(thenny < num) {
			return get_tax_id(gi, num, lens, q, hi);
		} else {
			cerr << "WHAT HAVE YOU DONE?! Puking on " << gi << endl;
			exit(99);
		}
	}
}

string get_parent_node(const char* taxid, size_t num, size_t lens, size_t lo, size_t hi) {
	size_t q = (lo + hi) / 2;
	// cout << q << "?" << endl;
	fseek(fp2, q, SEEK_SET);
	while(fgetc(fp2) != '\n' && q > 0) {
		q--;
		fseek(fp2, q, SEEK_SET);
	}

	// cout << q << " vs. " << lo << "!" << endl;

	char here[200];

	fgets(here, 199, fp2);

	if(strncmp(taxid, here, lens) == 0) {
		if(here[lens] == '\t') {
			// exact match
//			cerr << "Win!" << endl;
			//exit(66);
			return string(here);
		} else {
			// too high
			return get_parent_node(taxid, num, lens, lo, q);
		}
	} else {
		size_t sl = strlen(here);
		foreach(j, 199) {
			if(here[j] == '\t') here[j] = '\0';
		}
		size_t thenny = atoi(here);

		if(q == lo) {
			if(hi > lo + sl) {
				return get_parent_node(taxid, num, lens, lo + sl, hi);
			} else {
				return "";
			}
		}

		if(thenny > num) {
			return get_parent_node(taxid, num, lens, lo, q);
		} else if(thenny < num) {
			return get_parent_node(taxid, num, lens, q, hi);
		} else {
			cerr << "GORDON'S ALIVE?!" << endl;

			fseek(fp2, q, SEEK_SET);
			char logb[200];
			fread(logb, sizeof(char), 199, fp2);
			logb[199] = '\0';
			cerr << "Stopped at " << q << " while searching for " << num << " (" << taxid << "). Found: " << here << "." << endl;
			cerr << "Lower bound: " << lo << "." << endl;
			cerr << "Upper bound: " << hi << "." << endl;
			cerr << "Context: " << logb << endl;
			exit(663);

			exit(99);
		}
	}
}

string get_name_row_up(const char* taxid, size_t lens, size_t q) {
	q -= 1;
	fseek(fp3, q, SEEK_SET);
	while(fgetc(fp3) != '\n' && q > 0) {
		q--;
		fseek(fp3, q, SEEK_SET);
	}

	char here[200];

	fgets(here, 199, fp3);

	if((strncmp(taxid, here, lens) == 0) && (here[lens] == '\t')) {
		// right taxid

		if(strstr(here, "\t|\tscientific name\t|") == 0) {
			return get_name_row_up(taxid, lens, q);
		} else {
			return string(here);
		}
	} else {
		return string(""); // failed
	}
}

string get_name_row_down(const char* taxid, size_t lens, size_t q, size_t hi) {
	q += 1;
	fseek(fp3, q, SEEK_SET);
	while(fgetc(fp3) != '\n' && q < hi) {
		q++;
		fseek(fp3, q, SEEK_SET);
	}

	char here[200];

	fgets(here, 199, fp3);

	if((strncmp(taxid, here, lens) == 0) && (here[lens] == '\t')) {
		// right taxid

		if(strstr(here, "\t|\tscientific name\t|") == 0) {
			return get_name_row_down(taxid, lens, q, hi);
		} else {
			return string(here);
		}
	} else {
		return string(""); // failed
	}
}

string get_name_row(const char* taxid, size_t num, size_t lens, size_t lo, size_t hi) {
	size_t q = (lo + hi) / 2;
	// cout << q << "?" << endl;
	fseek(fp3, q, SEEK_SET);
	while(fgetc(fp3) != '\n' && q > 0) {
		q--;
		fseek(fp3, q, SEEK_SET);
	}

	// cout << q << " vs. " << lo << "!" << endl;

	char here[200];

	fgets(here, 199, fp3);

	if(strncmp(taxid, here, lens) == 0) {
		if(here[lens] == '\t') {
			// exact match ... but is it the RIGHT exact match?!

			if(strstr(here, "\t|\tscientific name\t|") == 0) {
				// keep digging!
				string guess = get_name_row_up(taxid, lens, q);
				if(guess == "") {
					guess = get_name_row_down(taxid, lens, q, hi);
				}

				if(guess == "") {
					cerr << "Couldn't look up name for taxid " << taxid << "." << endl;
				}

//				cerr << "((found " << guess << "))" << endl;

				return guess;
			} else {
				// we win!
//				cerr << "((got " << here << "))" << endl;
				return string(here);
			}
		} else {
			// too high
			return get_name_row(taxid, num, lens, lo, q);
		}
	} else {
		size_t sl = strlen(here);
		foreach(j, 199) {
			if(here[j] == '\t') here[j] = '\0';
		}
		size_t thenny = atoi(here);

		if(q == lo) {
			if(hi > lo + sl) {
				return get_name_row(taxid, num, lens, lo + sl, hi);
			} else {
				return "";
			}
		}

		if(thenny > num) {
			return get_name_row(taxid, num, lens, lo, q);
		} else if(thenny < num) {
			return get_name_row(taxid, num, lens, q, hi);
		} else {
			cerr << "Robots in disguise?!" << endl;

			fseek(fp3, q, SEEK_SET);
			char logb[200];
			fread(logb, sizeof(char), 199, fp3);
			logb[199] = '\0';
			cerr << "Stopped at " << q << " while searching for " << num << ". Found: " << here << "." << endl;
			cerr << "Lower bound: " << lo << "." << endl;
			cerr << "Upper bound: " << hi << "." << endl;
			cerr << "Findings: " << logb << endl;
			exit(663);

			exit(99);
		}
	}
}
