/*
 * main.cpp
 *
 *  Created on: 2013-08-08
 *      Author: rhetorica
 */

#include "main.h"

using namespace std;

size_t sz1, sz2, sz3;

bool exists(string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

bool is_file(string& name) {
	struct stat buffer;
	stat(name.c_str(), &buffer);
	return S_ISREG(buffer.st_mode);
//	return (buffer.st_mode & S_IFREG) > 0;
}

bool is_dir(string& name) {
	struct stat buffer;
	stat(name.c_str(), &buffer);
	return S_ISDIR(buffer.st_mode);
//	return (buffer.st_mode & S_IFREG) > 0;
}

void crawl_tree(string& dirname, vector<string>* filenames) {
	DIR *dir;
	struct dirent *ent = NULL;
	dir = opendir(dirname.c_str());
//	struct stat buf, bufg;

	unsigned long ccnt = 0;

	if(dir) {
		while((ent = readdir(dir))) {
			// if(ent->d_type == 4) continue;

			if(ent->d_name == string("..") || ent->d_name == string(".")) continue;

			string fname = string(dirname).append("/").append(string(ent->d_name));

			string ext4 = fname.substr(fname.length() - 5, 5);
			string ext3 = fname.substr(fname.length() - 4, 4);
			string ext2 = fname.substr(fname.length() - 3, 3);

			if(ext2 == string(".gz")) continue;
			if(ext2 == string(".sa")) continue;
			if(ext3 == string(".txt")) continue;
			if(ext3 == string(".amb")) continue;
			if(ext3 == string(".ann")) continue;
			if(ext3 == string(".bwt")) continue;
			if(ext3 == string(".pac")) continue;
			if(ext3 == string(".sam")) continue;
			if(ext3 == string(".bam")) continue;
			if(ext3 == string(".sai")) continue;
			if(ext3 == string(".lin")) continue;
			if(ext4 == string(".gist")) continue;

			// cerr << "-> " << fname << " (" << (int)(ent->d_type) << ")" << endl;

			if(exists(fname) && is_file(fname)) {
				cerr << "    " << ent->d_name << endl;
				filenames->push_back(fname);
				ccnt++;
			} else if(is_dir(fname)) {
				crawl_tree(fname, filenames);
			}

		}
	} else {
		cerr << "Could not open " << dirname << "." << endl;
		exit(85);
	}

	closedir(dir);
}

string gi_from_fasta(string filename) {
	ifstream e;
	e.open(filename.c_str(), ios::in);
	string l0;
	e >> l0;
	// cout << l0 << endl;
	int st = l0.find("|");
	int en = l0.find("|", st + 1);
	if(st == -1) return "";
	if(en == -1) return "";
	l0 = l0.substr(st + 1, en - st - 1);

	e.close();
	return l0;
}

FILE *fp1;

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

FILE *fp2;

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
			cerr << "Stopped at " << q << " while searching for " << num << ". Found: " << here << "." << endl;
			cerr << "Lower bound: " << lo << "." << endl;
			cerr << "Upper bound: " << hi << "." << endl;
			cerr << "Findings: " << logb << endl;
			exit(663);

			exit(99);
		}
	}
}

FILE *fp3;

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

void build_ancestry(ostream& linfile, string tax_id) {
	string ptaxid = tax_id;

	string cellular_organism = "131567";
	string root = "1";

	while(ptaxid != root && ptaxid != cellular_organism) {
		string prow = get_parent_node(ptaxid.c_str(), atoi(ptaxid.c_str()), ptaxid.length(), 0, sz2);
		int of = prow.find("\t|\t");
		int of2 = prow.find("\t|\t", of + 3);

		string ptaxid2 = prow.substr(of + 3, of2 - of - 3);
		string prank = prow.substr(of2 + 3, prow.find("\t|\t", of2 + 3) - of2 - 3);

		string pnamerow = get_name_row(ptaxid.c_str(), atoi(ptaxid.c_str()), ptaxid.length(), 0, sz3);
		int nof = pnamerow.find("\t|\t");
		int nof2 = pnamerow.find("\t|\t", nof + 3);
		if(nof == -1) { cerr << "!!" << pnamerow << "!!" << endl; exit(55); }
		if(nof2 == -1) { cerr << "!!2" << endl; exit(55); }
		string pname = pnamerow.substr(nof + 3, nof2 - nof - 3);

		linfile << prank << ": " << ptaxid << " (" << pname << ")" << endl;

		if(ptaxid2 == ptaxid) break;
		ptaxid = ptaxid2;
	}
}

int main(int argc, char** argv) {
	vector<string> filenames;

	string rootdir;
	if(argc > 1) {
		rootdir = string(argv[1]);

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

	if(argc == 3 && rootdir == string("-solo")) {
		string tax_id = string(argv[2]);
		build_ancestry(cout, tax_id);
		return 0;
	} else if(argc != 2) {
		cerr << "delin 1.1" << endl;
		cerr << "2013 rhetorica@cs.toronto.edu" << endl;
		cerr << endl;
		cerr << "Syntax: " << string(argv[0]) << " <directory tree path>" << endl;
		cerr << "     or " << string(argv[0]) << " -solo <taxid>" << endl;
		cerr << endl;
		cerr << "delin creates .ffn.lin files for all .ffn files it encounters. files must" << endl;
		cerr << "have standard NCBI sequence headers." << endl;
		cerr << endl;
		cerr << "all necessary support files (gi_taxid_nucl.dmp, nodes.dmp, names.dmp) must" << endl;
		cerr << "be located in the working directory." << endl;
		cerr << endl;
		cerr << "to get a specific hit, you can also use -solo to make delin print to stdout." << endl;
		return 2;
	}

	crawl_tree(rootdir, &filenames);

	foreach(i, filenames.size()) {
		ofstream linfile;

		linfile.open(string(filenames[i]).append(".lin").c_str(), ios::out);

		string gi = gi_from_fasta(filenames[i]);
		string tax_id = get_tax_id(gi.c_str(), atoi(gi.c_str()), gi.length(), 0, sz1);
		tax_id = tax_id.substr(0, tax_id.length() - 1);
		if(tax_id.length() == 0) {
			cerr << "No taxon ID: " << filenames[i] << endl;
			continue;
		}

		cout << filenames[i] << " (GI #" << gi << ") is taxon ID #" << tax_id << "." << endl;

		build_ancestry(linfile, tax_id);

		linfile.close();
	}
}
