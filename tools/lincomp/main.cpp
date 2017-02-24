/*
 * main.cpp
 *
 *  Created on: 2013-08-08
 *      Author: rhetorica
 */

#include "main.h"
#include "progress.h"
#include "fasta.h"
#include "dumpcrawl.h"
#include "parselin.h"
#include "dircrawl.h"
#include <stdexcept>

bool s_mode = false;
bool p_mode = false;
bool f_mode = false;

map<int, string> nrank;
map<int, int> parent;
map<int, string> nodename;

enum filter_type {
	NONE,
	TAXON,
	LABEL,
	SEQUENCE,
};

struct filter_rule {
	bool retain;
	filter_type type;
	int number;
	string fragment;
};

void load_ancestry(string tax_id) {
	string ptaxid = tax_id;

	string cellular_organism = "131567";
	string root = "1";

	while(ptaxid != root && ptaxid != cellular_organism) {
		string prow = get_parent_node(ptaxid.c_str(), atoi(ptaxid.c_str()), ptaxid.length(), 0, sz2);
		int of = prow.find("\t|\t");
		int of2 = prow.find("\t|\t", of + 3);

		string ptaxid2 = prow.substr(of + 3, of2 - of - 3);
		string prank = prow.substr(of2 + 3, prow.find("\t|\t", of2 + 3) - of2 - 3);

		string pname;

		if(ptaxid != "") {
			string pnamerow = get_name_row(ptaxid.c_str(), atoi(ptaxid.c_str()), ptaxid.length(), 0, sz3);

			int nof = pnamerow.find("\t|\t");
			int nof2 = pnamerow.find("\t|\t", nof + 3);
			if(nof == -1) { cerr << "!!" << pnamerow << "!!" << endl; exit(55); }
			if(nof2 == -1) { cerr << "!!2" << endl; exit(55); }

			pname = pnamerow.substr(nof + 3, nof2 - nof - 3);
		} else {
			pname = "(nothing)";
		}

		int tid = atoi(ptaxid.c_str());

		parent[tid] = atoi(ptaxid2.c_str());
		nrank[tid] = prank;
		nodename[tid] = pname;

		// linfile << prank << ": " << ptaxid << " (" << pname << ")" << endl;

		if(ptaxid2 == ptaxid) break; // circular?! awash in a sea of tears
		if(nodename.count(atoi(ptaxid2.c_str()))) break; // already known
		ptaxid = ptaxid2;
	}
}

void charstar_to_lower(char* s) {
	for( ; *s; ++s) *s = tolower(*s); // apparently this is idiomatic in C now (ugh)
}

int main(int argc, char** argv) {
	if(argc < 3) {
		cerr << "lincomp 1.42 taxonomic grep" << endl
		     << "2014-2015 rhetorica@cs.toronto.edu" << endl
		     << endl
		     << "Syntax: " << string(argv[0]) << " [prec|sens|filter {<taxon|label|seq> <+|-><id>, ...}] <guess> <reference>" << endl
		     << endl
		     << " prec: report precision (symmetrical between guess and reference files)" << endl
		     << " sens: report sensitivity of classifications in guess against truth in reference" << endl
		     << "     (default: symmetrical between both files; reports precision)" << endl
			 << " filter: send subset of reads to STDOUT" << endl
	         << "   taxon +<id>: remove all reads outside specified taxid" << endl
		     << "   taxon -<id>: remove all reads within specified taxid" << endl
			 << "   label +<id>: remove all reads lacking specified label substring" << endl
			 << "   label -<id>: remove all reads with specified label substring" << endl
			 << "   seq +<id>: remove all reads lacking specified subsequence" << endl
			 << "   seq -<id>: remove all reads including specified subsequence" << endl
	         << "   (taxon takes labels from the guess file and reads from the reference file)" << endl
             << "   (label and seq are case-insensitive and work just on the reference file)" << endl
             << "   (hint: use the same file for guess and reference when extracting from mixed FFNs)" << endl
		     << endl
		     << "lincomp reports the level of taxonomic consensus for two sequences. NCBI FFN" << endl
 		     << "files are recognized, as are two-column CSVs containing gene ID and taxid." << endl
			 << endl
			 << "filtering features require that the reference file be in FASTA format." << endl
		     << endl
		     << "all necessary support files (gi_taxid_nucl.dmp, nodes.dmp, names.dmp) must" << endl
		     << "be located in the working directory." << endl
		     << endl
		     << "output: a two-column CSV containing sequence identifiers and the level of consensus." << endl
		     << endl
		     << "based on delin 1.1, rhetorica@cs.toronto.edu" << endl;
		return 2;
	}

	load_db();

	string cwd(".");

	crawl_tree(cwd, &class_filenames);

	vector<string> filenames;

	vector<filter_rule> filters;

	if(string(argv[1]) == string("prec")) {
		if(argc != 4) {
			cerr << "Inappropriate arguments." << endl << endl
			     << "Precision testing requires two input files and supports no other parameters." << endl;
		}
		p_mode = true;
		filenames.push_back(string(argv[2]));
		filenames.push_back(string(argv[3]));
	} else if(string(argv[1]) == string("sens")) {
		if(argc != 4) {
			cerr << "Inappropriate arguments." << endl << endl
                 << "Sensitivity testing requires two input files and supports no other parameters." << endl;
			return 2;
		}
		s_mode = true;
		filenames.push_back(string(argv[2]));
		filenames.push_back(string(argv[3]));
	} else if(string(argv[1]) == string("filter")) {
		f_mode = true;

		if(argc < 5) {
			cerr << "Insufficient arguments." << endl << endl
				 <<  "At least one filter type, one filter, and two files are required." << endl;
			return 2;
		}

		filters.clear();

		filter_type ft = NONE;
		for(int i = 2; i < argc - 2; i++) {
			if(string(argv[i]) == "taxon") {
				ft = TAXON;
			} else if(string(argv[i]) == "label") {
				ft = LABEL;
			} else if(string(argv[i]) == "seq") {
				ft = SEQUENCE;
			} else {
				filter_rule new_rule;
				new_rule.type = ft;
				if(argv[i][0] == '+') {
					new_rule.retain = true;
				} else if(argv[i][0] == '-') {
					new_rule.retain = false;
				} else {
					cerr << "Bad arguments." << endl << endl
				       	 << "Unknown filter rule term: '" << string(argv[i]) << "'" << endl;
					return 2;
				}

				if(ft == TAXON) {
					new_rule.number = atoi(&(argv[i][1]));
					new_rule.fragment = "";
				} else {
					new_rule.fragment = string(argv[i]).substr(1);
					transform(new_rule.fragment.begin(), new_rule.fragment.end(), new_rule.fragment.begin(), ::tolower);
					new_rule.number = 0;
				}

				filters.push_back(new_rule);
/*				if(new_rule.type == TAXON) {
					fprintf(stderr, "[filter] rule %d: taxon %s %d\n", filters.size(),
						(new_rule.retain ? "retain" : "remove"),
						new_rule.number);
				} else {
					fprintf(stderr, "[filter] rule %d: %s %s %s\n", filters.size(),
						(new_rule.type == LABEL ? "label" : "sequence"),
						(new_rule.retain ? "retain" : "remove"),
						new_rule.fragment.c_str());
				} */
			}
		}
		filenames.push_back(string(argv[argc - 2]));
		filenames.push_back(string(argv[argc - 1]));
	} else {
		cerr << "Incorrect syntax!" << endl << endl
			 << "Run '" << string(argv[0]) << "' with no parameters for more information." << endl;
		return 2;
	}

	char** fasta[2];
	map<string, int> label[2]; // read -> taxid
	map<string, char*> sequence[2]; // read -> sequence
	
	size_t gi_fail_count = 0;
	int gdc = 0;

	size_t reads_kept = 0;
	size_t reads_killed = 0;

	foreach(i, filenames.size()) {
		int dc = readFasta(filenames[i], fasta[i]);

		if(dc > 0) {
			gdc += dc;
			prepclock(dc, "reads", "parse FASTA");
			foreach(j, (size_t)dc) {

				string tax_id = "-1";

				try {
					string gi;
					try {
						gi = gi_from_fasta(fasta[i][j * 2]);
					} catch(const void*& e) {
						if(!f_mode) {
							cerr << "Got stuck on line " << j << " while finding gi." << endl;
							cerr << fasta[i][j * 2];
						}
						goto breakpoint;
					}

					if(gi.size() == 0) {
						gi_fail_count++;
						int cl = load_class_id_from_fasta_label(fasta[i][j * 2]);
						if(cl != -1) tax_id = tax_id_from_lin(class_filenames[cl]);
						// cout << filenames[i] << " - " << fasta[i][j * 2] << "\n\t(cl #" << cl << ") is taxon ID #" << tax_id << "." << endl;
					} else {
						tax_id = get_tax_id(gi.c_str(), atoi(gi.c_str()), gi.length(), 0, sz1);
						if(tax_id == "") {
							if(!f_mode)
								cerr << "No taxon ID was retrieved for GI #" << gi << "; trying .lin lookup..." << endl;

							try {
								int cl = load_class_id_from_fasta_label(fasta[i][j * 2]);
								if(cl != -1) tax_id = tax_id_from_lin(class_filenames[cl]);
							} catch(const void*& e) {
								if(!f_mode) {
									cerr << "Got stuck on line " << j << " while finding tax_id." << endl;
									cerr << fasta[i][j * 2];
								}
								goto breakpoint;
							}
						}
						// cout << filenames[i] << " - " << fasta[i][j * 2] << "\n\t(GI #" << gi << ") is taxon ID #" << tax_id << "." << endl;
					}

					while(tax_id.find("\n") != string::npos) tax_id.replace(tax_id.find("\n"), 1, "");

					if(tax_id == "0" || tax_id == "-1" || tax_id == "") {
						if(!f_mode)
							cerr << "No taxon ID: " << filenames[i] << " -> " << fasta[i][j * 2] << endl;
						goto breakpoint;
					}

					try {
						load_ancestry(tax_id);
					} catch(const void*& e) {
						if(!f_mode) {
							cerr << "Got stuck on line " << j << " while decoding tax_id " << tax_id << "." << endl;
							cerr << fasta[i][j * 2];
						}
						goto breakpoint;
					}
				} catch(const out_of_range& e) {
					if(!f_mode) {
						cerr << "Failed to get content from line " << j << ":" << endl;
						cerr << fasta[i][j * 2] << endl << endl;
					}
					gi_fail_count++;
					goto breakpoint;
				}

				breakpoint:

				label[i][fasta[i][j * 2]] = atoi(tax_id.c_str());
				sequence[i][fasta[i][j * 2]] = fasta[i][j * 2 + 1];

				prog_status++;
				runClock();
			}
		} else {
			cerr << "Loading non-FASTA file " << filenames[i] << endl;
			ifstream fy;
			fy.open(filenames[i].c_str(), ios::in);
			if(fy.is_open()) {
				try {
					while(fy.good()) {
						gdc++;
						string read_name, tax_id, sbuff;
						istringstream line;
						getline(fy, sbuff);

						line.str(sbuff);

						line >> read_name;
						line >> tax_id;

						if(tax_id.c_str()[0] > '9' || tax_id.c_str()[0] < '0') {
							int cl = load_class_id_from_filename(tax_id.c_str());
							if(cl != -1)
								tax_id = tax_id_from_lin(class_filenames[cl]);
							else if(f_mode == false)
								cerr << "Could not determine label described by string: " << tax_id.c_str();
						}

						label[i][read_name] = atoi(tax_id.c_str());
					}
				} catch(ios_base::failure &e) {
					cerr << "Problem reading " << filenames[i] << "!" << endl;
					cerr << e.what();
					exit(6);
				}
			} else {
				cerr << "Couldn't open " << filenames[i] << "!" << endl;
				exit(5);
			}
		}
	}

	cerr << endl;

	cerr << "Match failure count: " << gi_fail_count << "/" << gdc << endl;

	cerr << "Counting common reads..." << endl;

	map<string, int> presences;

	size_t doublehits = 0;

	foreach(i, 2) {
		for(map<string,int>::iterator it = label[i].begin(); it != label[i].end(); ++it) {
			if(presences.count(it->first) == 0) {
				presences[it->first] = 1;
			} else {
				presences[it->first] = 2;
				doublehits++;
			}
		}
	}

	prepclock(doublehits, "reads", "compare");

	for(map<string,int>::iterator it = presences.begin(); it != presences.end(); ++it) {
		if(presences[it->first] == 1) {
			if(it->first != "")
				cerr << "Warning: did not find \"" << it->first << "\" in both files." << endl;
		} else {
			int Left = label[0][it->first];
			// int Right = label[1][it->first];
			// if(Left == Right) goto donehere;

			int Right = 0;

			bool keep = true;

			if(f_mode) {
				// who to keep in the filtering process?

				if(filters[0].retain) keep = false;

				for(int i = 0; i < filters.size(); i++) {
					char* tempstr;
					char* r;
					switch(filters[i].type) {
						case SEQUENCE:
							tempstr = new char[strlen(sequence[1][it->first])];
							tempstr = strcpy(tempstr, sequence[1][it->first]);
							charstar_to_lower(tempstr);
							r = strstr(tempstr, filters[i].fragment.c_str());

							if(r != NULL) keep = filters[i].retain;

							delete[] tempstr;
							break;
						case LABEL:
							tempstr = new char[strlen((it->first).c_str())];
							tempstr = strcpy(tempstr, (it->first).c_str());
							charstar_to_lower(tempstr);
							r = strstr(tempstr, filters[i].fragment.c_str());

							if(r != NULL) keep = filters[i].retain;

							delete[] tempstr;
							break;
						case TAXON:
							Left = label[0][it->first];
							Right = filters[i].number;
							while(Left != 0) {
								if(Left == Right) {
									keep = filters[i].retain;
									break;
								}
								if(parent.count(Left) == 0) break;
								Left = parent[Left];
							}

							break;
					}
				}

			} else if(s_mode) {
				// check for perfect inclusion with just top result
				// does it match anything in the other chain?
				// if so, the program was sensitive enough to return the top hit

				Right = label[1][it->first];
                while(Right != 0) {
                    if(Left == Right) {
                        Right = label[1][it->first];
                        Left = Right;
                        goto donehere;
                    }
                    Right = parent[Right];
                }

				// check for imperfect inclusion
				// figure out the smallest taxon that *does* match, and return it instead
				// the program failed to find the optimal hit

				while(Left != 0) {
					Right = label[1][it->first];
					while(Right != 0) {
						if(Left == Right) goto donehere;
						Right = parent[Right];
						if(parent.count(Right) == 0) break;
					}
					Left = parent[Left];
					if(parent.count(Left) == 0) break;
				}
			} else if(p_mode) {
				while(Left != 0) {
					Right = label[1][it->first];
					if(Left == Right) goto donehere;

					while(Right != 0) {
						//cerr << Left << " vs " << Right << endl;
						if(parent.count(Right) == 0) break;

						//cerr << " Parent of " << Right << " is " << parent[Right] << endl;
						Right = parent[Right];


						if(Left == Right) goto donehere;
					}

					//cerr << " break " << endl;

					if(parent.count(Left) == 0) break;
					Left = parent[Left];
				}
			}

			donehere:

			if(f_mode) {
				if(keep) {
					cout << ">" << it->first << endl << sequence[1][it->first] << endl;
					reads_kept++;
				} else {
					reads_killed++;
				}
			} else if(s_mode || p_mode) {
				cout << it->first << "\t" << Left << "\t" << label[0][it->first] << "\t" << label[1][it->first] << "\t";

				if(Left == Right) {
					cout << nrank[Right] << "\t" << nodename[Right] << endl;
				} else {
					cout << "no consensus for " << nodename[label[0][it->first]] << " and " << nodename[label[1][it->first]] << endl;
				}
			}
			prog_status++;
			runClock();
		}
	}

	if(f_mode) {
		cerr << endl << endl
			 << "Reads retained: " << reads_kept << endl
			 << "Reads removed: " << reads_killed << endl;
	}

	cerr << endl;

	return 0;
}
