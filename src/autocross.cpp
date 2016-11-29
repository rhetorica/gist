/*
 	Gist autocross.cpp

 		neural network training, memory management, and file I/O

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

#include "autocross.h"
#include "main.h"

#include <string>
#include <map>
#include <vector>

#include <libgen.h>

#include <cstring>

using namespace std;

bool autocross_generate = false;
bool autocross_use = false;

bool** readlabels = NULL; // [dc][cc] target scores for each read
number** autocross_MCW = NULL; // method-class weights [MAX_VOTE][cc]
number** autocross_best_MCW = NULL; // used during training to keep track of current best MCWs [MAX_VOTE][cc]
number** autocross_MC_error = NULL; // method-class accuracies [MAX_VOTE][cc]
number*** autocross_scoretable_MCW = NULL; // [dc][MAX_VOTE][cc]

int permissable_taxrank = GENUS; // include how many taxonomic ranks as correct label? (-1 = strain, 0 = species, 1 = genus, etc.)

number* autocross_rl_true = NULL; // [cc] how many reads are expected per class
number* autocross_rl_seen = NULL; // [cc] how many reads have been found per class

bool* use_rcs = NULL; // [dc] if true, read is in reverse complement orientation

// extern std::map<int,string> taxnames;
extern vector<string> *class_filenames;

// extern string RANK_NAMES[9];

string METHOD_NAMES[MAX_VOTE] = {
	"nt.bayes", // 0
	"aa.bayes", // 1
	"nt_rc.bayes", // 2
	"nt.1nn", // 3
	"aa.1nn", // 4
	"nt_rc.1nn", // 5
	"nt.codelta", // 6
	"aa.codelta", // 7
	"nt_rc.codelta", // 8
	"nt.gmm", // 9
	"aa.gmm", // 10
	"nt_rc.gmm", // 11
	"nt.sparse", // 12
	"aa.sparse", // 13
	"nt_rc.sparse", // 14
	"bwa", // 15
	"16s", // 16
	"prior", // 17
};

extern bool* has_protein; // whether or not a peptide could be extrapolated for a given read

extern number **scoretable[MAX_VOTE];

void autocross_prep() {
	cerr << "[head] Preparing Autocross fields." << endl;

	autocross_MCW = new number*[MAX_VOTE];

	foreach(k, MAX_VOTE) {
		autocross_MCW[k] = new number[cc];
		foreach(i, cc) {
			autocross_MCW[k][i] = 0;
		}
	}

	if(autocross_generate) {
		autocross_MC_error = new number*[MAX_VOTE];
		autocross_best_MCW = new number*[MAX_VOTE];

		foreach(k, MAX_VOTE) {
			autocross_best_MCW[k] = new number[cc];
			autocross_MC_error[k] = new number[cc];
			foreach(i, cc) {
				autocross_MC_error[k][i] = 0;
				autocross_best_MCW[k][i] = 0;
			}
		}

		autocross_scoretable_MCW = new number**[dc];

		foreach(j, dc) {
			autocross_scoretable_MCW[j] = new number*[MAX_VOTE];
			foreach(k, MAX_VOTE) {
				autocross_scoretable_MCW[j][k] = new number[cc];
				foreach(i, cc) {
					autocross_scoretable_MCW[j][k][i] = 0;
				}
			}
		}
	}

	cerr << "[head] Autocross prep complete." << endl;
}

void autocross_cleanup() {
	foreach(k, MAX_VOTE) {
		delete[] autocross_MCW[k];
	}
	delete[] autocross_MCW;

	if(autocross_generate) {
		if(use_rcs != NULL) {
			delete[] use_rcs;
			use_rcs = NULL;
		}

		foreach(j, dc) {
			foreach(k, MAX_VOTE) {
				delete[] autocross_scoretable_MCW[j][k];
			}
			delete[] autocross_scoretable_MCW[j];
		}

		foreach(k, MAX_VOTE) {
			delete[] autocross_MC_error[k];
			delete[] autocross_best_MCW[k];
		}
		delete[] autocross_MC_error;

		delete[] autocross_scoretable_MCW;
		delete[] autocross_best_MCW;

		delete[] autocross_rl_seen;
		delete[] autocross_rl_true;
	}
}

// Returns a boolean indicating if the specified classifier is a reverse complement classifier.
bool method_is_rc(VOTETYPE x) {
	if(x == NT_NB_RC || x == NT_NN_RC || x == NT_MM_RC || x == NT_CD_RC || x == NT_SB_RC) return true;
	return false;
}

// Returns a boolean indicating if the specified classifier is a non-reverse-complement nucleotide classifier.
bool method_has_rc_partner(VOTETYPE x) {
	if(x == NT_NB || x == NT_NN || x == NT_MM || x == NT_CD || x == NT_SB) return true;
	return false;
}

void autocross_load_weights(string Wfn) {
	std::ifstream f;

	cerr << "[head] Reading autocross weights from file: " << Wfn << endl;

	f.open(Wfn.c_str(), ios::in);

	if(f.good()) {
		while(f.good()) {
			size_t taxon_id;
			f >> taxon_id;

			bool this_done = false;

			foreach(i, cc) {

				if(strain_taxon[i] == taxon_id) {
					string classifier_type;
					f >> classifier_type;

					forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
						if(METHOD_NAMES[k] == classifier_type) {
							if(autocross_MCW[k][i] == 0) {
								f >> autocross_MCW[k][i];
								this_done = true;
							}

							// cout << "Loaded autocross weight: " << METHOD_NAMES[k] << " weight for strain " << taxnames[taxon_id] << " = " << autocross_MCW[k][i] << endl;
							break;
						}
					}

					if(this_done) break; // multiple classes may have the same taxon ID, so we'll reload them in order
				}
			}
		}

		f.close();
	} else {
		cerr << "[head] ERROR: Could not read autocross weights from " << Wfn << endl;
	}

}

void autocross_save_weights(string Wfn) {
	std::ofstream f;

	f.open(Wfn.c_str(), ios::out);

	if(f.good()) {
		foreach(i, cc)
			forclassifiers(k)
				if(!method_is_rc((VOTETYPE)k))
					f << strain_taxon[i] << "\t" << METHOD_NAMES[k] << "\t" << autocross_MCW[k][i] << endl;

		f.close();

		cerr << "[head] Done writing autocross weights." << endl;
	} else {
		cerr << "[head] ERROR: Could not write autocross weights to " << Wfn << endl;
	}
}

void autocross_load_labels_from_fasta() {
	// assume genome filenames up to first period start the gene name up to first period
	// (e.g. NCBI accession numbers)

	readlabels = new bool*[dc];

	foreach(j, dc) {
		readlabels[j] = new bool[cc];
		foreach(i, cc) readlabels[j][i] = false;
	}

	autocross_rl_seen = new number[cc];
	autocross_rl_true = new number[cc];

	foreach(i, cc) {
		autocross_rl_seen[i] = 0;
		autocross_rl_true[i] = 0;
	}

	foreach(j, dc) {
		size_t ji = j * 2;
		char* endpos = strchr(dv[ji], '|');
		if(endpos == 0) endpos = strchr(dv[ji], ':');
		if(endpos == 0) endpos = strchr(dv[ji], '.');
		size_t endset = (size_t)endpos - (size_t)dv[ji];
		char* dm = new char[endset+1];
		strncpy(dm, dv[ji], endset);
		dm[endset] = '\0';
		size_t period = (size_t)strchr(dm, '.');
		if(period != 0) dm[(size_t)period - (size_t)dm] = '\0';
		size_t colon = (size_t)strchr(dm, ':');
		if(colon != 0) dm[(size_t)colon - (size_t)dm] = '\0';
		size_t pipe = (size_t)strchr(dm, '|');
		if(pipe != 0) dm[(size_t)pipe - (size_t)dm] = '\0';
		// strcpy()
		// string dm = dv[ji]->substr(0, endset);
		//cout << "[allff] Read " << dv[j * 2] << endl
		//	 << "        belongs to " << dm << endl;


		bool matched = false;

		foreach(i, cc) {
			size_t lastslash = class_filenames->at(i).rfind("/");

			string bn = class_filenames->at(i).substr(lastslash + 1);

			string bns = bn.substr(0, bn.find('.'));

			// cout << "BNS: " << bns << endl;

			if(strcmp(bns.c_str(), dm) == 0) {
//				cout << "Matched " << class_filenames->at(i) << endl; // << " (" << taxnames[strain_taxon[i]] << ")" << endl;
				matched = true;

				readlabels[j][i] = true;
				autocross_rl_true[i]++;

				// add friends and family:
				if(permissable_taxrank != -1) {
					foreach(k, cc) {
						if(taxa[k][permissable_taxrank] == taxa[i][permissable_taxrank]) {
							// taxa[k][permissable_taxrank] = 1;
							readlabels[j][i] = true;

							/*cout << "Added " << class_filenames->at(k) << " (" << taxnames[strain_taxon[k]] << ")" << endl;
							cout << " (because taxa(" << k << ", " << permissable_taxrank << ") equals taxa(" << i << ", " << permissable_taxrank << ")): "
									<< "taxon " << taxa[k][permissable_taxrank] << " == taxon " <<  taxa[i][permissable_taxrank] << endl;
							cout << " (same " << RANK_NAMES[permissable_taxrank] << ")" << endl;*/
						}
					}
				}

				break;
			}
		}

		if(!matched) {
			cerr << "[allff] Couldn't find class for read " << j << endl
			     << "        read name: " << dv[j*2] << endl
			     << "        apparent label: " << dm << endl;

		}

		delete[] dm;
	}

	cerr << "[head] Autocross labels loaded." << endl;
}

// returns the re-weighting (from console switches) for the given method.
number simpleweight(VOTETYPE a) {
	switch(a) {
	case NT_NB:
	case NT_NB_RC:
		return bn_weight;
	case NT_MM:
	case NT_MM_RC:
		return mn_weight;
	case NT_CD:
	case NT_CD_RC:
		return cn_weight;
	case NT_NN:
	case NT_NN_RC:
		return nn_weight;
	case NT_SB:
	case NT_SB_RC:
		return sn_weight;
	case AA_NB:
		return bp_weight;
	case AA_MM:
		return mp_weight;
	case AA_CD:
		return cp_weight;
	case AA_NN:
		return np_weight;
	case AA_SB:
		return sp_weight;
	case NT_BWA:
		return bwa_weight;
	case PRIOR_16S:
		return p16s_weight;
	case PRIOR_AP:
		return ap_weight;
	default:
		cerr << "Unexpected weight inquiry during training: " << a << endl;
		exit(98);
	}
}

size_t classification_hits;

number autocross_score(bool mode) {
	number overall_error = 0;

	classification_hits = 0;

	if(use_rcs == NULL) {
		use_rcs = new bool[dc];

		foreach(j, dc) {
			if(n_prefix_length > 0) {
				number fw_score = 0;
				number rc_score = 0;

				foreach(i, cc) {
					if(USE_METHOD[NT_NB]) {
						fw_score += scoretable[NT_NB][j][i] * bn_weight;
						rc_score += scoretable[NT_NB_RC][j][i] * bn_weight;
					}
					if(USE_METHOD[NT_MM]) {
						fw_score += scoretable[NT_MM][j][i] * mn_weight;
						rc_score += scoretable[NT_MM_RC][j][i] * mn_weight;
					}
					if(USE_METHOD[NT_CD]) {
						fw_score += scoretable[NT_CD][j][i] * cn_weight;
						rc_score += scoretable[NT_CD_RC][j][i] * cn_weight;
					}
					if(USE_METHOD[NT_NN]) {
						fw_score += scoretable[NT_NN][j][i] * nn_weight;
						rc_score += scoretable[NT_NN_RC][j][i] * nn_weight;
					}
					if(USE_METHOD[NT_SB]) {
						fw_score += scoretable[NT_SB][j][i] * sn_weight;
						rc_score += scoretable[NT_SB_RC][j][i] * sn_weight;
					}
				}

				if(rc_score > fw_score) {
					use_rcs[j] = true;
				} else {
					use_rcs[j] = false;
				}
			}
		}
	}

	foreach(i, cc)
		autocross_rl_seen[i] = 0;

	for(size_t j = 0; j < (size_t)dc; j++) {
		bool use_rc = use_rcs[j];
		number votes[cc];
		number ovotes[cc];

		// foreach(k, MAX_VOTE) votes[k] = 0;

		foreach(i, cc) {
			votes[i] = 0;

			if(USE_METHOD[NT_BWA]) votes[i] += scoretable[NT_BWA][j][i] * bwa_weight;

			if(n_prefix_length > 0) {
				if(use_rc) {
					if(USE_METHOD[NT_NB_RC]) votes[i] += scoretable[NT_NB_RC][j][i] * bn_weight * autocross_MCW[NT_NB_RC][i];
					if(USE_METHOD[NT_MM_RC]) votes[i] += scoretable[NT_MM_RC][j][i] * mn_weight * autocross_MCW[NT_MM_RC][i];
					if(USE_METHOD[NT_CD_RC]) votes[i] += scoretable[NT_CD_RC][j][i] * cn_weight * autocross_MCW[NT_CD_RC][i];
					if(USE_METHOD[NT_NN_RC]) votes[i] += scoretable[NT_NN_RC][j][i] * nn_weight * autocross_MCW[NT_NN_RC][i];
					if(USE_METHOD[NT_SB_RC]) votes[i] += scoretable[NT_SB_RC][j][i] * sn_weight * autocross_MCW[NT_SB_RC][i];
				} else {
					if(USE_METHOD[NT_NB]) votes[i] += scoretable[NT_NB][j][i] * bn_weight * autocross_MCW[NT_NB][i];
					if(USE_METHOD[NT_MM]) votes[i] += scoretable[NT_MM][j][i] * mn_weight * autocross_MCW[NT_MM][i];
					if(USE_METHOD[NT_CD]) votes[i] += scoretable[NT_CD][j][i] * cn_weight * autocross_MCW[NT_CD][i];
					if(USE_METHOD[NT_NN]) votes[i] += scoretable[NT_NN][j][i] * nn_weight * autocross_MCW[NT_NN][i];
					if(USE_METHOD[NT_SB]) votes[i] += scoretable[NT_SB][j][i] * sn_weight * autocross_MCW[NT_SB][i];
				}
			}

			if(p_prefix_length > 0 && has_protein[j]) {
				if(USE_METHOD[AA_NB]) votes[i] += scoretable[AA_NB][j][i] * bp_weight * autocross_MCW[AA_NB][i];
				if(USE_METHOD[AA_MM]) votes[i] += scoretable[AA_MM][j][i] * mp_weight * autocross_MCW[AA_MM][i];
				if(USE_METHOD[AA_CD]) votes[i] += scoretable[AA_CD][j][i] * cp_weight * autocross_MCW[AA_CD][i];
				if(USE_METHOD[AA_NN]) votes[i] += scoretable[AA_NN][j][i] * np_weight * autocross_MCW[AA_NN][i];
				if(USE_METHOD[AA_SB]) votes[i] += scoretable[AA_SB][j][i] * sp_weight * autocross_MCW[AA_SB][i];
			}

			if(USE_METHOD[PRIOR_16S]) votes[i] += class_16s[i] * p16s_weight * autocross_MCW[PRIOR_16S][i];
			if(USE_METHOD[PRIOR_AP]) votes[i] += class_ap[i] * ap_weight * autocross_MCW[PRIOR_AP][i];
		}

		// votes are now final log scores.

		// let's get some final scores!

		foreach(i, cc) {
			if(isnan(votes[i])) votes[i] = 0;
		}

		number vmax = logsumexp(votes, cc);

		foreach(i, cc) {
			ovotes[i] = exp(votes[i] - vmax); // * autocross_CW[i];
		}

		// pick new best vote:

		int maxi = -1;
		number maxvote = 0;

		foreach(i, cc) {
			if(ovotes[i] > maxvote || maxi == -1) {
				maxi = i;
				maxvote = ovotes[i];
			}

			if(mode == false)
				overall_error += pow(readlabels[j][i] - ovotes[i], 2.0);

			forclassifiers(ki) {
				VOTETYPE k = (VOTETYPE)ki;
				if(use_rcs[j]) {
					if(method_is_rc(k))
						k = (VOTETYPE)(k - 2);
					else if(method_has_rc_partner(k))
						continue;
				} else {
					if(method_is_rc(k))
						continue;
				}

				autocross_MC_error[k][i] += ((float)(readlabels[j][i]) - ovotes[i]) * (isnan(scoretable[k][j][i]) ? 0 : scoretable[k][j][i]) * autocross_MCW[k][i] * simpleweight((VOTETYPE)k);
			}
		}

		if(readlabels[j][maxi]) classification_hits++;
		autocross_rl_seen[maxi]++;
	}

	forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) foreach(i, cc) autocross_MC_error[k][i] /= (number)dc;

	// RMSE:
	return sqrt(overall_error / cc / dc);
}

number learning_rate = 0.02;

#define AC_MAX_ITER 1000
#define MAX_FAILURES 2

extern number** overall_scoretable[MAX_RUNLEVEL];

void autocross_train() {

	int hits[MAX_VOTE];
	forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) hits[k] = 0;

	foreach(j, dc) {
		forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
			number maxv = 0;
			int maxi = -1;
			foreach(i, cc) {
				if(scoretable[k][j][i] > maxv || maxi == -1) {
					maxv = scoretable[k][j][i];
					maxi = i;
				}
				if((method_has_rc_partner((VOTETYPE)k) && scoretable[k+2][j][i] > maxv) || maxi == -1) {
					maxv = scoretable[k+2][j][i];
					maxi = i;
				}
			}
			if(readlabels[j][maxi]) hits[k]++;
		}
	}

	cout << "Method performance:" << endl;
	forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
		cout << "Method " << k << " scored " << hits[k] << "/" << dc << endl;
	}

	int overall_hits = 0;

	foreach(j, dc) {
		number maxv = 0;
		int maxi = -1;
		foreach(i, cc) {
			if(overall_scoretable[0][j][i] > maxv || maxi == -1) {
				maxv = overall_scoretable[0][j][i];
				maxi = i;
			}
		}
		if(readlabels[j][maxi]) overall_hits++;
	}

	cout << "Standard runlevel scored " << overall_hits << "/" << dc << endl;

	srand(time(NULL));

	foreach(k, MAX_VOTE) foreach(i, cc) autocross_MCW[k][i] = 0;
	forclassifiers(k) {
		if(method_is_rc((VOTETYPE)k))
			foreach(i, cc) autocross_MCW[k][i] = 0;
		else {
			foreach(i, cc) autocross_MCW[k][i] = ((number)(rand()) / (number)RAND_MAX) / 4.0 + 0.375;
		}
	}

	number last_op;

	cerr << endl;

	size_t fails = 0;

	size_t best_class_hits = 0;

	foreach(n, AC_MAX_ITER) {
		cerr << "[autocross] Starting iteration n = " << n;
		number op = autocross_score(false);

		cerr << "\x1b[2K\r";
		// cerr << "\r[" << n << "] " << "op: " << op;
		cerr << "[" << n << "] " << "op: " << op << "; " << (int)(((double)classification_hits / (double)dc) * 100) << "% classified correctly";
		cerr.flush();
		cout << "performance = " << classification_hits << "/" << dc << endl;

		if(classification_hits >= best_class_hits) {
			forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
				foreach(i, cc) {
					autocross_best_MCW[k][i] = autocross_MCW[k][i];
				}
			}

			best_class_hits = classification_hits;
		}

		if(n != 0 && op >= last_op) {
			cout << "Not optimizing in iteration " << n << "." << endl;
			fails++;

			learning_rate = -learning_rate / 2;

			cout << "last op: " << last_op << endl;

			if(fails > MAX_FAILURES) {
				cout << "Final performance = " << classification_hits << "/" << dc << endl;
				break;
			}
		} else {
			if(fails > 0) fails--;
		}

		number MCWsum = 0;
		foreach(i, cc) {
//			number rl_error = (autocross_rl_true[i] - autocross_rl_seen[i]) / dc; // ranging from -1 to +1
			forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
				autocross_MCW[k][i] += autocross_MC_error[k][i] * learning_rate; // * (1 + rl_error);
				if(autocross_MCW[k][i] < 0) autocross_MCW[k][i] = 0.0001;
				if(autocross_MCW[k][i] > 1) autocross_MCW[k][i] = 0.9999;
				MCWsum += autocross_MCW[k][i];
			}
		}

		forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
			foreach(i, cc) {
				autocross_MCW[k][i] /= MCWsum;
			}
		}

		last_op = op;
	}

#ifdef SCORES_BEFORE_LOWS
	if(classification_hits < best_class_hits) {
		forclassifiers(k) if(!method_is_rc((VOTETYPE)k)) {
			foreach(i, cc) {
				autocross_MCW[k][i] = autocross_best_MCW[k][i];
			}
		}

		cerr << "[autocross] rejected performance of " << classification_hits << "/" << dc << "." << endl
			 << "            (already achieved " << best_class_hits << "/" << dc << "; using that instead)" << endl;
	}
#endif

	cerr << endl;

	cout << "[autocross] Stopped training." << endl;
}
