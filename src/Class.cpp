/*
 	Gist Class.cpp

 		defines an individual reference genome and the math required to sustain it
 		see also Table.cpp and DifferenceTable.cpp

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

#include "Class.h"
#include "main.h"
#include <time.h>
#include <sstream>
#include "bwa.h"
#include <map>
#include <cstring>

unsigned int class_read_length = 100;

number Class::gradeCodelta(Table* dT) {
	DifferenceTable* cod = NULL;

	//Table* mu;

	size_t rows, cols, base;

	// size_t prefix_length;

	if(dT->tt == NUCLEOTIDE) {
		cod = this->ncodelta;
		//mu = this->avgtn;
		// prefix_length = n_prefix_length;
		if(cod == NULL) {
			cerr << "No nucleotide codelta matrix!" << endl;
			return 0;
		}
	} else if(dT->tt == PROTEIN) {
		cod = this->pcodelta;
		//mu = this->avgtp;
		// prefix_length = p_prefix_length;
		if(cod == NULL) {
			cerr << "No protein codelta matrix!" << endl;
			return 0;
		}
	}

	rows = dT->rows;
	cols = dT->cols;
	base = dT->base;

	// size_t xcoord, ycoord;

	size_t halfsize = cod->halfsize;

	// number exfreq = 1 / pow(cols, prefix_length + 1);

	number sum = 0;
	//number deltasum = 0;
	//number expected_deltasum = 0;

	#ifdef DOG_EAR
		size_t rcons = (rows / base) * (cols / 2);
	#endif

	size_t xcoord, ycoord;

	foreach(i1, rows) {
		if(dT->tt != NUCLEOTIDE) adjustRow(i1);
		if(i1 >= rows) break;
		foreach(j1, cols) {

			ycoord = Coord(i1, j1);
			if(ycoord >= halfsize) goto stoppit_gc;

			foreach(i2, rows) {
				if(dT->tt != NUCLEOTIDE) adjustRow(i2);
				if(i2 >= rows) break;
				foreach(j2, cols) {

					xcoord = Coord(i2, j2);

					#ifdef DOG_EAR
						bool alt = (ycoord > xcoord);

						number yel = dT->mydata[i1 + (alt ? rcons : 0) ][j1];
						number xel = dT->mydata[i2 + (alt ? rcons : 0) ][j2];
					#else
						number yel = dT->mydata[i1][j1];
						number xel = dT->mydata[i2][j2];
					#endif

					// number yel = dT->mydata[i1][j1];// - mu->mydata[i1][j1];
					// number xel = dT->mydata[i2][j2];// - mu->mydata[i2][j2];

					number delta = pow(yel - xel, 2) * sgn(yel - xel);

					//deltasum += abs(delta);
					//expected_deltasum += abs(cod->mydata[xcoord][ycoord]);

					number heart = (delta) - (cod->mydata[xcoord][ycoord]);

					number val = abs(heart);

					// cout << "heart: " << heart << " ";

					// cout << "val: " << val << " ";

					/*
					if((cod->mydata[xcoord][ycoord]) < 0) {
						cerr << "negative codelta after normalization!" << endl;
						return 0;
					}
					*/

					sum += val;
				}
			}
		}
	}

stoppit_gc:

/*
	cout << "Done codelta grading. Sum = " << sum;
	cout << ", data delta sum = " << deltasum;
	cout << ", class expected delta sum = " << expected_deltasum << endl;
	exit(55);*/

	return 2 * sum;
}

void Class::makeMM(tableType tty) {
	Table** means;
	Table** vars;
	Table** data;
	bool* validity;

	size_t K;
	size_t prefix_length;

	srand(time(NULL));

	if(tty == PROTEIN) {
		this->mmpMeans = new Table*[mmpK];
		this->mmpVars = new Table*[mmpK];
		this->mmpValid = new bool[mmpK];

		means = this->mmpMeans;
		vars = this->mmpVars;
		data = this->p_tr[mp_prefix_length];
		K = mmpK;
		prefix_length = mp_prefix_length;
		validity = this->mmpValid;
	} else if(tty == NUCLEOTIDE) {
		this->mmnMeans = new Table*[mmnK];
		this->mmnVars = new Table*[mmnK];
		this->mmnValid = new bool[mmnK];

		means = this->mmnMeans;
		vars = this->mmnVars;
		data = this->n_tr[mn_prefix_length];
		K = mmnK;
		prefix_length = mn_prefix_length;
		validity = this->mmnValid;
	} else {
		cerr << "[EM] Nonsense class type? Probable data corruption." << endl;
		exit(70);
	}

	number resps[K][this->recc]; // how much does each latent variable govern a data element? resps(:,i) sums to 1
	number w[K]; // total prior of each latent variable (Sr / i)
	number Sr[K]; // sum of class k's responsibilities; "nk" in original MATLAB version

	// init means and vars:

	string emptystring = string("");
	string smoothstring = string("smooth");
	foreach(k, K) {
		means[k] = new Table(&smoothstring, prefix_length, tty);
		vars[k] = new Table(&emptystring, prefix_length, tty);
	}

	const size_t MAX_ATTEMPTS = 120;

	bool sr0;

	size_t final_t = 0;

	size_t attempto = 0;

	number bestllh = -1;
	Table* bestmeans[K];
	Table* bestvars[K];

	foreach(k, K) {
		bestmeans[k] = new Table(&smoothstring, prefix_length, tty);
		bestvars[k] = new Table(&emptystring, prefix_length, tty);
	}

	size_t steps_so_far = 0;

	foreach(attempts, MAX_ATTEMPTS) {
		number llh = 0, last_llh = 0; // previous loglikelihood value.

		bool failure = false;
		sr0 = false;

		// initialize responsibilities randomly and normalize:
		foreach(i, (size_t)recc) {
			if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
			number localsum = 0;
			foreach(k, K) {
				resps[k][i] = (((number)rand() / (number)RAND_MAX)) * 0.98 + 0.02;
				localsum += resps[k][i];
			}

			foreach(k, K) {
				resps[k][i] /= localsum;
			}
		}

		if(debuggmm) {
			foreach(k, K) {
				foreach(i, (size_t)recc) {
					if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
					cout << "<div title='" << (number)(resps[k][i]) << "' ";
					printf("style='background: rgb(%7.0f, %7.0f, %7.0f); width: 4px; height: 4px; display: inline-block'></div>", resps[k][i] * 255.0, (number)(resps[k][i]) * 255.0, (number)(resps[k][i]) * 255.0);
				}
				cout << "<br>\n";
			}
		}

		// initialize Sr; used as denominator when scaling for means and variances
		foreach(k, K) {
			Sr[k] = 0;
			foreach(i, (size_t)recc) {
				if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
				Sr[k] += resps[k][i];
			}
			if(tty == PROTEIN)
				w[k] = Sr[k] / (recc - dead_recc_p);
			else
				w[k] = Sr[k] / recc;
		}

		size_t t = 0;

		for(; t < MM_ITER; t++) {
			steps_so_far++;

			if(debuggmm) cerr << t << ". llh: " << llh << endl;

			// calculate normal distributions from responsibilities:
			// (maximization step)

			// means:

			foreach(k, K) {
				means[k]->clear();

				foreach(j, (size_t)recc) {
					if(tty == PROTEIN && p_tr[mp_prefix_length][j] == NULL) continue;
					Table* tempo = (*(data[j])) * resps[k][j];
					//if(j == 20) {
						//cout << "<h2>Cluster " << k << "</h2>";
						//cout << "Data: <br>";
						//data[j]->dump();
						//cout << "Scale: " << resps[k][j] << "<br>";
						//cout << "Tempo: <br>";
						//tempo->dump();
					//}

					*(means[k]) += *tempo;
					delete tempo;
				}
				if(Sr[k] == 0) {
					if(sr0) {
						failure = true;
						cerr << "[EM] Cluster " << k << " eliminated during iteration " << t << "." << endl;
						break;
					} else {
						sr0 = true;
						failure = true;
						break;
					}
				} else if(debuggmm) {
					cerr << "Sr " << k << ": " << Sr[k] << endl;
				}
				*(means[k]) /= Sr[k];
			}

			// variance:

			foreach(k, K) {
				vars[k]->clear();

				foreach(j, (size_t)recc) {
					if(tty == PROTEIN && p_tr[mp_prefix_length][j] == NULL) continue;
					Table* tempo = (*(data[j])) - *(means[k]);
					//cout << "Step 1:<br>";
					//tempo->dump();
					Table* tempo2 = *tempo * *tempo;
					delete tempo;
					//cout << "Step 2:<br>";
					//tempo->dump();
					tempo = *tempo2 * resps[k][j];
					delete tempo2;
					//delete &tempo2;
					//cout << "Step 3:<br>";
					//tempo->dump();

					*(vars[k]) += *tempo;
					//cout << "Step 4:<br>";
					//tempo->dump();

					delete tempo;
					//exit(99);
				}
				*(vars[k]) /= Sr[k];
			}

			// update weights (and denominators for maximization):

			foreach(k, K) {
				Sr[k] = 0;
				foreach(i, (size_t)recc) {
					if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
					Sr[k] += resps[k][i];
				}
				if(Sr[k] == 0) {
					if(sr0) {
						failure = true;
						cerr << "    Cluster " << k << " eliminated after iteration " << t << "." << endl;
						break;
						if(debuggmm) {
							foreach(k, K) {
								cerr << k << "======================================================" << endl;
								foreach(i, (size_t)recc) {
									if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
									cerr << i << ": " << resps[k][i] << "\t";
								}
							}
						}
					} else {
						sr0 = true;
						failure = true;
						break;
					}
				}
				if(isnan(Sr[k])) {
					cerr << "Sr " << k << " is NaN!" << endl;
					exit(71);
				}
				if(tty == PROTEIN)
					w[k] = Sr[k] / (recc - dead_recc_p);
				else
					w[k] = Sr[k] / recc;
				if(isnan(w[k]) || isinf(w[k])) {
					cerr << "[EM] Magic weight!" << endl;
					exit(72);
				}
			}

			if(sr0) break;

			// calculate probability of each element under each component and re-normalize:
			// (expectation step)

			number T[recc];
			number preT[K];

			foreach(i, (size_t)recc) {
				if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
				foreach(k, K) {
					resps[k][i] = data[i]->loggaussprob(*(means[k]), *(vars[k]));
					if(isnan(resps[k][i]) || isinf(resps[k][i])) {
						cerr << "[EM] " << myname << ": distribution collapsed in step " << t << " of attempt " << (attempto + 1) << "/" << MAX_ATTEMPTS << endl;
						failure = true;
						break;
					}
					resps[k][i] = resps[k][i] + w[k]; // apply prior
					preT[k] = resps[k][i];

					if(isinf(preT[k])) {
						cerr << "[EM] Responsibility inf." << endl;
						cerr << "[EM] k = " << k << endl;
						failure = true;
						break;
					}

					if(isnan(preT[k])) {
						cerr << "[EM] Responsibility nan." << endl;
						cerr << "[EM] k = " << k << endl;
						failure = true;
						break;
					}
				}

				if(failure) break;

				T[i] = logsumexp(preT, (long int)K);
				if(isnan(T[i]) || isinf(T[i])) {
					cerr << "[EM] Logsumexp error." << endl;
					failure = true;
					break;
				}

				foreach(k, K) {
					resps[k][i] = exp(resps[k][i] - T[i]);
	/*
					if(resps[k][i] == 0) {
						cerr << "Responsibility zero." << endl;
						cerr << "k=" << k << endl;
						cerr << "i=" << i << endl;
						cerr << "preT=" << preT[k] << endl;
						cerr << "T=" << T[i] << endl;
						exit(73);
					}*/
				}
			}

			if(failure) break;

			if(debuggmm) {
				foreach(k, K) {
					cout << "Sr: " << Sr[k];
					cout << k << ". mean. ";
					means[k]->dump();
					cout << endl;

					cout << k << ". variance. ";
					vars[k]->dump();
					cout << endl;
		/*
					cout << k << ". resps. ";
					foreach(i, (size_t)recc) {
						if(tty == PROTEIN && ptransitions[i] == NULL) continue;
						cout << "<div title='" << (number)(resps[k][i]) << "' ";
						printf("style='background: rgb(%7.0f, %7.0f, %7.0f); width: 4px; height: 4px; display: inline-block'></div>", resps[k][i] * 255.0, (number)(resps[k][i]) * 255.0, (number)(resps[k][i]) * 255.0);
					}
					cout << "<br>";*/
				}

				cout << "\n<h1>Step " << t << "</h1>";

				foreach(o, 8) {
					int idx1 = o / 4;
					int idx2 = o % 4;
					int idy1 = o / 4 + 2;
					int idy2 = o % 4;
					cout << "<div style='background: black; width: 500px; height: 500px; position: relative; display: inline-block;'>\n";

					const number inflate = 1;
					foreach(i, (size_t)recc) {
						if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
						cout << "<div style='position: absolute; left: " << data[i]->mydata[idx1][idx2] * 500 * inflate << "px; background: ";
						printf("rgb(%d, %d, %d)", int(resps[0][i] * 255), int(resps[1][i] * 255), int(resps[2][i] * 255));
						cout << "; display: block;";
						cout << "top: " << data[i]->mydata[idy1][idy2] * 500 * inflate << "px; width: 2px; height: 2px; opacity: 1'></div>\n";
					}
					foreach(k, K) {
						cout << "<div style='position: absolute; left: " << means[k]->mydata[idx1][idx2] * 500 * inflate << "px; color: black; background: ";
						if(k == 0) cout << "red";
						if(k == 1) cout << "green";
						if(k == 2) cout << "blue";
						cout << ";";
						cout << "top: " << means[k]->mydata[idy1][idy2] * 500 * inflate << "px;'>" << k << "</div>\n";
					}

					cout << "</div>\n";
				}
			}

			llh = 0;

			foreach(i, (size_t)recc) {
				if(tty == PROTEIN && p_tr[mp_prefix_length][i] == NULL) continue;
				llh += T[i];
			}

			if(tty == PROTEIN)
				llh /= (recc - dead_recc_p);
			else
				llh /= recc;

			// test for convergence:
			if((llh - last_llh) < 10e-8 * llh && (llh - last_llh) > -10e-8 * llh)
				break;
			else
				last_llh = llh;

			if(isnan(llh)) {
				cerr << "LLH error." << endl;
				failure = true;
				break;
			}
		}

		if(!failure) {
			if(llh > bestllh || bestllh == -1) {
				cerr << "[EM] Suggesting model with llh = " << llh << " in " << t << " steps. (" << steps_so_far << " step(s) spent on " << myname << ")" << endl;
				bestllh = llh;
				foreach(k, K) {
					bestmeans[k]->clear();
					*(bestmeans[k]) += *(means[k]);
					bestvars[k]->clear();
					*(bestvars[k]) += *(vars[k]);
				}

				final_t = t;
			}
		}
		attempto = attempts;
	}

	size_t modes_found = 0;

	foreach(k, K) {
		delete means[k];
		delete vars[k];
		means[k] = bestmeans[k];
		vars[k] = bestvars[k];

		if(Sr[k] < 10e-8) {
			// validity[k] = true; // testing
			validity[k] = false;
		} else {
			validity[k] = true;
			modes_found++;
		}
	}

	if(bestllh == -1) {
		cerr << "[EM] " << myname << ": "
		     << "Failed to generate " << (tty == PROTEIN ? "AA" : "NT" ) << " Gaussian mixture model." << endl
		     << (attempto + 1) << " attempts and " << steps_so_far << " total steps." << endl;
	} else {
		cerr << "[EM] " << myname << ": "
			 << "Generated " << (tty == PROTEIN ? "AA" : "NT" ) << " Gaussian mixture model; log-likelihood " << bestllh << " after " << final_t << " steps and " << (attempto+1) << " attempt(s). (" << steps_so_far << " steps total)" << endl
			 << "[EM] Distribution has " << modes_found << " valid mode(s)." << endl;
	}
}

void Class::makeCodelta(tableType tty) {

	DifferenceTable* cod;
	Table *mu, *va;

	size_t prefix_length;

	size_t rows, cols, base;

	if(tty == PROTEIN) {
		prefix_length = cp_prefix_length;
		mu = avgtp;
		va = vartp;

		if(debugdump) cout << "Protein codelta matrix" << endl;
	} else if(tty == NUCLEOTIDE) {
		prefix_length = cn_prefix_length;
		mu = avgtn;
		va = vartn;

		if(debugdump) cout << "Nucleotide codelta matrix" << endl;
	} else {
		return;
	}

	cod = new DifferenceTable(prefix_length, tty);

	cols = mu->cols;
	base = mu->base;
	rows = (int)pow(base, prefix_length);

	size_t xcoord, ycoord;
	number xel, yel;
	// number xmu, ymu;

	derands = 0;

	number maxchance = mu->loggaussprob(*mu, *va);

	// cout << "Maxchance: " << maxchance << endl;

	size_t halfsize = cod->halfsize;

	#ifdef DOG_EAR
		size_t rcons = (rows / base) * (cols / 2);
	#endif

	foreach(el, (size_t)recc) {

		Table* bfocus = NULL;
		Table* cfocus = NULL;

		if(tty == PROTEIN) bfocus = p_tr[bp_prefix_length][el];
		if(tty == NUCLEOTIDE) bfocus = n_tr[bn_prefix_length][el];
		if(tty == PROTEIN) cfocus = p_tr[cp_prefix_length][el];
		if(tty == NUCLEOTIDE) cfocus = n_tr[cn_prefix_length][el];
		if(bfocus == NULL || cfocus == NULL) continue;

		number elchance;

		elchance = bfocus->loggaussprob(*mu, *va) / maxchance;

		// cout << "Elchance: " << elchance << endl;

		foreach(i1, rows) {
			if(tty != NUCLEOTIDE) adjustRow(i1);
			if(i1 >= rows) break;
			foreach(j1, cols) {

				ycoord = Coord(i1, j1);
				if(ycoord >= halfsize) goto stoppit_mc1;

				foreach(i2, rows) {
					if(tty != NUCLEOTIDE) adjustRow(i2);
					if(i2 >= rows) break;
					foreach(j2, cols) {

						xcoord = Coord(i2, j2);

						#ifdef DOG_EAR
							bool alt = (ycoord > xcoord);

							yel = cfocus->mydata[i1 + (alt ? rcons : 0) ][j1];
							xel = cfocus->mydata[i2 + (alt ? rcons : 0) ][j2];
						#else
							yel = cfocus->mydata[i1][j1];
							xel = cfocus->mydata[i2][j2];
						#endif

						cod->mydata[xcoord][ycoord] += pow(xel - yel, 2) * sgn(yel - xel) * elchance;
					}
				}
			}
		}
		stoppit_mc1:
		continue;
	}

	if(derands > 0) cerr << "[mcd] Warning: derandomized " << derands << " ambiguous nucleotides during class construction." << endl;

	// generated un-normalized codelta table. need to divide by recc for actual cod table.

	number maxcod = 0;
	number meancod = 0;
	number mincod = 0;
	size_t hits = 0;

	foreach(xcoord, cod->size) {
		foreach(ycoord, halfsize) {
			if(tty == PROTEIN)
				cod->mydata[xcoord][ycoord] = cod->mydata[xcoord][ycoord] / (number)(recc - dead_recc_p);
			else
				cod->mydata[xcoord][ycoord] = cod->mydata[xcoord][ycoord] / (number)(recc);

			if(cod->mydata[xcoord][ycoord] > maxcod) maxcod = cod->mydata[xcoord][ycoord];
			if(cod->mydata[xcoord][ycoord] < mincod) mincod = cod->mydata[xcoord][ycoord];
			meancod = meancod + cod->mydata[xcoord][ycoord];
			// if(debugdump) printf("%1.4f ", cod->mydata[xcoord][ycoord]);
			hits++;
		}
	}

	// generated actual cod table with some stats.

	// if(debugdump) cout << endl;

	meancod /= (number)hits;

	//cout << "Maximum codelta: " << maxcod << endl;
	//cout << "Minimum codelta: " << mincod << endl;
	//cout << "Mean codelta: " << meancod << endl;

	if(tty == PROTEIN) this->pcodelta = cod;
	if(tty == NUCLEOTIDE) this->ncodelta = cod;
}

void Class::closestDistance(Table *nT, Table* nrcT, Table *pT, number& mindeltaN, number& mindeltaNRC, number& mindeltaP) {
	mindeltaN = -1;
	mindeltaNRC = -1;
	mindeltaP = -1;
	number thisdelta;

	if(ADD_METHOD[NT_NN]) foreach(i, (size_t)recc) {
		thisdelta = n_tr[nn_prefix_length][i]->delta(nT);
		if(thisdelta < mindeltaN || mindeltaN == -1) mindeltaN = thisdelta;
		if(thisdelta == 0) break;
	}

	if(ADD_METHOD[NT_NN_RC]) foreach(i, (size_t)recc) {
		thisdelta = n_tr[nn_prefix_length][i]->delta(nrcT);
		if(thisdelta < mindeltaNRC || mindeltaNRC == -1) mindeltaNRC = thisdelta;
		if(thisdelta == 0) break;
	}

	if(ADD_METHOD[AA_NN] && pT != NULL) foreach(i, (size_t)recc) {
		thisdelta = p_tr[np_prefix_length][i]->delta(pT);
		if(thisdelta < mindeltaP || mindeltaP == -1) mindeltaP = thisdelta;
		if(thisdelta == 0) break;
	}

	mindeltaN = sqrt(mindeltaN);
	mindeltaNRC = sqrt(mindeltaNRC);
	mindeltaP = sqrt(mindeltaP);
}

number* Class::score(size_t row) {
	number* returns = new number[MAX_VOTE];
	number nscore = 0, pscore = 0, nrcscore = 0, ndist = 0, pdist = 0, nrcdist = 0, ncovs = 0, pcovs = 0, nrccovs = 0,
			nsparsescore = 0, nrcsparsescore =0, psparsescore = 0;
/*
	cerr << "BN PREFIX LENGTH: " << bn_prefix_length << endl;
	cerr << "PREFIX LENGTH OF READ in BN: " << data_n[bn_prefix_length][row]->prefix_length << endl;
	cerr << "AVGTN PREFIX LENGTH: " << avgtn->prefix_length << endl;
*/

	if(ADD_METHOD[NT_NB]) if(data_n == NULL) throw string("The whole nucleotide read table is missing.");
	if(ADD_METHOD[NT_NB]) if(data_n[bn_prefix_length] == NULL) throw string("Nucleotide Bayes read table is missing.");
	if(ADD_METHOD[NT_NB]) if(data_n[bn_prefix_length][row] == NULL) throw string("Nucleotide Bayes read table for this read is missing.");

	if(ADD_METHOD[NT_NB]) nscore = data_n[bn_prefix_length][row]->loggaussprob(*avgtn, *vartn);
	if(ADD_METHOD[NT_CD]) ncovs = gradeCodelta(data_n[cn_prefix_length][row]);

	if(ADD_METHOD[AA_NB] && data_p[bp_prefix_length][row] != NULL) pscore = data_p[bp_prefix_length][row]->loggaussprob(*avgtp, *vartp);
	if(ADD_METHOD[AA_CD] && data_p[cp_prefix_length][row] != NULL) pcovs = gradeCodelta(data_p[cp_prefix_length][row]);

	if(ADD_METHOD[NT_NB_RC]) nrcscore = data_n_rc[bn_prefix_length][row]->loggaussprob(*avgtn, *vartn);
	if(ADD_METHOD[NT_CD_RC]) nrccovs = gradeCodelta(data_n_rc[cn_prefix_length][row]);

	if(ADD_METHOD[AA_NN] || ADD_METHOD[NT_NN] || ADD_METHOD[NT_NN_RC])
		closestDistance(data_n[nn_prefix_length][row], data_n_rc[nn_prefix_length][row], data_p[np_prefix_length][row], ndist, nrcdist, pdist);

	if(ADD_METHOD[NT_SB]) nsparsescore = sparse_data_n[row]->logmultiprob(nsparse);
	if(ADD_METHOD[NT_SB_RC]) nrcsparsescore = sparse_data_n_rc[row]->logmultiprob(nsparse);
	if(ADD_METHOD[AA_SB] && sparse_data_p[row] != NULL) psparsescore = sparse_data_p[row]->logmultiprob(psparse);

	returns[NT_NB] = nscore;
	returns[AA_NB] = pscore;
	returns[NT_NB_RC] = nrcscore;
	returns[NT_NN] = ndist;
	returns[AA_NN] = pdist;
	returns[NT_NN_RC] = nrcdist;
	returns[NT_CD] = ncovs;
	returns[AA_CD] = pcovs;
	returns[NT_CD_RC] = nrccovs;
	returns[NT_MM] = 0;
	returns[AA_MM] = 0;
	returns[NT_MM_RC] = 0;
	returns[NT_SB] = nsparsescore;
	returns[AA_SB] = psparsescore;
	returns[NT_SB_RC] = nrcsparsescore;
	returns[NT_BWA] = 0;
	returns[PRIOR_16S] = 0;
	returns[PRIOR_AP] = 0;

	return returns;
}

number* Class::getMMscores(Table* datum) {
	number* returns;

	Table** means;
	Table** vars;
	bool* validity;

	size_t K;

	tableType tty = datum->tt;

	if(tty == PROTEIN) {
		means = this->mmpMeans;
		vars = this->mmpVars;
		K = mmpK;
		validity = this->mmpValid;
	} else if(tty == NUCLEOTIDE) {
		means = this->mmnMeans;
		vars = this->mmnVars;
		K = mmnK;
		validity = this->mmnValid;
	} else {
		return NULL;
	}

	returns = new number[K];

	foreach(k, K) {
		//cerr << "Calling loggaussprob for GMM cluster " << k << " in class of " << this->recc << " genes." << endl;
		if(validity[k]) {
			returns[k] = datum->loggaussprob(*(means[k]), *(vars[k]));
			if(isinf(returns[k]) && returns[k] > 0) returns[k] = 1;
			//if(returns[k] < 0) returns[k] = -10e10;
		} else {
			returns[k] = 0;
		}
		//cerr << "Called loggaussprob for GMM cluster " << k << " in class of " << this->recc << " genes." << endl;
	}

	return returns;
}

void Class::saveClass(string filename) {
	ofstream f;
	f.open(filename.c_str(), ios::out);

	f << "gist class definition " << GIST_VERSION << endl;
	f << "nucleotide-naive-bayes-prefix-length: " << bn_prefix_length << endl;
	f << "protein-naive-bayes-prefix-length: " << bp_prefix_length << endl;
	f << "nucleotide-codelta-prefix-length: " << cn_prefix_length << endl;
	f << "protein-codelta-prefix-length: " << cp_prefix_length << endl;
	f << "nucleotide-em-prefix-length: " << mn_prefix_length << endl;
	f << "protein-em-prefix-length: " << mp_prefix_length << endl;
	f << "nucleotide-row-prefix-length: " << nn_prefix_length << endl;
	f << "protein-row-prefix-length: " << np_prefix_length << endl;
	f << "gene-count: " << recc << endl;

	if(n_prefix_length > 0) {
		if(ENABLE_METHOD[NT_NB]) {
			f << "nucleotide-naive-bayes-means:";
			this->avgtn->write(f);
			f << endl;

			f << "nucleotide-naive-bayes-variance:";
			this->vartn->write(f);
			f << endl;
		}

		if(ENABLE_METHOD[NT_MM]) {
			f << "nucleotide-em-count: ";
			f << mmnK;
			f << endl;
		}

		if(ncodelta != NULL) {
			f << "nucleotide-codelta:";
			this->ncodelta->write(f);
			f << endl;
		}

		if(ENABLE_METHOD[NT_NN]) {
			foreach(i, (size_t)recc) {
				f << "nucleotide-row:";
				this->n_tr[nn_prefix_length][i]->write(f);
				f << endl;
			}
		}

		if(mmnK > 0 && ENABLE_METHOD[NT_MM]) {
			foreach(i, mmnK) {
				if(mmnValid[i]) {
					f << "nucleotide-em-means:";
					this->mmnMeans[i]->write(f);
					f << endl;

					f << "nucleotide-em-variance:";
					this->mmnVars[i]->write(f);
					f << endl;
				}
			}
		}
	}

	if(p_prefix_length > 0) {
		if(ENABLE_METHOD[AA_NB]) {
			f << "protein-naive-bayes-means:";
			this->avgtp->write(f);
			f << endl;

			f << "protein-naive-bayes-variance:";
			this->vartn->write(f);
			f << endl;
		}

		if(ENABLE_METHOD[AA_MM]) {
			f << "protein-em-count: ";
			f << mmpK;
			f << endl;
		}

		if(pcodelta != NULL) {
			f << "protein-codelta:";
			this->pcodelta->write(f);
			f << endl;
		}

		if(ENABLE_METHOD[AA_NN]) {
			foreach(i, (size_t)recc) {
				f << "protein-row:";
				if(p_tr[np_prefix_length][i] == NULL) {
					f << " null";
				} else {
					this->p_tr[np_prefix_length][i]->write(f);
				}
				f << endl;
			}
		}

		if(mmpK > 0 && ENABLE_METHOD[AA_MM]) {
			foreach(i, mmpK) {
				if(mmpValid[i]) {
					f << "protein-em-means:";
					this->mmpMeans[i]->write(f);
					f << endl;

					f << "protein-em-variance:";
					this->mmpVars[i]->write(f);
					f << endl;
				}
			}
		}
	}

	f.close();
}

void Class::loadField(string fieldname, string value) {

	static std::map<string, int> prefixes;
	if(prefixes.size() == 0) {
		prefixes[string("n-prefix-length")] = n_prefix_length;
		prefixes[string("p-prefix-length")] = p_prefix_length;
		prefixes[string("nucleotide-naive-bayes-prefix-length")] = bn_prefix_length;
		prefixes[string("protein-naive-bayes-prefix-length")] = bp_prefix_length;
		prefixes[string("nucleotide-codelta-prefix-length")] = cn_prefix_length;
		prefixes[string("protein-codelta-prefix-length")] = cp_prefix_length;
		prefixes[string("nucleotide-em-prefix-length")] = mn_prefix_length;
		prefixes[string("protein-em-prefix-length")] = mp_prefix_length;
		prefixes[string("nucleotide-row-prefix-length")] = nn_prefix_length;
		prefixes[string("protein-row-prefix-length")] = np_prefix_length;
	}

	if(prefixes.count(fieldname) == 1) {
		if(prefixes[fieldname] != atoi(value.c_str())) {
			if(prefixes[fieldname] != 0) {
				string* fruit = new string("Wrong ");
				fruit->append(fieldname);
				fruit->append(" value in class. Got ");

				ostringstream e;
				e << value;

				fruit->append(e.str());
				fruit->append(", expected ");

				e.str("");
				e << n_prefix_length;
				fruit->append(e.str());
				fruit->append(".");
				cerr << *fruit << endl;
				throw fruit;
			}
		}
	} else if(fieldname == string("gene-count")) {
		recc = atoi(value.c_str());
		//cerr << "Genes: " << recc << endl;
		if(ENABLE_METHOD[NT_NN]) this->n_tr[nn_prefix_length] = new Table*[recc];
		if(ENABLE_METHOD[AA_NN]) this->p_tr[np_prefix_length] = new Table*[recc];
	} else if(fieldname == string("nucleotide-em-count")) {
		this->mmnK = atoi(value.c_str());
		if(ENABLE_METHOD[NT_MM]) {
			mmnValid = new bool[mmnK];
			foreach(i, mmnK) {
				mmnValid[i] = false;
			}
			this->mmnMeans = new Table*[mmnK];
			this->mmnVars = new Table*[mmnK];
		}
		//cerr << "Nucleotide EM clusters: " << mmnK << endl;
	} else if(fieldname == string("protein-em-count")) {
		this->mmpK = atoi(value.c_str());
		if(ENABLE_METHOD[AA_MM]){
			mmpValid = new bool[mmpK];
			foreach(i, mmpK) {
				mmpValid[i] = false;
			}
			this->mmpMeans = new Table*[mmpK];
			this->mmpVars = new Table*[mmpK];
		}
		//cerr << "Protein EM clusters: " << mmpK << endl;
	} else if(fieldname == string("nucleotide-naive-bayes-means") && ENABLE_METHOD[NT_NB]) {
		if(bn_prefix_length != 0) this->avgtn = new Table(value, bn_prefix_length, NUCLEOTIDE);

	} else if(fieldname == string("nucleotide-naive-bayes-variance") && ENABLE_METHOD[NT_NB]) {
		if(bn_prefix_length != 0) this->vartn = new Table(value, bn_prefix_length, NUCLEOTIDE);

	} else if(fieldname == string("protein-naive-bayes-means") && ENABLE_METHOD[AA_NB]) {
		if(bp_prefix_length != 0) this->avgtp = new Table(value, bp_prefix_length, PROTEIN);

	} else if(fieldname == string("protein-naive-bayes-variance") && ENABLE_METHOD[AA_NB]) {
		if(bp_prefix_length != 0) this->vartp = new Table(value, bp_prefix_length, PROTEIN);

	} else if(fieldname == string("nucleotide-em-means") && ENABLE_METHOD[NT_MM]) {
		//cerr << "Loading Nucleotide EM mean #" << mmnKi << endl;
		if(mn_prefix_length != 0) {
			this->mmnMeans[mmnKi] = new Table(value, mn_prefix_length, NUCLEOTIDE);
			mmnValid[mmnKi] = true;
		}

	} else if(fieldname == string("nucleotide-em-variance") && ENABLE_METHOD[NT_MM]) {
		//cerr << "Loading Nucleotide EM var #" << mmnKi << " of " << mmnK << endl;
		if(mn_prefix_length != 0)
			this->mmnVars[mmnKi] = new Table(value, mn_prefix_length, NUCLEOTIDE);
		mmnKi++;

	} else if(fieldname == string("protein-em-means") && ENABLE_METHOD[AA_MM]) {
		//cerr << "Loading Protein EM mean #" << mmpKi << endl;
		if(mp_prefix_length != 0) {
			this->mmpMeans[mmpKi] = new Table(value, mp_prefix_length, PROTEIN);
			mmpValid[mmpKi] = true;
		}

	} else if(fieldname == string("protein-em-variance") && ENABLE_METHOD[AA_MM]) {
		//cerr << "Loading Protein EM var #" << mmpKi << endl;
		if(mp_prefix_length != 0)
			this->mmpVars[mmpKi] = new Table(value, mp_prefix_length, PROTEIN);
		mmpKi++;

	} else if(fieldname == string("nucleotide-codelta") && ENABLE_METHOD[NT_CD]) {
		if(cn_prefix_length != 0) this->ncodelta = new DifferenceTable(value, cn_prefix_length, NUCLEOTIDE);

	} else if(fieldname == string("protein-codelta") && ENABLE_METHOD[AA_CD]) {
		if(cp_prefix_length != 0) this->pcodelta = new DifferenceTable(value, cp_prefix_length, PROTEIN);

	} else if(fieldname == string("nucleotide-row") && ENABLE_METHOD[NT_NN]) {
		if(nn_prefix_length != 0) this->n_tr[nn_prefix_length][reci_n] = new Table(value, nn_prefix_length, NUCLEOTIDE);
		reci_n++;
	} else if(fieldname == string("protein-row") && ENABLE_METHOD[AA_NN]) {
		if(np_prefix_length != 0)  {
			if(value == string("null")) {
				dead_recc_p++;
				this->p_tr[np_prefix_length][reci_p] = NULL;
			} else {
				this->p_tr[np_prefix_length][reci_p] = new Table(value, np_prefix_length, PROTEIN);
			}
		}
		reci_p++;
	} else {
		string* fruit = new string("Unrecognized field " + fieldname);
		cerr << "WARNING: " << *fruit << endl;
		delete fruit;
		// throw fruit;
	}
}

void Class::rebuildClass(string filename) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	if(ENABLE_METHOD[AA_SB] || ENABLE_METHOD[NT_SB] || ENABLE_METHOD[NT_SB_RC]) {
		char** recv;
		//backupclock();
		recc = readFasta(filename.substr(0, filename.length() - 5), recv, false);


		if(recc == -1) {
			cerr << endl << "[rc] Couldn't open " << filename << ". Source FASTA files are required for sparse Bayes methods." << endl << endl;
			string* fruit = new string("[rc] Couldn't load class file " + filename);
			throw fruit;
		}

		for(size_t i = 0; i < (size_t)recc; i++) {
			if(ENABLE_METHOD[NT_SB] && sn_prefix_length > 0) {
				if(nsparse == NULL) {
					nsparse = new SparseTable(recv[i * 2 + 1], sn_prefix_length, NUCLEOTIDE);
				} else {
					nsparse->add_sequence(recv[i * 2 + 1]);
				}
			}

			if(ENABLE_METHOD[AA_SB] && sp_prefix_length > 0) {
				char* trans = translate(recv[i * 2 + 1], true);
				size_t translen = strlen(trans);

				if(translen > sp_prefix_length) {
					if(psparse == NULL) {
						psparse = new SparseTable(trans, sp_prefix_length, PROTEIN);
					} else {
						psparse->add_sequence(trans);
					}
				}

				delete trans; trans = NULL;
			}
		}

		foreach(i, (size_t)recc) {
			delete[] recv[i]; recv[i] = NULL;
		}
		delete[] recv; recv = NULL;

		//restoreclock();
	}

	if(f.is_open()) {

		if(max_n_prefix_length > 0) n_tr = new Table**[max_n_prefix_length + 1];
		if(max_p_prefix_length > 0) p_tr = new Table**[max_p_prefix_length + 1];

		foreach(K, max_n_prefix_length + 1)
			n_tr[K] = NULL;

		foreach(K, max_p_prefix_length + 1)
			p_tr[K] = NULL;

		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find(":");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					string value = s.substr(offset + 1);
					//cerr << "Loading " << filename << " --> " << fieldname << " (" << value.length() << " bytes)" << endl;
					this->loadField(fieldname, value);
				} else if(!disk_mode) {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}
	} else {
		string* fruit = new string("[rc] Couldn't load class file " + filename);
		cerr << *fruit << endl;
		throw fruit;
	}

}

Class::Class(string filename) { // from fasta
	avgtn = avgtp = vartn = vartp = NULL;
	ncodelta = pcodelta = NULL;
	nsparse = psparse = NULL;
	reci_n = 0;
	reci_p = 0;
	mmnKi = 0;
	mmpKi = 0;
	mmnValid = NULL;
	mmpValid = NULL;
	myname = filename;
	dead_recc_p = 0;

	indexBWA(myname);

	// Read in data:

	char** recv;
	recc = readFasta(filename, recv, false);

	if(recc == -1) {
		cerr << "[cc] Couldn't open " << filename << "." << endl;
	} else if(recc == -2) {
		this->rebuildClass(filename);
		return;
	} else {
		cerr << "[cc] " << recc << " records in " << filename << endl;
	}

	cerr << "[cc] Allocating transition tables for " << filename << endl;

	if(max_n_prefix_length > 0) n_tr = new Table**[max_n_prefix_length + 1];
	if(max_p_prefix_length > 0) p_tr = new Table**[max_p_prefix_length + 1];

	size_t ntc = 0, ptc = 0;

	foreach(i, max_n_prefix_length + 1) {
		if(n_prefix_enabled[i]) {
			n_tr[i] = new Table*[recc];
			ntc++;
		} else {
			n_tr[i] = NULL;
		}
	}

	foreach(i, max_p_prefix_length + 1) {
		if(p_prefix_enabled[i]) {
			p_tr[i] = new Table*[recc];
			ptc++;
		} else {
			p_tr[i] = NULL;
		}
	}

	cerr << "[cc] Allocating " << ntc << " nucleotide tableset(s) and " << ptc << " protein tableset(s)." << endl;

	//if(n_prefix_length > 0) { ntransitions = new Table*[recc]; } else { ntransitions = NULL; }
	//if(p_prefix_length > 0) { ptransitions = new Table*[recc]; } else { ptransitions = NULL; }

	// Create Markov tables for each sequence:

	for(size_t i = 0; i < (size_t)recc; i++) {

		string rs = string(recv[i * 2 + 1]);

		try {
			foreach(K, max_n_prefix_length + 1) {
				if(n_tr[K] == NULL) continue;

				n_tr[K][i] = new Table(&rs, K, NUCLEOTIDE);
				if(n_tr[K][i]->length == 0) {
					cerr << "[cc] Record " << i << " in " << filename << " has 0 length." << endl;
					exit(74);
				}
				*(n_tr[K][i]) /= n_tr[K][i]->length;
				n_tr[K][i]->sqrt_tf();
			}

			char* trans = translate(recv[i * 2 + 1], true);
			size_t translen = strlen(trans);

			foreach(K, max_p_prefix_length + 1) {
				if(p_tr[K] == NULL) continue;

				if(translen > (size_t)K) {
					p_tr[K][i] = new Table(trans, K, PROTEIN);
					if(p_tr[K][i]->length == 0) {
						cerr << "[cc] Protein record " << i << " in " << filename << " has 0 length but amino acids were recovered. I'm confused!" << endl;
						cerr << "[cc] Translated sequence: " << trans << endl;
						cerr << "[cc] Original sequence: " << rs << endl;
						exit(75);
					}
					*(p_tr[K][i]) /= p_tr[K][i]->length;
					p_tr[K][i]->sqrt_tf();
				} else {
					p_tr[K][i] = NULL;
					dead_recc_p++;
				}
			}

			if(ENABLE_METHOD[NT_SB] && sn_prefix_length > 0) {
				if(nsparse == NULL) {
					nsparse = new SparseTable(recv[i * 2 + 1], sn_prefix_length, NUCLEOTIDE);
				} else {
					nsparse->add_sequence(recv[i * 2 + 1]);
				}
			}

			if(ENABLE_METHOD[AA_SB] && sp_prefix_length > 0 && translen > sp_prefix_length) {
				if(psparse == NULL) {
					psparse = new SparseTable(trans, sp_prefix_length, PROTEIN);
				} else {
					psparse->add_sequence(trans);
				}
			}

			delete trans; trans = NULL;

		} catch(string* tantrum) {
			cerr << "[cc] " << *tantrum << " in record #" << i << ":" << endl;
			cerr << "[cc] " << *(recv[i * 2]) << endl;
			cerr << "[cc] Length: " << rs.length() << " characters" << endl;
			throw tantrum;
		}
	}

	cerr << "[cc] Generated gene transition tables for " << filename << endl;

	// init for avg and var calcs:

	string emptystring = "";
	string smoothstring = "smooth";

	Table *smoothtn, *smoothtp;

	if((ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]) && bn_prefix_length > 0) {
		avgtn = new Table(&emptystring, bn_prefix_length, NUCLEOTIDE);
		vartn = new Table(&emptystring, bn_prefix_length, NUCLEOTIDE);
		smoothtn = new Table(&smoothstring, bn_prefix_length, NUCLEOTIDE);
	} else {
		avgtn = NULL;
		vartn = NULL;
		smoothtn = NULL;
	}

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) && bp_prefix_length > 0) {
		avgtp = new Table(&emptystring, bp_prefix_length, PROTEIN);
		vartp = new Table(&emptystring, bp_prefix_length, PROTEIN);
		smoothtp = new Table(&smoothstring, bp_prefix_length, PROTEIN);
	} else {
		avgtp = NULL;
		vartp = NULL;
		smoothtp = NULL;
	}

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) || (ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]))
		cerr << "[cc] Initialized mean and variance transition tables for " << filename << endl;

	// calculate averages:

	for(size_t i = 0; i < (size_t)recc; i++) {
		// cout << i << endl;
		if(bn_prefix_length > 0) *avgtn += *(n_tr[bn_prefix_length][i]);
		if(bp_prefix_length > 0 && p_tr[bp_prefix_length][i] != NULL) *avgtp += *(p_tr[bp_prefix_length][i]);
	}

	// if(n_prefix_length > 0) avgtn->dump();
	// if(p_prefix_length > 0) avgtp->dump();

	if((ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]) && bn_prefix_length > 0) {
		*avgtn += *smoothtn;
		*avgtn /= (float)recc;
	}

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) && bp_prefix_length > 0) {
		*avgtp += *smoothtp;
		if(recc > dead_recc_p) {
			*avgtp /= (float)(recc - dead_recc_p);
		} else {
			cerr << "[cc] WARNING: No protein sequences could be recovered from " << filename << "." << endl;
			exit(77);
		}
	}

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) || (ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]))
		cerr << "[cc] Calculated means for " << filename << endl;

	// calculate variance:

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) || (ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]))
		foreach(i, (size_t)recc) {
			if(bn_prefix_length > 0 && (ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD])) {
				vartn->padd(*(n_tr[bn_prefix_length][i]), 2.0);
				/*vartn->dump();
				string e; cin >> e;*/
			}

			if(bp_prefix_length > 0 && p_tr[bp_prefix_length][i] != NULL && (ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD])) {
				vartp->padd(*(p_tr[bp_prefix_length][i]), 2.0);
			}
		}

	cerr << "[cc] Calculated variances for " << filename << endl;

	if(debuggmm) {
		cout << "AVG TN:<br>";
		avgtn->dump();
		cout << "VAR TN:<br>";
		vartn->dump();
	}

	Table* avgtn2 = NULL;
	Table* avgtp2 = NULL;

	if(n_prefix_length > 0) {
		avgtn2 = avgtn->powt(2.0);

		cerr << "[cc] Calculated nucleotide mean squared for " << filename << endl;

		vartn->padd(*smoothtn, 1.0);

		*vartn /= (float)recc;

		Table *vartn_new = NULL;
		vartn_new = *vartn - *avgtn2;
		delete vartn; vartn = vartn_new;

	}

	if(p_prefix_length > 0) {
		avgtp2 = avgtp->powt(2.0);

		cerr << "[cc] Calculated protein mean squared for " << filename << endl;

		vartp->padd(*smoothtp, 1.0);

		if(recc > dead_recc_p) {
			*vartp /= (float)(recc - dead_recc_p);
		} else {
			cerr << "[cc] WARNING: No protein sequences could be recovered from " << filename << "." << endl;
			// exit(76);
		}

		Table *vartp_new = NULL;
		vartp_new = *vartp - *avgtp2;
		delete vartp; vartp = vartp_new;
	}

	foreach(i, (size_t)recc) {
		delete[] recv[i]; recv[i] = NULL;
	}
	delete[] recv; recv = NULL;

	delete avgtn2;
	delete avgtp2;
	delete smoothtn;
	delete smoothtp;

	// COD:

	if(ENABLE_METHOD[NT_CD]) cerr << "[cc] Starting NT codelta for " << filename << endl;
	if(ENABLE_METHOD[NT_CD]) this->makeCodelta(NUCLEOTIDE);
	if(ENABLE_METHOD[AA_CD]) cerr << "[cc] Starting AA codelta for " << filename << endl;
	if(ENABLE_METHOD[AA_CD]) this->makeCodelta(PROTEIN);
	if(ENABLE_METHOD[AA_CD] || ENABLE_METHOD[NT_CD]) cerr << "[cc] Finished codelta calculations for " << filename << endl;

	// MM:

	mmnVars = NULL;
	mmpVars = NULL;
	mmnMeans = NULL;
	mmpMeans = NULL;
	mmnK = NT_MM_K;
	mmpK = AA_MM_K;

	if(ENABLE_METHOD[NT_MM]) cerr << "[cc] Starting NT EM of GMM for " << filename << endl;
	if(ENABLE_METHOD[NT_MM]) this->makeMM(NUCLEOTIDE);

	if(ENABLE_METHOD[AA_MM]) cerr << "[cc] Starting AA EM of GMM for " << filename << endl;
	if(ENABLE_METHOD[AA_MM]) this->makeMM(PROTEIN);

	cerr << "[cc] Post-generation cleanup for class " << filename << endl;

	foreach(K, max_n_prefix_length + 1) {
		if(K != (size_t)nn_prefix_length && n_tr[K] != NULL && n_prefix_enabled[K]) {
			foreach(j, recc) {
				delete n_tr[K][j];
				n_tr[K][j] = NULL;
			}
			delete[] n_tr[K];
			n_tr[K] = NULL;
		}
	}

	foreach(K, max_p_prefix_length + 1) {
		if(K != (size_t)np_prefix_length && p_tr[K] != NULL && n_prefix_enabled[K]) {
			foreach(j, recc) {
				delete p_tr[K][j];
				p_tr[K][j] = NULL;
			}
			delete[] p_tr[K];
			p_tr[K] = NULL;
		}
	}

	cerr << "[cc] Done constructing class " << filename << endl << endl;
}

void Class::info() {
	if(avgtn != NULL) {
		cout << "Nucleotide avg: " << endl; avgtn->dump();
		cout << "Nucleotide var: " << endl; vartn->dump();
	}

	if(avgtp != NULL) {
		cout << "Nucleotide avg: " << endl; avgtp->dump();
		cout << "Nucleotide var: " << endl; vartp->dump();
	}
}

Class::~Class() {
	foreach(K, max_n_prefix_length + 1) {
		if(n_tr[K] != NULL) {
			foreach(i, (size_t)recc) {
				delete n_tr[K][i];
				n_tr[K][i] = NULL;
			}
		}
		delete[] n_tr[K];
	}

	foreach(K, max_p_prefix_length + 1) {
		if(p_tr[K] != NULL) {
			foreach(i, (size_t)recc) {
				delete p_tr[K][i];
				p_tr[K][i] = NULL;
			}
		}
		delete[] p_tr[K];
	}

	delete[] n_tr; n_tr = NULL;
	delete[] p_tr; p_tr = NULL;

	if((ENABLE_METHOD[NT_NB] || ENABLE_METHOD[NT_CD]) && bn_prefix_length > 0) {
		delete avgtn; avgtn = NULL;
		delete vartn; vartn = NULL;
	}

	if((ENABLE_METHOD[AA_NB] || ENABLE_METHOD[AA_CD]) && bp_prefix_length > 0) {
		delete avgtp; avgtp = NULL;
		delete vartp; vartp = NULL;
	}

	delete ncodelta; ncodelta = NULL;
	delete pcodelta; pcodelta = NULL;

	if(mn_prefix_length != 0 && ENABLE_METHOD[NT_MM]) {
		foreach(i, mmnK) {
			if(mmnValid[i]) delete mmnMeans[i];
			mmnMeans[i] = NULL;
			if(mmnValid[i]) delete mmnVars[i];
			mmnVars[i] = NULL;
		}
		delete[] mmnMeans; mmnMeans = NULL;
		delete[] mmnVars; mmnVars = NULL;
	}

	if(mp_prefix_length != 0 && ENABLE_METHOD[AA_MM]) {
		foreach(i, mmpK) {
			if(mmpValid[i]) delete mmpMeans[i];
			mmpMeans[i] = NULL;
			if(mmpValid[i]) delete mmpVars[i];
			mmpVars[i] = NULL;
		}
		delete[] mmpMeans; mmpMeans = NULL;
		delete[] mmpVars; mmpVars = NULL;
	}

	delete[] mmnValid; mmnValid = NULL;
	delete[] mmpValid; mmpValid = NULL;

	delete nsparse; nsparse = NULL;
	delete psparse; psparse = NULL;
}

