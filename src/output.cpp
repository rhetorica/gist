/*
 	Gist output.cpp

 		final-layer stuff: taxonomic incorporation and actual output

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

#include "output.h"
#include "alglib/statistics.h"
#include "alglib/ap.h"
#include "classify.h"
#include <string>
#include <algorithm>

extern vector<string>* class_filenames;
extern std::map<int,string> taxnames;

std::map<int,number>* taxscores;

string RANK_NAMES[9] = {
	"species",
	"genus",
	"family",
	"order",
	"class",
	"phylum",
	"superkingdom",
	"life",
	"n/a"
};

/*
 * Given a strong hit, look at its parent taxon. If that taxon is indistinct according to a studentttest1(), then include the whole taxon.
 * Repeat for higher taxonomic levels. Take the strongest score from each child taxon when performing this integration if successful,
 * but the taxon mean if not.
 *
 * If a taxon only has one member, report just the child; do not report the taxon.
 *
 * This recurses over the tree of life.
 */

const TAXON_RANKS last_level = GENUS;

number getTaxonMeanScore(TAXON_RANKS rank, size_t taxid, size_t pcount, size_t j) {
	// cout << "gTMS call: " << RANK_NAMES[rank] << ", taxid=" << taxid << ", read=" << j << ", rl=" << RUNLEVEL << "." << endl;

	if(taxid == 0) return 0;

	size_t* pool = new size_t[pcount];
	size_t pool_count = 0;

	number meany = 0;
	size_t mc = 0;

	if(rank > last_level) {

		//cout << "All children: " << endl;

		foreach(i, cc) {
			if(rank != MAX_RANK) {
				if(taxa[i][rank] != taxid) continue; // skip non-children
			}

			if(rank > last_level) { // until we get to the last step, at which point everything is considered a child of equal level.
				foreach(k, pool_count) {
					if(taxa[i][rank-1] == taxa[pool[k]][rank-1]) goto contea; // avoid repeats
				}
			}

			// collect children:
			pool[pool_count] = i;
			pool_count++;

			/*
			if(rank > last_level)
				cout << "   " << taxa[i][rank-1] << " (" << taxnames[taxa[i][rank-1]] << ")" << endl;
			else
				cout << "   class " << i << " (" << class_filenames->at(i) << "), with score: " << overall_scoretable[RUNLEVEL-1][j][i] << endl;
*/
			meany += getTaxonMeanScore((TAXON_RANKS)(rank-1), taxa[i][rank-1], pcount, j);
			mc++;

			contea:
			continue;
		}

	} else {
		foreach(al, cc) {
			if(taxa[al][rank] == taxid) {
				meany += overall_scoretable[RUNLEVEL-1][j][al];

				mc++;
			}
		}

		if(mc == 0 || isnan(meany / (number)mc)) {
			cerr << "Agh, not again!" << endl;

			foreach(al, cc) {
				if(taxa[al][rank] == taxid) {
					meany += overall_scoretable[RUNLEVEL-1][j][al];

					mc++;
				}
			}
		}
	}

	delete[] pool; pool = NULL;

	return meany / (double)mc;

}

// never subsume a parent taxon if the mean is below this value:
#define SIGNAL_FLOOR 0

number calc_taxon_inclusion(number* return_scores, unsigned long taxid, TAXON_RANKS rank, size_t pcount, size_t j, number bloom) {
	if(pcount == 0) pcount = cc;

	number retval = 0;

//	if(rank >= MAX_RANK - 1)
		// cout << "(taxid #" << taxid << ", level " << RANK_NAMES[rank] << ", for read " << j << ", runlevel " << RUNLEVEL << "):" << endl;

	// select children to recurse over:

	//cout << "Creating..." << endl;

	size_t* child_taxon = new size_t[pcount];
	size_t* child_taxon_index = new size_t[pcount];
	size_t child_taxon_count = 0;

	size_t* pool = new size_t[pcount];
	size_t pool_count = 0;

	//cout << "All children: " << endl;

	foreach(i, cc) {
		if(rank != MAX_RANK) {
			if(taxa[i][rank] != taxid) continue; // skip non-children
		}

		if(rank > last_level) { // until we get to the last step, at which point everything is considered a child of equal level.
//			cout << "Checking for repeats..." << endl;
			foreach(k, pool_count) {
				if(taxa[i][rank-1] == taxa[pool[k]][rank-1]) goto conte; // avoid repeats
			}
		}

		// collect children:
		pool[pool_count] = i;
		pool_count++;
/*
		if(rank > last_level)
			cout << "   " << taxa[i][rank-1] << " (" << taxnames[taxa[i][rank-1]] << ")" << endl;
		else
			cout << "   class " << i << " (" << class_filenames->at(i) << "), with score: " << overall_scoretable[RUNLEVEL-1][j][i] << endl;
*/
		conte:
		continue;
	}

//	cout << "Pool size: " << pool_count << "." << endl;

	// cout << "Interesting children: " << endl;

	foreach(i, cc) {
		if(shortlist[i][j] == false) continue;

		if(rank != MAX_RANK) {
			if(taxa[i][rank] != taxid) continue; // skip non-children
		}

		if(rank > last_level) { // until we get to the last step, at which point everything is considered a child of equal level.
			foreach(k, child_taxon_count) {
				if(taxa[i][rank-1] == child_taxon[k]) goto cont; // avoid repeats
			}
		}

		child_taxon[child_taxon_count] = taxa[i][rank-1];
		child_taxon_index[child_taxon_count] = i;
		child_taxon_count++;
/*
		if(rank > last_level)
			cout << "   " << taxa[i][rank-1] << " (" << taxnames[taxa[i][rank-1]] << ")" << endl;
		else
			cout << "   class " << i << " (" << class_filenames->at(i) << "), with score: " << overall_scoretable[RUNLEVEL-1][j][i] << endl;
*/
		cont:
		continue;
	}

//	cout << "Child taxa of interest: " << child_taxon_count << "." << endl;
//	if(child_taxon_count >= 1)
//		cout << "First child taxon of interest: " << RANK_NAMES[rank-1] << " " << taxnames[child_taxon[0]] << ". (" << child_taxon[0] << ")" << endl;

	// child_taxa now contains actual children taxa.

	number* probsplus = new number[pool_count]; // overall_scoretable hits for all taxonomic subunits
	double* probs = new double[pool_count-1]; // overall_scoretable hits for all-but-one taxonomic subunits
	bool* hitflag = new bool[pool_count]; // true if child taxon contains hits
	alglib::real_1d_array rprobs;

	number localmean = 0;

	if(rank <= last_level) {
		// base case, we're looking strains. Start with shortlist results.

		number r = getTaxonMeanScore(last_level, taxid, pool_count, j);

		if(isnan(r)) {
			foreach(jl, pool_count) {
				int i = pool[jl];
				cout << "   class " << i << " (" << class_filenames->at(i) << "), with score: " << overall_scoretable[RUNLEVEL-1][j][i] << endl;
			}
			cerr << "Stop, drop, and roll: average of " << RANK_NAMES[rank] << " members is " << r << "." << endl;
		}

		size_t blocker = 0;
		if(pool_count > 1) foreach(i, pool_count) {
			if(! shortlist[pool[i]][j]) {
//				cout << "Not a candidate blocker at terminal level: class #" << i << "." << endl;
				continue;
			}

			blocker = i;

			size_t k = 0;
			foreach(eye, pool_count) {
				if(blocker == eye) {
					k--;
				} else {
					probs[k] = overall_scoretable[RUNLEVEL-1][j][pool[eye]];
				}
				k++;
			}

			double bs, ls, rs;

			rprobs.setlength(pool_count - 1);

			rprobs.setcontent(pool_count - 1, probs);

			number scoreval = overall_scoretable[RUNLEVEL-1][j][pool[blocker]];
			// alglib::ae_int_t sample_count = (alglib::ae_int_t)(child_taxon_count - 1);
			alglib::ae_int_t sample_count = (alglib::ae_int_t)(pool_count);

			alglib::studentttest1(rprobs, sample_count, scoreval, bs, ls, rs);

//			cout << "Probability that distinct hit " << pool[blocker]
//				 << " (" << class_filenames->at(pool[blocker]) << ", species #" << taxa[pool[blocker]][SPECIES] << ")"
//				 << " is distinct from parent " << RANK_NAMES[rank] << " #" << taxid << ": " << rs << endl;

			if(rs < bloom && r > SIGNAL_FLOOR ) {
				// go inclusive:
//				cout << "Subsuming parent last level taxon." << endl;
				foreach(eye, pool_count) {
					shortlist[pool[eye]][j] = true;
				}

		//		double r = overall_scoretable[RUNLEVEL-1][j][child_taxon_index[blocker]];
		//		if(r > retval) retval = r;

				break;
			} else {
				// stay exclusive:
//				cout << "Not including parent last level taxon." << endl;
			}
		}

		//if(r > retval)
		retval = r;


		if(RUNLEVEL == MAX_RUNLEVEL) {
//			cout << "Tax score for terminal group " << taxid << ": " << retval << endl;
			taxscores[j][taxid] = retval;
		} else {
//			cout << "RL = " << RUNLEVEL << endl;
		}

//		cout << "Non-CS exit of terminal CTI taxon... " << endl;

	} else {
		// recursive step, we're looking at the contents of a larger taxon. Start with return scores.

		foreach(i, pool_count) {
			bool detcalc = false;
			foreach(k, child_taxon_count) {
				if(taxa[pool[i]][rank-1] == child_taxon[k]) {
					detcalc = true;
					break;
				}
			}

			if(detcalc) {
//				cout << "Detcalc true; taxon " << taxa[pool[i]][rank-1] << " contains hits." << endl;
				probsplus[i] = calc_taxon_inclusion(NULL, taxa[pool[i]][rank-1], (TAXON_RANKS)(rank - 1), pcount, j, bloom);
			} else {
//				cout << "Detcalc false; taxon " << pool[i] << " doesn't contain hits." << endl;
				probsplus[i] = getTaxonMeanScore((TAXON_RANKS)(rank-1), taxa[pool[i]][rank-1], pcount, j);
			}

			// recursively assemble weighted average for tree branch known not to contain any hits:

			hitflag[i] = detcalc;
		}

		if(rank >= SUPERKINGDOM) {
//			cout << "Aborting due to superkingdom rank." << endl;
			goto clearstate; // Let's not get cheesy.
		}

		foreach(i, pool_count) localmean += probsplus[i];
		localmean /= (number)pool_count;

		if(isnan(localmean)) {
			cerr << "Problem in " << taxnames[taxid] << "." << endl;
			cerr << "Hits:" << endl;
			foreach(i, pool_count) cerr << probsplus[i] << " (" << taxnames[taxa[pool[i]][rank-1]] << ")" << endl;
			cerr.flush();
			cerr << "Stop, drop, and roll: average of " << RANK_NAMES[rank] << " members is " << localmean << "." << endl;
		}

		// okay, now as before:

		size_t blocker = 0;
		foreach(i, child_taxon_count) {
			if(! hitflag[i]) continue;

			blocker = i;

			size_t k = 0;
			foreach(eye, pool_count) {
				if(blocker == eye) {
					k--;
				} else {
					probs[k] = probsplus[eye];
				}
				k++;
			}

			double bs, ls, rs;

			rprobs.setcontent(pool_count - 1, probs);

			// TODO: lots of arguments about using all siblings of blocker here, or if rprobs should only contain hits that aren't blocker.
			// Pretty sure it's supposed to be pool... but y'never know.
			alglib::studentttest1(rprobs, (alglib::ae_int_t)(pool_count - 1), probsplus[blocker], bs, ls, rs);

//			cout << "Probability that " << RANK_NAMES[rank-1] << " #" << taxa[pool[blocker]][rank-1]
//				 << " (" << taxnames[taxa[pool[blocker]][rank-1]] << ")"
//				 << " is distinct from parent " << RANK_NAMES[rank] << " #" << taxid << ": " << rs << endl;

			if(rs < bloom && taxid != 0 && localmean > SIGNAL_FLOOR ) {
				// go inclusive:
//				cout << "Subsuming parent " << RANK_NAMES[rank] << "." << endl;
				foreach(al, cc) {
					if(taxa[al][rank] == taxid) shortlist[al][j] = true;
				}

				number r = probsplus[blocker];
				if(r > retval) retval = r;
				if(localmean > retval) retval = localmean;

				if(RUNLEVEL == MAX_RUNLEVEL) {
//					cout << "Storing taxscores " << taxid << " = " << retval << " (?" << r << ")" << endl;
					taxscores[j][taxid] = retval;
//					cout << "Retrieving taxscores " << taxid << " = " << taxscores[j][taxid] << " (?" << r << ")" << endl;
				} else {
//					cout << "RL = " << RUNLEVEL << endl;
				}

				break;
			} else if(taxid == 0) {
				cerr << "Not subsuming " << RANK_NAMES[rank] << " #" << taxid << ", because that is illegal." << endl;
				cerr << "Clean your data! (Found in lineage of class " << class_filenames->at(child_taxon_index[i]) << ")" << endl;

				cout << "Not subsuming " << RANK_NAMES[rank] << " #" << taxid << ", because that is illegal." << endl;
				cout << "Clean your data! (Found in lineage of class " << class_filenames->at(child_taxon_index[i]) << ")" << endl;
			} else {
				// stay exclusive:
//				cout << "Not subsuming parent " << RANK_NAMES[rank] << "." << endl;

				if(localmean > retval) retval = localmean;

				if(RUNLEVEL == MAX_RUNLEVEL) {
//					cout << "Storing taxscores " << taxid << " = " << retval << " (?" << localmean << ")" << endl;
					taxscores[j][taxid] = retval;
					// cout << "Retrieving taxscores " << taxid << " = " << taxscores[taxid] << " (?" << localmean << ")" << endl;
				} else {
//					cout << "RL = " << RUNLEVEL << endl;
				}
			}
		}

//		cout << "Non-CS exit of non-terminal CTI taxon... " << endl;

		if(isnan(localmean) || isnan(retval)) {
			cerr << "Problem in " << taxnames[taxid] << "." << endl;
			cerr << "Hits:" << endl;
			foreach(i, pool_count) { cerr << i << ": " << probsplus[i] << " (class #" << pool[i]; foreach(k, MAX_RANK) cerr << "; " << RANK_NAMES[k] << " #" << taxa[pool[i]][k] << " = " << taxnames[taxa[pool[i]][k]]; cerr << ")" << endl; }
			cerr.flush();
			cerr << "Stop, drop, and roll: average of " << RANK_NAMES[rank] << " members is " << localmean << ", retval " << retval << "." << endl;
		}
	}

//	cout << "Non-CS exit of CTI... " << endl;

	clearstate:

	if(isnan(localmean) || isnan(retval)) {
		cerr << "Problem in " << taxnames[taxid] << "." << endl;
		cerr << "Hits:" << endl;
		foreach(i, pool_count) cerr << probsplus[i] << " (" << taxnames[taxa[pool[i]][rank-1]] << ")" << endl;
		cerr.flush();
		cerr << "Stop, drop, and roll: average of " << RANK_NAMES[rank] << " members is " << localmean << ", retval " << retval << "." << endl;
	}

//	cout << "CTI cleanup... " << endl;

	delete[] child_taxon; child_taxon = NULL;
	delete[] child_taxon_index; child_taxon_index = NULL;
	delete[] pool; pool = NULL;
	delete[] probs; probs = NULL;
	delete[] probsplus; probsplus = NULL;
	delete[] hitflag; hitflag = NULL;

//	cout << "CTI cleanup done." << endl;

	return retval;
}

struct output_obj {
	size_t taxid;
	int rank; // -1 = strain
	number score;
	size_t classid;
	output_obj(size_t taxid, int rank, number score, size_t classid) : taxid(taxid), rank(rank), score(score), classid(classid) { }
};

bool sorter(const output_obj &lhs, const output_obj &rhs) { return lhs.score > rhs.score; }

void write_output() {

	bool* tokill = new bool[cc];

	if(output_mode == OUTPUT_CSV) {
		cout << "read label\tbest hit?\tstrain\tstrain score";
		for(int k = last_level; k < MAX_RANK-1; k++)
			cout << "\t" << RANK_NAMES[k] << "\t" << RANK_NAMES[k] << " score";
		cout << endl;
	}

	foreach(j, dc) {
		if(output_mode == OUTPUT_READABLE)
			cout << ">" << string(dv[j * 2]) << endl;

		bool first_hit = true;

		vector<output_obj> hits;

//		cout << "Shortlist subreport: ";
		foreach(i, cc) {
			tokill[i] = shortlist[i][j];

	//		if(tokill[i]) cout << taxnames[taxa[i][SPECIES]] << "#" << i << " (final " << overall_scoretable[RUNLEVEL-1][j][i] << ")" << ", ";
		}
		//cout << "!" << endl;

		for(int rank = MAX_RANK - 1; rank >= last_level; rank--) {
			foreach(i, cc) {
				if(tokill[i] == false) continue;

				size_t taxid = taxa[i][rank];
				size_t subt = 0; // number of first child taxon
				size_t hcount = 0; // count of active children taxa
				if(rank > last_level) subt = taxa[i][rank-1];

				bool variety = false;

				for(size_t k = i + 1; k < cc; k++) {

					if(taxid == taxa[k][rank]) {
						if(rank > last_level && taxa[i][rank-1] != subt)
							variety = true; // taxon has more than one child

						if(tokill[k] == false) // not full
							goto contol;
						else
							hcount++;

					}
				} // fallthrough: it's full

				if(rank > last_level && variety == false) continue; // don't list parent taxon if only one child taxon exists

				if(rank <= last_level && hcount <= 1) continue; // don't list parent taxon if only one strain exists

				for(size_t k = i; k < cc; k++) {
					tokill[k] = false;
				}

				hits.push_back(output_obj(taxid, rank, taxscores[j][taxid], i));

				contol:
				continue;
			}
		}

		// stragglers are unique strain hits:
		foreach(i, cc) {
			if(tokill[i]) {
				hits.push_back(output_obj(strain_taxon[i], -1, overall_scoretable[RUNLEVEL-1][j][i], i));
			}
		}

		// sort hits by score:
		sort(hits.begin(), hits.end(), sorter);

		// output time:
		foreach(q, hits.size()) {
			// number score = hits[q].score;
			int rank = hits[q].rank;
			size_t taxid = hits[q].taxid;
			size_t i = hits[q].classid;

			if(rank == -1) {
				if(output_mode == OUTPUT_READABLE) {
					// cout << "  strain \"" << class_filenames->at(i) << "\"";
					cout << "  strain #" << strain_taxon[i] << " \"" << taxnames[strain_taxon[i]] << "\"";
					cout << "\t" << overall_scoretable[RUNLEVEL-1][j][i];

					//foreach(k, MAX_RUNLEVEL)
					//	cout << "\t" << overall_scoretable[k][j][i];

					// for(int k = last_level; k < MAX_RANK-1; k++)
					//	cout << "\t" << RANK_NAMES[k] << " = " << taxscores[j][taxa[i][k]];

					cout << endl;
					// tokill[i] = false; // unnecessary
				} else if(output_mode == OUTPUT_CSV) {
					if(first_hit) {
						cout << string(dv[j * 2]) << "\t1\t";
						first_hit = false;
					} else {
						cout << string(dv[j * 2]) << "\t0\t";
					}

					cout << strain_taxon[i] << " " << taxnames[strain_taxon[i]] << "\t" << overall_scoretable[RUNLEVEL-1][j][i];

					for(int k = last_level; k < MAX_RANK-1; k++)
						cout << "\t" << taxa[i][k] << " " << taxnames[taxa[i][k]] << "\t" << taxscores[j][taxa[i][k]];

					cout << endl;
				}
			} else {
				if(output_mode == OUTPUT_READABLE) {
					cout << "  " << RANK_NAMES[rank] << " #" << taxa[i][rank] << " \"" << taxnames[taxid] << "\"";
					//for(int k = rank; k < MAX_RANK-1; k++)
					//	cout << "\t" << RANK_NAMES[k] << " = " << taxscores[j][taxa[i][k]];
					cout << "\t" << taxscores[j][taxid];
					cout << endl;
				} else if(output_mode == OUTPUT_CSV) {
					if(first_hit) {
						cout << string(dv[j * 2]) << "\t1\t";
						first_hit = false;
					} else {
						cout << string(dv[j * 2]) << "\t0\t";
					}

					cout << "\t"; // no strain

					for(int k = last_level; k < rank; k++) {
						cout << "\t\t";
					}

					for(int k = rank; k < MAX_RANK-1; k++)
						cout << "\t" << taxa[i][k] << " " << taxnames[taxa[i][k]] << "\t" << taxscores[j][taxa[i][k]];

					cout << endl;
				}
			}
		}

		hits.clear();
	}

	delete[] tokill; tokill = NULL;
}
