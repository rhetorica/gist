/*
 	Gist classify.cpp

 		determines class labels, based on raw scores generated previously
 		see score.cpp for management of raw score generation
 		and Class.cpp for most of the math

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

#include "main.h"
#include "Table.h"
#include "Class.h"
#include "classify.h"
#include <cstring>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <pthread.h>
#include <sstream>
#include "autocross.h"

number** scoretable[MAX_VOTE]; // scoretable[method][data][class]
number** overall_scoretable[MAX_RUNLEVEL]; // [runlevel-1][data][class]
number*** scoretable_NT_MM; // [cluster][data][class]
number*** scoretable_NT_MM_RC; // [cluster][data][class]
number*** scoretable_AA_MM; // [cluster][data][class]
bool* has_protein; // [data]

extern bool autocross_generate;

size_t classify_status = 0;

// extern Class **classes;
// extern vector<string>* class_filenames;

/*
// Simple linked list for candidate classes:
struct ClassList {
	ClassList* next;
	Class* data;
	number* scorelist;
	void append(Class* d, number* sl) {
		next = new ClassList(d, sl);
	}
	ClassList(Class* d, number* sl) {
		data = d;
		next = NULL;
		this->scorelist = sl;
	}
	~ClassList() {
		delete next;
		delete scorelist;
		next = NULL;
		data = NULL;
		scorelist = NULL;
	}
};*/

//extern bool debuglse;

bool clock_running = false;

// prevent all display of progress clock (classify.cpp)
bool no_clock_display = false;

struct threadparams {
	size_t start_i;
	size_t end_i;
};

void* runclassifythread(void* t) {
	threadparams* ts = (threadparams*)t;
	try {
		classify(ts->start_i, ts->end_i);
	} catch(string* e) {
		cerr << *e << endl;
		delete ts;
		return (void*)3;
	}
	// std::cerr << endl << "[thread] Scoring thread finished." << std::endl;

	delete ts;
	return NULL;
}

void classifythreads(size_t threads) {
	void* status;

	// Threads must be "joinable" in order to wait on them to finish:
	pthread_t threadids[threads];
	pthread_attr_t attribs;
	pthread_attr_init(&attribs);
	pthread_attr_setdetachstate(&attribs, PTHREAD_CREATE_JOINABLE);

	for(size_t i = 0; i < threads; i++) {
		threadparams* thrps = new threadparams();
		thrps->start_i = ((float)dc / (float)threads) * (float)(i);
		thrps->end_i = ((float)dc / (float)threads) * (float)(i + 1);
		void* thrpsy = (void*)thrps;

		int rc = pthread_create(&threadids[i], NULL, runclassifythread, thrpsy);
		if(rc) {
			std::cerr << "[head] Unable to create thread " << i << ". Condition " << rc << std::endl;
			exit(100);
		}
	}

	pthread_attr_destroy(&attribs);

	// wait on each thread to finish:
	for(size_t i = 0; i < threads; i++) {
		int rc = pthread_join(threadids[i], &status);
		if(rc) {
			std::cerr << "[head] Unable to monitor scoring thread " << i << ". Condition " << rc << std::endl;
			exit(101);
		}
		if((unsigned long)status > 0) {
			std::cerr << "[head] Scoring thread " << i << " exiting with status " << (unsigned long)status << std::endl;
			exit(102);
		}
	}
}

void runClock(size_t prog_status, size_t prog_target, time_t& starttime, string& prog_unit, string& prog_tag) {
	if(clock_running) return;
	clock_running = true;

	// prevent all display of progress clock (classify.cpp)
	if(no_clock_display) return;

	ostringstream rlog;

	try {
		number progfrac = (number)prog_status / (number)prog_target;
		number oldprogfrac = (number)(prog_status - 1) / (number)prog_target;

		if((int)(progfrac * 1000) > (int)(oldprogfrac * 1000)) {
			time_t elapsed = time(NULL) - starttime;
			number speed = 0;
			/*static time_t lasttime;
			if((int)(lasttime) == (int)(elapsed)) return;
			lasttime = elapsed;*/

			if(elapsed > 0) speed = ((number)prog_status / (number)elapsed); // * (number)((number)dc / (number)prog_target);
			if(speed > 10) speed = int(speed);
			else if(speed > 1) speed = (number)((int)(speed * 10.0)) / 10.0;
			else if(speed > 0.1) speed = (number)((int)(speed * 100.0)) / 100.0;
			else speed = (number)((int)(speed * 1000.0)) / 1000.0;

			time_t total = 0;
			if(prog_status > 0) total = elapsed / progfrac;

			time_t remaining = total * (1 - progfrac);

			std::streamsize ss = rlog.precision();

			rlog << "\x1b[2K\r";
			rlog.precision(1);

			if(remaining < 0 || prog_status < 1) {
				rlog << "(calculating ETA for " << prog_tag << ")...";
			} else {
				rlog << "[" << (int)(progfrac * prog_target) << "/" << prog_target << "]" << (prog_tag.length() > 0 ? " " : "") << prog_tag << " (";
				rlog << std::fixed << ((float)((int)(progfrac * 1000)) / 10.0);
				rlog.unsetf(ios_base::floatfield);
				rlog.precision(ss);
				rlog << "%) " << displayTime(elapsed) << " elapsed, " << displayTime(remaining) << " left (" << speed << " " << prog_unit << "/sec)";
			}
		}

	} catch(void* e) {
		rlog << "\r(working)...";
	}

	cerr << rlog.str();

	clock_running = false;
}

number classify_read(std::stringstream& lout, size_t j) {

	number llh = 0;

	static number max_nn_dist = sqrt(n_prefix_length);
	static number max_cn_dist = sqrt(n_prefix_length); // ?
	static number max_np_dist = sqrt(p_prefix_length);
	static number max_cp_dist = sqrt(p_prefix_length); // ?

	number votes[cc];

	// receive results from per-class processor,
	// and finish determining scores

	number maxvote = 0;
	int maxi = -1;

	// Assemble final result set:

	number** resultset = new number*[cc];

	foreach(i, cc) {
		resultset[i] = new number[MAX_VOTE];
	}

	try {
		//lout << "Scores (i^ x k>): ";
		//forclassifiers(k) lout << "\t" << k;
		//lout << endl;
		foreach(i, cc) {
			forclassifiers(k) {
				resultset[i][k] = scoretable[k][j][i];
		//		lout << "\t" << resultset[i][k];
			}
		//	lout << endl;
		}
	} catch (string& e) {
		cerr << e;
		exit(10);
	}

	// scrub step for non-Gaussian methods.

	foreach(i, cc) {
		if(USE_METHOD[NT_CD]) if(resultset[i][NT_CD] < MINV) resultset[i][NT_CD] = MINV;
		if(USE_METHOD[NT_CD_RC]) if(resultset[i][NT_CD_RC] < MINV) resultset[i][NT_CD_RC] = MINV;
		if(USE_METHOD[AA_CD]) if(resultset[i][AA_CD] < MINV) resultset[i][AA_CD] = MINV;
		resultset[i][NT_CD] = 1 / resultset[i][NT_CD] / max_cn_dist;
		resultset[i][NT_CD_RC] = 1 / resultset[i][NT_CD_RC] / max_cn_dist;
		resultset[i][AA_CD] = 1 / resultset[i][AA_CD] / max_cp_dist;

		if(USE_METHOD[NT_NN]) if(resultset[i][NT_NN] < MINV) resultset[i][NT_NN] = MINV;
		if(USE_METHOD[NT_NN_RC]) if(resultset[i][NT_NN_RC] < MINV) resultset[i][NT_NN_RC] = MINV;
		if(USE_METHOD[AA_NN]) if(resultset[i][AA_NN] < MINV) resultset[i][AA_NN] = MINV;
		resultset[i][NT_NN] = 1 / resultset[i][NT_NN] / max_nn_dist;
		resultset[i][NT_NN_RC] = 1 / resultset[i][NT_NN_RC] / max_nn_dist;
		resultset[i][AA_NN] = 1 / resultset[i][AA_NN] / max_np_dist;

		// todo: add SB here? may be necessary as -infs are much more likely than for NB or GMM
	}

	if(autocross_generate) {
		forclassifiers(k) foreach(i, cc) scoretable[k][j][i] = resultset[i][k];
	}

	bool use_rc = false;

	if(n_prefix_length > 0) {
		number fw_score = 0;
		number rc_score = 0;

		foreach(i, cc) {
			if(USE_METHOD[NT_NB]) {
				fw_score += resultset[i][NT_NB] * bn_weight;
				rc_score += resultset[i][NT_NB_RC] * bn_weight;
			}
			if(USE_METHOD[NT_MM]) {
				fw_score += resultset[i][NT_MM] * mn_weight;
				rc_score += resultset[i][NT_MM_RC] * mn_weight;
			}
			if(USE_METHOD[NT_CD]) {
				fw_score += resultset[i][NT_CD] * cn_weight;
				rc_score += resultset[i][NT_CD_RC] * cn_weight;
			}
			if(USE_METHOD[NT_NN]) {
				fw_score += resultset[i][NT_NN] * nn_weight;
				rc_score += resultset[i][NT_NN_RC] * nn_weight;
			}
			if(USE_METHOD[NT_SB]) {
				fw_score += resultset[i][NT_SB] * sn_weight;
				rc_score += resultset[i][NT_SB_RC] * sn_weight;
			}
		}

		if(rc_score > fw_score) {
			use_rc = true;
			// lout << "Using reverse complement (llh = " << rc_score << ") over forward sequence (llh = " << fw_score << ")" << endl;
		} else {
			use_rc = false;
			// lout << "Using forward sequence (llh = " << fw_score << ") over reverse complement (llh = " << rc_score << ")" << endl;
		}
	}

	foreach(i, cc) {
		// if(!(RUNLEVEL == 1 || shortlist[i][j])) continue;

		votes[i] = 0;
		// lout << "Doing calculations for vote " << i << "\n";

		if(USE_METHOD[NT_BWA]) votes[i] += resultset[i][NT_BWA] * bwa_weight;

		if(n_prefix_length > 0) {
			if(use_rc) {
				if(USE_METHOD[NT_NB_RC]) votes[i] += resultset[i][NT_NB_RC] * bn_weight * method_weights[RUNLEVEL-1][NT_NB_RC] * (autocross_use ? autocross_MCW[NT_NB][i] : 1.0);
				if(USE_METHOD[NT_MM_RC]) votes[i] += resultset[i][NT_MM_RC] * mn_weight * method_weights[RUNLEVEL-1][NT_MM_RC] * (autocross_use ? autocross_MCW[NT_MM][i] : 1.0);
				if(USE_METHOD[NT_CD_RC]) votes[i] += resultset[i][NT_CD_RC] * cn_weight * method_weights[RUNLEVEL-1][NT_CD_RC] * (autocross_use ? autocross_MCW[NT_CD][i] : 1.0);
				if(USE_METHOD[NT_NN_RC]) votes[i] += resultset[i][NT_NN_RC] * nn_weight * method_weights[RUNLEVEL-1][NT_NN_RC] * (autocross_use ? autocross_MCW[NT_NN][i] : 1.0);
				if(USE_METHOD[NT_SB_RC]) votes[i] += resultset[i][NT_SB_RC] * sn_weight * method_weights[RUNLEVEL-1][NT_SB_RC] * (autocross_use ? autocross_MCW[NT_SB][i] : 1.0);
			} else {
				if(USE_METHOD[NT_NB]) votes[i] += resultset[i][NT_NB] * bn_weight * method_weights[RUNLEVEL-1][NT_NB] * (autocross_use ? autocross_MCW[NT_NB][i] : 1.0);
				if(USE_METHOD[NT_MM]) votes[i] += resultset[i][NT_MM] * mn_weight * method_weights[RUNLEVEL-1][NT_MM] * (autocross_use ? autocross_MCW[NT_MM][i] : 1.0);
				if(USE_METHOD[NT_CD]) votes[i] += resultset[i][NT_CD] * cn_weight * method_weights[RUNLEVEL-1][NT_CD] * (autocross_use ? autocross_MCW[NT_CD][i] : 1.0);
				if(USE_METHOD[NT_NN]) votes[i] += resultset[i][NT_NN] * nn_weight * method_weights[RUNLEVEL-1][NT_NN] * (autocross_use ? autocross_MCW[NT_NN][i] : 1.0);
				if(USE_METHOD[NT_SB]) votes[i] += resultset[i][NT_SB] * sn_weight * method_weights[RUNLEVEL-1][NT_SB] * (autocross_use ? autocross_MCW[NT_SB][i] : 1.0);
			}
		}

		if(p_prefix_length > 0 && has_protein[j]) {
			if(USE_METHOD[AA_NB]) votes[i] += resultset[i][AA_NB] * bp_weight * method_weights[RUNLEVEL-1][AA_NB] * (autocross_use ? autocross_MCW[AA_NB][i] : 1.0);
			if(USE_METHOD[AA_MM]) votes[i] += resultset[i][AA_MM] * mp_weight * method_weights[RUNLEVEL-1][AA_MM] * (autocross_use ? autocross_MCW[AA_MM][i] : 1.0);
			if(USE_METHOD[AA_CD]) votes[i] += resultset[i][AA_CD] * cp_weight * method_weights[RUNLEVEL-1][AA_CD] * (autocross_use ? autocross_MCW[AA_CD][i] : 1.0);
			if(USE_METHOD[AA_NN]) votes[i] += resultset[i][AA_NN] * np_weight * method_weights[RUNLEVEL-1][AA_NN] * (autocross_use ? autocross_MCW[AA_NN][i] : 1.0);
			if(USE_METHOD[AA_SB]) votes[i] += resultset[i][AA_SB] * sp_weight * method_weights[RUNLEVEL-1][AA_SB] * (autocross_use ? autocross_MCW[AA_SB][i] : 1.0);
		}

		if(USE_METHOD[PRIOR_16S]) votes[i] += class_16s[i] * p16s_weight * method_weights[RUNLEVEL-1][PRIOR_16S] * (autocross_use ? autocross_MCW[PRIOR_16S][i] : 1.0);
		if(USE_METHOD[PRIOR_AP]) votes[i] += class_ap[i] * ap_weight * method_weights[RUNLEVEL-1][PRIOR_AP] * (autocross_use ? autocross_MCW[PRIOR_AP][i] : 1.0);
		if(USE_METHOD[PRIOR_BIAS]) votes[i] += (autocross_use ? autocross_MCW[PRIOR_BIAS][i] : 0.0);

		llh += votes[i];
	}

	// votes are now final log scores.

	foreach(i, cc) {
		delete[] resultset[i];
	}

	delete[] resultset; resultset = NULL;

	// let's get some final scores!

	// pick new best vote:

	maxi = -1;
	maxvote = 0;

	foreach(i, cc) {
		// if(!(RUNLEVEL == 1 || shortlist[i][j])) continue;

		if(isnan(votes[i])) votes[i] = 0;

		if(votes[i] > maxvote || maxi == -1) {
			maxi = i;
			maxvote = votes[i];
		}
		overall_scoretable[RUNLEVEL-1][j][i] = votes[i];
	}

	if(maxi == -1) {
		cerr << "Couldn't generate class guess for sample " << j << "." << endl;
		exit(9);
	} else {
		// lout << "best = " << class_filenames->at(maxi) << " (class #" << maxi << ")\n";
	}

	return llh;
}

void classify(size_t start_i, size_t end_i) {
	// Classification phase:

	string clock_unit = "reads";
	string clock_tag = "classifying data";
	time_t clock_time = time(NULL);

	if(debugdump) {
		cout << endl << "-------------------------------------" << endl;
		cout << "Thread from " << start_i << " to " << end_i << endl;
		cout << endl << "-------------------------------------" << endl;
	}

	// prepare data for analysis:
	for(size_t j = start_i; j < end_i; j++) {
		std::stringstream lout;

		classify_status++;
		classify_read(lout, j);
		runClock(classify_status, dc, clock_time, clock_unit, clock_tag);

		// cout << lout.str();
	}

	// cerr << endl;
}
