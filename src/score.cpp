/*
 	Gist score.cpp

 		manages computation of class scores for individual reads
 		see classify.cpp for final class label assignments
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

#include "score.h"
#include "Table.h"
#include "Class.h"
#include "main.h"
#include "classify.h"
#include <sstream>
#include <cstring>

extern Class **classes;
extern vector<string> *class_filenames;

struct threadparams {
	size_t start_i;
	size_t end_i;
	size_t thread_id;
	size_t max_threads;
};

size_t score_prog = 0;
size_t score_goal = 0;

void apply_MM_scores(size_t j) {
	// MM calculations:

	number* n_mm_scores_unwound = NULL;
	number* nrc_mm_scores_unwound = NULL;
	number* p_mm_scores_unwound = NULL;

	if(ADD_METHOD[NT_MM]) {
		n_mm_scores_unwound = new number[cc * NT_MM_K];

		foreach(i, cc) {
			foreach(k, NT_MM_K) {
				n_mm_scores_unwound[i * NT_MM_K + k] = scoretable_NT_MM[k][j][i];
			}
		}

		number bottom = logsumexp(n_mm_scores_unwound, cc * NT_MM_K);

		foreach(i, cc) {
			foreach(k, NT_MM_K) {
				n_mm_scores_unwound[i * NT_MM_K + k] = exp(n_mm_scores_unwound[i * NT_MM_K + k] - bottom);
				scoretable[NT_MM][j][i] += n_mm_scores_unwound[i * NT_MM_K + k];
			}
		}
	} else {
		foreach(i, cc) scoretable[NT_MM][j][i] = 0;
	}

	if(ADD_METHOD[NT_MM_RC]) {
		nrc_mm_scores_unwound = new number[cc * NT_MM_K];

		foreach(i, cc) {
			foreach(k, NT_MM_K) {
				nrc_mm_scores_unwound[i * NT_MM_K + k] = scoretable_NT_MM_RC[k][j][i];
			}
		}

		number bottom = logsumexp(nrc_mm_scores_unwound, cc * NT_MM_K);

		foreach(i, cc) {
			foreach(k, NT_MM_K) {
				nrc_mm_scores_unwound[i * NT_MM_K + k] = exp(nrc_mm_scores_unwound[i * NT_MM_K + k] - bottom);
				scoretable[NT_MM_RC][j][i] += nrc_mm_scores_unwound[i * NT_MM_K + k];
			}
		}
	} else {
		foreach(i, cc) scoretable[NT_MM_RC][j][i] = 0;
	}

	if(ADD_METHOD[AA_MM] && has_protein[j]) {
		p_mm_scores_unwound = new number[cc * AA_MM_K];

		foreach(i, cc) {
			foreach(k, AA_MM_K) {
				p_mm_scores_unwound[i * AA_MM_K + k] = scoretable_AA_MM[k][j][i];
			}
		}

		number bottom = logsumexp(p_mm_scores_unwound, cc * AA_MM_K);

		foreach(i, cc) {
			foreach(k, AA_MM_K) {
				p_mm_scores_unwound[i * AA_MM_K + k] = exp(p_mm_scores_unwound[i * AA_MM_K + k] - bottom);
				scoretable[AA_MM][j][i] += p_mm_scores_unwound[i * AA_MM_K + k];
			}
		}
	} else {
		foreach(i, cc) scoretable[AA_MM][j][i] = -1;
	}

	delete[] p_mm_scores_unwound;
	delete[] nrc_mm_scores_unwound;
	delete[] n_mm_scores_unwound;
}

void score_read(stringstream& lout, size_t j, size_t k) { // j = read, k = class
	if(classes[k] != NULL) {
		//bool use_backups[MAX_VOTE];
		//memcpy(use_backups, ADD_METHOD, MAX_VOTE); // TODO: check this - risky with bools?

		/*
		forclassifiers(i)
			if(scoretable[i][j][k] != -1 && scoretable[i] != 0)
				ADD_METHOD[i] = false; // don't repeat work!
		*/

		number* scores = classes[k]->score(j);

		forclassifiers(i) if(ADD_METHOD[i]) {
			if(i != NT_BWA)
				scoretable[i][j][k] = scores[i];
/*			lout << "Scores: " << endl;
			forclassifiers(k) {
				lout << scoretable[i][k];
			}
			lout << endl;*/
		} else scoretable[i][j][k] = -1;

		//memcpy(ADD_METHOD, use_backups, MAX_VOTE);

		delete[] scores; scores = NULL;

		if(ADD_METHOD[AA_MM] && data_p[mp_prefix_length][j] != NULL) {
			number* scs = classes[k]->getMMscores(data_p[mp_prefix_length][j]);
			foreach(i, AA_MM_K) scoretable_AA_MM[i][j][k] = scs[i];
			delete[] scs;
		}

		if(ADD_METHOD[NT_MM]) {
			number* scs = classes[k]->getMMscores(data_n[mn_prefix_length][j]);
			foreach(i, NT_MM_K) scoretable_NT_MM[i][j][k] = scs[i];
			delete[] scs;
		}

		if(ADD_METHOD[NT_MM_RC]) {
			number* scs = classes[k]->getMMscores(data_n_rc[mn_prefix_length][j]);
			foreach(i, NT_MM_K) scoretable_NT_MM_RC[i][j][k] = scs[i];
			delete[] scs;
		}

		apply_MM_scores(j);
	}
}

void score(size_t start_i, size_t end_i) {
	if(end_i > (size_t)dc) end_i = dc;
	// cout << endl << "calc from " << start_i << " to " << end_i << endl << endl;

	score_goal = dc * cc;

	string clock_unit = "calcs";
	string clock_tag = "scoring data";
	time_t clock_time = time(NULL);

	if(debugdump) {
		cout << endl << "-------------------------------------" << endl;
		cout << "Thread from " << start_i << " to " << end_i << endl;
		cout << endl << "-------------------------------------" << endl;
	}

	for(size_t j = start_i; j < end_i; j++) {
		std::stringstream lout;

		lout.str("");

		// lout << ">" << *(dv[j*2]) << "\n";

		foreach(k, cc) if(classes[k] != NULL) {
			//if(disk_mode) prog_status++;
			// if(RUNLEVEL == 1 || shortlist[k][j])
			score_read(lout, j, k);
			runClock(score_prog, score_goal, clock_time, clock_unit, clock_tag);
		}

//		if(!disk_mode) prog_status++;

		if(lout.str().length() > 0) cout << lout.str();
	}
}

bool* classmap;
size_t classmap_locked_by = 0;
bool* thread_alive;

void* runscorethread(void* t) {
	threadparams* ts = (threadparams*)t;

	string clock_unit = "calcs";
	string clock_tag = "scoring";
	time_t clock_time = time(NULL);

	// cout << endl << "calc thread from " << ts->start_i << " to " << ts->end_i << " with thread ID " << ts->thread_id << endl << endl;

	try {
		if(disk_mode) {

			// target selection ==================================================================================

			next_target:
			while(classmap_locked_by != ts->thread_id && thread_alive[classmap_locked_by]) usleep(100);
			if(classmap_locked_by != ts->thread_id) {
				classmap_locked_by++;
				if(classmap_locked_by >= ts->max_threads) classmap_locked_by = 0;
				goto next_target;
			}

			/*cout << endl << "[" << ts->thread_id << "] Entered critical section." << endl;
			cout.flush();

			cout << "[" << ts->thread_id << "] Status report: ";
			foreach(ia, cc) {
				cout << (classmap[ia] ? "T" : "F");
				cout.flush();
			}
			cout << "." << endl;*/

			size_t k;
			for(k = 0; k < cc; k++) {
				if(classmap[k] == false) {
					//cout << endl << "[" << k << " is free]" << endl;
					break;
				} else {
					//cout << endl << "[" << k << " was taken]" << endl;
				}
			}
			/*
			size_t k = ((number)rand() / (number)RAND_MAX) * cc;
			if(classmap[k] == true) {

			}*/

			if(k < cc) {
				classmap[k] = true;
				//cout << endl << "[" << ts->thread_id << "] Locking class " << k << "/" << cc << endl;
				//cout.flush();
				classmap_locked_by++;
				if(classmap_locked_by >= ts->max_threads) classmap_locked_by = 0;
			} else {
				//cout << endl << "[" << ts->thread_id << "] Couldn't find work; dying. k = " << k << endl;
				//cout.flush();
				classmap_locked_by++;
				if(classmap_locked_by >= ts->max_threads) classmap_locked_by = 0;
				thread_alive[ts->thread_id] = false;
				delete ts;
				return NULL;
			}

			// begin execution ==================================================================================

			if(classes[k] == NULL) {
				//cout << endl << "[" << ts->thread_id << "] Starting class " << k << endl;
				//cout.flush();

				classes[k] = new Class(class_filenames->at(k));

				for(size_t j = ts->start_i; j < (size_t)ts->end_i; j++) {
					std::stringstream lout;

					lout.str();
					//lout << "[" << ts->thread_id << "] r" << j << ", c" << k << "/" << cc << endl;

					if(classmap_locked_by == ts->thread_id) classmap_locked_by++;
					if(classmap_locked_by >= ts->max_threads) classmap_locked_by = 0;

					score_prog++;
					// if(RUNLEVEL == 1 || shortlist[k][j])
					score_read(lout, j, k);
					runClock(score_prog, score_goal, clock_time, clock_unit, clock_tag);

//					if(!disk_mode) prog_status++;

					if(lout.str().length() > 0) cout << lout.str();
				}

				//cout << endl << "[" << ts->thread_id << "] Killing class " << k << endl;
				//cout.flush();
				delete classes[k];
				classes[k] = NULL;
			} else {
				/*cerr << endl << "[" << ts->thread_id << "] Aborted work on class " << k << " due to thread collision!" << endl;
				cout << endl << "[" << ts->thread_id << "] Aborted work on class " << k << " due to thread collision!" << endl;
				cerr.flush();
				cout.flush();*/
			}

			// end execution ==================================================================================

			goto next_target;
		} else {
			score(ts->start_i, ts->end_i);
		}

	} catch(string* e) {
		cerr << *e << endl;
		delete ts;
		return (void*)3;
	}
	// std::cerr << std::endl << "[thread] Score calculation thread finished." << std::endl;
	delete ts;
	return NULL;
}

void scorethreads(size_t threads) {
	void* status;

	classmap_locked_by = 0;

	if(disk_mode) {
		classmap = new bool[cc];
		thread_alive = new bool[threads];
		foreach(icc, cc) classmap[icc] = false;
		foreach(ith, threads) thread_alive[ith] = true;
	}

	// Threads must be "joinable" in order to wait on them to finish:
	pthread_t threadids[threads];
	pthread_attr_t attribs;
	pthread_attr_init(&attribs);
	pthread_attr_setdetachstate(&attribs, PTHREAD_CREATE_JOINABLE);

	for(size_t i = 0; i < threads; i++) {
		threadparams* thrps = new threadparams();
		if(batch_size == 0) {
			if(!disk_mode) {
				thrps->start_i = ((float)dc / (float)threads) * (float)(i);
				thrps->end_i = ((float)dc / (float)threads) * (float)(i + 1);
			} else {
				thrps->start_i = 0;
				thrps->end_i = dc;
			}
		} else {
			if(!disk_mode) {
				//thrps->start_i = ((float)dc / (float)threads) * (float)(i);
				//thrps->end_i = ((float)dc / (float)threads) * (float)(i + 1);

				size_t local_zero = current_batch * batch_size;
				size_t local_max = (current_batch + 1) * batch_size;
				size_t local_count = local_max - local_zero;

				thrps->start_i = ( (float)local_count / (float)threads ) * (float)(i) + local_zero;
				thrps->end_i = ( (float)local_count / (float)threads ) * (float)(i + 1) + local_zero;

			} else {
				thrps->start_i = current_batch * batch_size;
				thrps->end_i = (current_batch + 1) * batch_size;
				if(thrps->end_i > (size_t)dc) thrps->end_i = dc;
				//if(thrps->start_i >= thrps->end_i)
			}
		}
		thrps->thread_id = i;
		thrps->max_threads = threads;
		void* thrpsy = (void*)thrps;

		int rc = pthread_create(&threadids[i], NULL, runscorethread, thrpsy);
		if(rc) {
			std::cerr << "[head] Unable to create score calculation thread " << i << ". Condition " << rc << std::endl;
			exit(110);
		}
		usleep(100);
	}

	pthread_attr_destroy(&attribs);

	// wait on each thread to finish:
	for(size_t i = 0; i < threads; i++) {
		int rc = pthread_join(threadids[i], &status);
		if(rc) {
			std::cerr << "[head] Unable to monitor score calculation thread " << i << ". Condition " << rc << std::endl;
			exit(111);
		}
		if((unsigned long)status > 0) {
			std::cerr << "[head] Score calculation thread " << i << " exiting with status " << (unsigned long)status << std::endl;
			exit(112);
		}
	}

	delete classmap;
	classmap = NULL;
	delete thread_alive;
	thread_alive = NULL;
}

