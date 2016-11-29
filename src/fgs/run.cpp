#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hmm.h"
#include "regionals.h"

int FGS_THRESHOLD = 70;

HMM hmm;
TRAIN train;

int i, j, c, num_seq;
int wholegenome;
char mystring[1000] = "";

int count = 0;
int bp_count; /* count the length of each line in input file */

int format = 0;

int FragGeneScan_setup(const char* train_dir, const char* train_file, int wgn) {

	strcpy(mstate_file, train_dir);
	strcat(mstate_file, "gene");
	strcpy(rstate_file, train_dir);
	strcat(rstate_file, "rgene");
	strcpy(nstate_file, train_dir);
	strcat(nstate_file, "noncoding");
	strcpy(sstate_file, train_dir);
	strcat(sstate_file, "start");
	strcpy(pstate_file, train_dir);
	strcat(pstate_file, "stop");
	strcpy(s1state_file, train_dir);
	strcat(s1state_file, "stop1");
	strcpy(p1state_file, train_dir);
	strcat(p1state_file, "start1");
	strcpy(dstate_file, train_dir);
	strcat(dstate_file, "pwm");

	wholegenome = wgn;
	// format = 1;

	strcpy(hmm_file, train_dir);
	strcat(hmm_file, train_file);

	if (access(hmm_file, F_OK) == -1) {
		fprintf(stderr,
				"ERROR: The file for model parameters [%s] does not exist\n",
				hmm_file);
		//print_usage();
		exit(EXIT_FAILURE);
	}

	/* check whether the specified files exist */
	if (access(mstate_file, F_OK) == -1) {
		fprintf(stderr, "Forward prob. file [%s] does not exist\n",
				mstate_file);
		exit(1);
	}
	if (access(rstate_file, F_OK) == -1) {
		fprintf(stderr, "Backward prob. file [%s] does not exist\n",
				rstate_file);
		exit(1);
	}
	if (access(nstate_file, F_OK) == -1) {
		fprintf(stderr, "noncoding prob. file [%s] does not exist\n",
				nstate_file);
		exit(1);
	}
	if (access(sstate_file, F_OK) == -1) {
		fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
		exit(1);
	}
	if (access(pstate_file, F_OK) == -1) {
		fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
		exit(1);
	}
	if (access(s1state_file, F_OK) == -1) {
		fprintf(stderr, "start1 prob. file [%s] does not exist\n",
				s1state_file);
		exit(1);
	}
	if (access(p1state_file, F_OK) == -1) {
		fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
		exit(1);
	}
	if (access(dstate_file, F_OK) == -1) {
		fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
		exit(1);
	}
	if (access(hmm_file, F_OK) == -1) {
		fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
		exit(1);
	}

	/* read all initial model */
	hmm.N = NUM_STATE;
	get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file,
			sstate_file, pstate_file, s1state_file, p1state_file, dstate_file,
			&train);

	return 0;
}

char* FragGeneScan(char* obs_seq) {
	get_prob_from_cg(&hmm, &train, obs_seq);

	char* results;

	if (strlen(obs_seq) > (size_t)FGS_THRESHOLD) {
		viterbi(&hmm, obs_seq, &results, wholegenome, format);
		return results;
	} else {
		return NULL;
	}


/*
	// count reads
	while (fgets(mystring, sizeof mystring, fp)) {
		if (mystring[0] == '>') {
			count++;
		}
	}
	num_seq = count;
	obs_seq_len = (int *) malloc(num_seq * sizeof(int));
	printf("no. of seqs: %d\n", num_seq);

	// determine the length of each read
	i = 0;
	count = 0;
	rewind(fp);
	while (fgets(mystring, sizeof mystring, fp)) {
		if (mystring[0] == '>') {
			if (i > 0) {
				obs_seq_len[count] = i;
				count++;
			}
			i = 0;
		} else {
			bp_count = strlen(mystring) - 1;
			while (mystring[bp_count - 1] == 10 || mystring[bp_count - 1] == 13) {
				bp_count--;
			}

			i += bp_count;
		}
	}
	obs_seq_len[count] = i;

	// do the magic:
	count = -1;
	rewind(fp);
	j = 0;

	while (fgets(mystring, sizeof mystring, fp)) {

		if (mystring[0] == '>') {

			// process last-built read:
			if (count >= 0 && j > 0) {

				get_prob_from_cg(&hmm, &train, obs_seq);

				if (strlen(obs_seq) > 70) {
					viterbi(&hmm, obs_seq, fp_out, fp_aa, fp_dna, obs_head,
							wholegenome, format);
				}
			}

			bp_count = strlen(mystring) - 1;
			while (mystring[bp_count - 1] == 10 || mystring[bp_count - 1] == 13) {
				bp_count--;
			}

			obs_head = (char *) malloc((bp_count + 1) * sizeof(char));
			memset(obs_head, 0, (bp_count + 1) * sizeof(char));
			memcpy(obs_head, mystring, bp_count);

			if (count == -1 || (count >= 0 && j > 0)) {
				count++;
				obs_seq = (char *) malloc(
						obs_seq_len[count] * sizeof(char) + 1);
				memset(obs_seq, 0, obs_seq_len[count] * sizeof(char) + 1);
			}
			j = 0;

		} else {
			// add to the current read:
			bp_count = strlen(mystring) - 1;
			while (mystring[bp_count - 1] == 10 || mystring[bp_count - 1] == 13) {
				bp_count--;
			}
			memcpy(obs_seq + j, mystring, bp_count);
			j += bp_count;
		}
	}

	// pop off the last read:

	if (count >= 0) {

		get_prob_from_cg(&hmm, &train, obs_seq);

		if (strlen(obs_seq) > 70) {
			viterbi(&hmm, obs_seq, fp_out, fp_aa, fp_dna, obs_head, wholegenome,
					format);
		}
	}
*/
}

