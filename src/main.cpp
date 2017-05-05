/*
 	Gist main.cpp

 		initialization, file i/o, and primary session management

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

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <limits>
#include <unistd.h>

#include <cstring>
#include <ctime>

#include <sstream>
#include <iostream>

#include <string>
#include <map>

#include "main.h"

// classes:
#include "Table.h"
#include "SparseTable.h"
#include "Class.h"

// program steps:
#include "build.h"
#include "score.h"
#include "classify.h"
#include "output.h"
#include "autocross.h"

// reused algorithms:
#include "fgs/use.h"
#include "bwa.h"

const bool APPEND_TYPE = false; // don't append "/genus" etc. markers on taxonomic units when loaded.

output_t output_mode;

std::map<int,string> taxnames;

bool** shortlist; // [cc][dc]
size_t* shortlist_histo; // [cc]

size_t RUNLEVEL = 0;

const char* codons = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const char* nucleotides = "TCAG";
const char* peptides = "FLSYCWPHQRIMTNKVADEG";

bool USE_METHOD[MAX_VOTE];			// methods that score() needs to pay attention to
bool ENABLE_METHOD[MAX_VOTE];		// all methods this session will consider
bool ADD_METHOD[MAX_VOTE];			// methods that classify() needs to pay attention to
bool UNWEIGHTED_METHOD[MAX_VOTE];	// priors that autocross can't optimize

bool outclasses = false;

bool disk_mode = false;
size_t batch_size = 0;
size_t current_batch = 0;

bool no_sqrt = true;

bool use_class_threads = false;

bool debugdump = false;
bool debuggmm = false;
bool debuglse = false; // used for logsumexp!

size_t n_prefix_length = 3;
size_t p_prefix_length = 1;
int num_threads = 0;

size_t bn_prefix_length = 4;
size_t bp_prefix_length = 1;

size_t sn_prefix_length = 15;
size_t sp_prefix_length = 6;

size_t cn_prefix_length = 3;
size_t cp_prefix_length = 1;

size_t nn_prefix_length = 3;
size_t np_prefix_length = 1;

size_t mn_prefix_length = 3;
size_t mp_prefix_length = 1;

size_t max_n_prefix_length = 4;
size_t max_p_prefix_length = 1;

bool* n_prefix_enabled = NULL;
bool* p_prefix_enabled = NULL;

// permanent weights:

number bwa_weight = 1;

number bn_weight = 1;
number bp_weight = 1;

number cn_weight = 1;
number cp_weight = 1;

number nn_weight = 1;
number np_weight = 1;

number mn_weight = 1;
number mp_weight = 1;

number sn_weight = 1;
number sp_weight = 1;

number ap_weight = 1;
number p16s_weight = 1;

// loaded weights:

number *class_16s = NULL;
number *class_ap = NULL;

// config file weights (per-runlevel):

number method_weights[MAX_RUNLEVEL][MAX_VOTE];

size_t NT_MM_K = 3;
size_t AA_MM_K = 2;

size_t derands = 0;
number noise_rate = 0; // synthetic noise level

size_t cc; // class count

std::string *dataname = NULL; // data file filename

int dc; // data record count
char** dv = NULL; // data

Table*** data_n = NULL;
Table*** data_n_rc = NULL;
Table*** data_p = NULL;

SparseTable** sparse_data_n = NULL;
SparseTable** sparse_data_n_rc = NULL;
SparseTable** sparse_data_p = NULL;

Class **classes;
vector<string> *class_filenames;
size_t **taxa; // [cc][rank]
size_t *strain_taxon; // [cc]

extern unsigned int class_read_length;

size_t s1_include_quota = 4; // how many strains should score phase 1 always return?
number s1_include_threshold = 0.95; // if strains have log-p above this % threshold, add them on anyway
number s1_bloom_threshold = 0.8; // if taxa have t-test p-scores below this threshold, include the rest of the taxon

size_t s2_include_quota = 4; // how many strains should score phase 2 always return?
number s2_include_threshold = 0.95; // if strains have log-p above this % threshold, add them on anyway
number s2_bloom_threshold = 0.9; // if taxa have t-test p-scores below this threshold, include the rest of the taxon

// autocross settings:

string autocross_weights_fn;

bool performance_assessment = false;

// resumption:
bool resume_mode = false;

// prevent all display of progress clock (classify.cpp)
extern bool no_clock_display;

string my_to_string(number x) {
	ostringstream ai;
	ai << x;
	return ai.str();
}

void waiter() {
	cerr << "press enter to continue";
	cin.ignore();
	cerr << endl;
}

int sgn(number val) {
    return (0 < val) - (val < 0);
}

string getexepath() {
	char result[PATH_MAX];
	ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
	string path = string(result, (count > 0) ? count : 0);
	path = path.substr(0, path.find_last_of('/'));
	return path;
}

char twobits(char e, int i) {
	if(noise_rate > 0) {
		number r = ((number)rand()) / (number)RAND_MAX;

		if(r < noise_rate) {
			return rand() % 4;
		}
	}

	switch(e) {
	case 't':
	case 'T':
	case 'u':
	case 'U':
		return 0;
	case 'c':
	case 'C':
		return 1;
	case 'a':
	case 'A':
		return 2;
	case 'g':
	case 'G':
		return 3;
	case 'y':
	case 'Y':
		derands++;
		return 0 + rand() % 2;
	case 'r':
	case 'R':
		derands++;
		return 2 + rand() % 2;
	case 'w':
	case 'W':
		derands++;
		return 0 + 2 * (rand() % 2);
	case 's':
	case 'S':
		derands++;
		return 1 + 2 * (rand() % 2);
	case 'k':
	case 'K':
		derands++;
		return 0 + 3 * (rand() % 2);
	case 'm':
	case 'M':
		derands++;
		return 1 + 1 * (rand() % 2);

	case 'v':
	case 'V':
		derands++;
		return 1 + 1 * (rand() % 3);
	case 'h':
	case 'H':
		derands++;
		return 0 + 1 * (rand() % 3);
	case 'd':
	case 'D':
		derands++;
		{
		int a = 0 + (rand() % 3);
		if(a > 0) a++;
		return a;
		}
	case 'b':
	case 'B':
		derands++;
		{
		int a = 0 + (rand() % 3);
		if(a > 1) a++;
		return a;
		}

	case 'n':
	case 'N':
		derands++;
		return rand() % 4;
	default:
		cerr << endl << endl << "Found " << (int)e << " (" << (char)e << ") during nucleotide sequence construction at position " << i << "." << endl << endl;
		string* tantrum = new string("Illegal character ");
		tantrum->append(1, e);
		throw tantrum;
	}
}

char fivebits(char e) {

	// FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG

	switch(e) {
	case 'F': return 0;
	case 'L': return 1;
	case 'S': return 2;
	case 'Y': return 3;
	case 'C': return 4;
	case 'W': return 5;
	case 'P': return 6;
	case 'H': return 7;
	case 'Q': return 8;
	case 'R': return 9;
	case 'I': return 10;
	case 'M': return 11;
	case 'T': return 12;
	case 'N': return 13;
	case 'K': return 14;
	case 'V': return 15;
	case 'A': return 16;
	case 'D': return 17;
	case 'E': return 18;
	case 'G': return 19;
	default: return -1;
	}
}

bool trust_source = false;
bool trust_class_source = true;

char* translate(char* nts, bool class_data) {
	//static size_t hit_time = 0;

	char* output;

	size_t ntslen = strlen(nts);

	if((trust_source && !class_data) || (trust_class_source && class_data)) {
		output = new char[ntslen / 3 + 1];
		output[ntslen / 3] = 0;
		for(size_t i = 0; i < ntslen - 2; i += 3) {
			char a = nts[i];
			char b = nts[i+1];
			char c = nts[i+2];

			int index = (twobits(a,i) << 4) + (twobits(b,i+1) << 2) + twobits(c, i+2);
			char aa = codons[index];
			if(aa != '*') {
				output[i / 3] = aa;
			} else {
				output[i / 3] = 0;
				break;
			}
		}
	} else {
		//char* newsy = new char[ntslen + 1];
		//strcpy(newsy, nts);
		char* outsy = FragGeneScan(nts);
		if(outsy == NULL) {
			output = new char[1];
			output[0] = 0;
		} else if(outsy[0] == '_') {
			size_t outsylen = strlen(outsy);
			output = new char[outsylen];
			strcpy(output, outsy + 1); // skip the first byte
			delete[] outsy;
		} else {
			output = outsy;
		}
		//delete[] newsy;
	}
	//cout << hit_time << " Translated " << *nts << endl << "into " << *output << endl;
	//hit_time++;
	return output;
}

char* reverse_complement(char* nrs) {
	size_t nrslen = strlen(nrs);
	char* outbuffer = new char[nrslen + 1];
	outbuffer[nrslen] = 0;
	for(size_t i = 0; i < nrslen; i++) {
		int j = nrslen - i - 1;
		switch(nrs[i] & 0b01011111) { // force uppercase
		case 'A': outbuffer[j] = 'T'; break;
		case 'C': outbuffer[j] = 'G'; break;
		case 'G': outbuffer[j] = 'C'; break;
		case 'T': outbuffer[j] = 'A'; break;
		case 'W': outbuffer[j] = 'W'; break;
		case 'S': outbuffer[j] = 'S'; break;
		case 'Y': outbuffer[j] = 'R'; break;
		case 'R': outbuffer[j] = 'Y'; break;
		case 'N': outbuffer[j] = 'N'; break;
		case 'U': outbuffer[j] = 'A'; break;
		default: outbuffer[j] = 'N'; break;
		}
	}
	//string* result = new string(outbuffer);
	//delete[] outbuffer; outbuffer = NULL;
	return outbuffer;
}

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

bool maxdump = false;

// output format: an array of C strings
// every even-numbered string is a FASTA label
// every odd-numbered string is a FASTA sequence
// return value: the number of sequences read, or a negative number if loading failed
// skim_mode: count number of records in fasta without actually loading sequence data into dv; just get sequence names
int readFasta(string filename, char**& records, bool skim_mode) {
	ifstream f;
	size_t recordcount = 0;
	string buffer = string("");
	f.open(filename.c_str(), ios::in | ios::binary);

	if(f.is_open()) {

		char linebufpre[256];

		bool in_header = false;

		// first line:

		if(f.good()) {
			f.getline(linebufpre, 256, '\n');

			if(linebufpre[0] == '>') {
				in_header = true;

				buffer.append(linebufpre);
				if(!f.fail()) {
					buffer.append("\n");
					in_header = false;
				}

				recordcount++;

				if(f.fail() && !f.eof() && !f.bad()) f.clear();
			} else {
				return -2; // Not a FASTA file
			}
		}

		// subsequent lines:

		while(!f.eof() && !f.bad()) {
			f.getline(linebufpre, 256, '\n');

			if(linebufpre[0] == '>') {
				buffer.append("\n").append(linebufpre);
				in_header = true;

				if(!f.fail()) {
					buffer.append("\n");
					in_header = false;
				}

				recordcount++;
			} else {
				if(in_header) {
					buffer.append(linebufpre);
					if(!f.fail()) {
						buffer.append("\n");
						in_header = false;
					}
				} else {
					string linebuf = linebufpre;
					scrub(linebuf);
					buffer.append(linebuf);
				}
			}

			if(f.fail() && !f.eof() && !f.bad()) f.clear();
		}
	} else {
		return -1;
	}

	time_t clock_start_time = time(NULL);
	string clock_unit = "records";
	string clock_tag = "loading data";
	size_t prog_status = 0;

	//backupclock();
	//prepclock(recordcount, "records", "loading data");

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

			if(buff2[start] == '>' || !skim_mode) {
				if(buff2[start] == '>')
					start++;

				records[t] = new char[end - start + 2];
				strncpy(records[t], (char*)(buff2 + start), end - start + 1);
				records[t][end - start + 1] = '\0';
			} else {
				records[t] = NULL; // no read data
			}

			start = i + 1;

			prog_status = t / 2;
			runClock(prog_status, recordcount, clock_start_time, clock_unit, clock_tag);
		}
	}

	f.close();

	cerr << endl << "[load] Found " << recordcount << " sequences in FASTA file " << filename << "." << endl;

	// exit(55);

	return recordcount;
}

#define override_params false
// #define forced_params "gist -data /home/rhetorica/gist-ws/bs/baby.fasta -cnw 0.5 -nnw 0 -classes /home/rhetorica/gist-ws/classes -n 3"
// #define forced_params "gist -classes /home/rhetorica/gist-mm-ws/classes -data /home/rhetorica/gist-mm-ws/bs/baby.fasta -n 1 -mnw 1"
#define forced_params "gist -noise 0 -classes /home/rhetorica/gist-mm-ws/classes -t 0 -n 1 -nnw 0 -npw 1 -cpw 1 -cnw 0 -mnw 0 -p 1 -data /home/rhetorica/gist-mm-ws/bs/baby.fasta"

int override_args(char**& args, int argc, char** argv) {
	args = argv;

	if(argc == 1 && override_params) {
		string newargs = forced_params;

		argc = 1;
		for(size_t zz = 0; zz < newargs.length(); zz++)
			if(newargs[zz] == ' ') argc++;

		// cerr << "Argc: " << argc << endl;

		args = new char*[argc];

		int argn = 0, wordstart = 0;

		for(size_t i = 0; i < newargs.length(); i++) {
			if(i == newargs.length() - 1 || newargs[i + 1] == ' ') {

				int wordlength = i - wordstart + 1;

				args[argn] = new char[wordlength + 1];
				args[argn][wordlength] = 0;
				strcpy(args[argn], newargs.substr(wordstart, wordlength).c_str());
				wordstart = i + 2;

				// cerr << "input " << argn << " received: " << args[argn] << endl;

				argn++;
			}
		}

		return argn;
	}

	return argc;
}

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

void get_paths(string& dirname, vector<string>* filenames) {
	DIR *dir;
	struct dirent *ent = NULL;
	dir = opendir(dirname.c_str());
//	struct stat buf, bufg;

	unsigned long ccnt = 0;

	if(dir) {
		while((ent = readdir(dir))) {
			if(ent->d_type == 4) continue;

			string fname = string(dirname).append("/").append(string(ent->d_name));

			string ext3 = fname.substr(fname.length() - 4, 4);
			string ext2 = fname.substr(fname.length() - 3, 3);

			if(ext2 == string(".sa")) continue;
			if(ext3 == string(".amb")) continue;
			if(ext3 == string(".ann")) continue;
			if(ext3 == string(".bwt")) continue;
			if(ext3 == string(".pac")) continue;
			if(ext3 == string(".sam")) continue;
			if(ext3 == string(".bam")) continue;
			if(ext3 == string(".sai")) continue;
			if(ext3 == string(".lin")) continue;
			if(ext3 == string(".log")) continue;
			if(ext3 == string(".txt")) continue;
			if(ext3 == string(".lst")) continue;
			if(ext3 == string(".mcw")) continue;

			// cerr << "-> " << fname << " (" << (int)(ent->d_type) << ")" << endl;

			string fname_g = string(fname).append(".gist");
			if(exists(fname_g) && is_file(fname_g)) {
				continue;
			} else if(exists(fname) && is_file(fname)) {
				// cerr << "    " << ent->d_name << endl;
				filenames->push_back(fname);
				ccnt++;
			}

		}
	} else {
		cerr << "Could not open " << dirname << "." << endl;
		exit(20);
	}

	closedir(dir);

	cerr << "Retrieved " << ccnt << " classes from directory \"" << dirname << "\"." << endl;
}

// from logspace to logspace, and dust to dust
number logsumexp(number* xs, size_t count) {
	number y = 0;
	number nxs[count];

	memcpy(nxs, xs, sizeof(number) * count);

	if(debuglse) foreach(i, count) {
		cout << "logsumexp working with " << nxs[i] << endl;
	}

	foreach(i, count) {
		if(y < nxs[i]) y = nxs[i];
	}

	if(debuglse) cout << "y = " << y << endl;

	foreach(i, count) {
		nxs[i] -= y;
		if(debuglse) cout << "nxsi = " << nxs[i] << endl;
		if(debuglse) cout << "exp nxsi = " << exp(nxs[i]) << endl;
	}

	number sum = 0;

	foreach(i, count) {
		sum += exp(nxs[i]);
	}

	sum = y + log(sum);

	if(isinf((number)sum)) {
		if(debuglse) cout << "logsum barfed: returning " << y << "!" << endl;
		sum = y;
	} else if(debuglse) {
		cout << "logsumexp returning " << sum << endl;
	}

	return sum;
}

void loadConfField(string fieldname, string value) {

	const char* c = value.c_str();

	// cout << "'" << fieldname << "' = " << "'" << value << "'" << endl;

	if(fieldname == string("#"))
		return;

	else if(fieldname == string("pass1.threshold"))
		s1_include_threshold = atof(c);
	else if(fieldname == string("pass1.quota"))
		s1_include_quota = atoi(c);
	else if(fieldname == string("pass1.bloom"))
		s1_bloom_threshold = atof(c);

	else if(fieldname == string("pass2.threshold"))
		s2_include_threshold = atof(c);
	else if(fieldname == string("pass2.quota"))
		s2_include_quota = atoi(c);
	else if(fieldname == string("pass2.bloom"))
		s2_bloom_threshold = atof(c);

	else if(fieldname == string("pass1.prior"))
		method_weights[0][PRIOR_AP] = atof(c);
	else if(fieldname == string("pass1.16s"))
		method_weights[0][PRIOR_16S] = atof(c);
	else if(fieldname == string("pass1.bwa"))
		method_weights[0][NT_BWA] = atof(c);

	else if(fieldname == string("pass1.nt.1nn"))
		method_weights[0][NT_NN] = atof(c);
	else if(fieldname == string("pass1.nt.bayes"))
		method_weights[0][NT_NB] = atof(c);
	else if(fieldname == string("pass1.nt.gmm"))
		method_weights[0][NT_MM] = atof(c);
	else if(fieldname == string("pass1.nt.codelta"))
		method_weights[0][NT_CD] = atof(c);
	else if(fieldname == string("pass1.nt.sparse"))
		method_weights[0][NT_SB] = atof(c);

	else if(fieldname == string("pass1.aa.1nn"))
		method_weights[0][AA_NN] = atof(c);
	else if(fieldname == string("pass1.aa.bayes"))
		method_weights[0][AA_NB] = atof(c);
	else if(fieldname == string("pass1.aa.gmm"))
		method_weights[0][AA_MM] = atof(c);
	else if(fieldname == string("pass1.aa.codelta"))
		method_weights[0][AA_CD] = atof(c);
	else if(fieldname == string("pass1.aa.sparse"))
		method_weights[0][AA_SB] = atof(c);


	else if(fieldname == string("pass2.prior"))
		method_weights[1][PRIOR_AP] = atof(c);
	else if(fieldname == string("pass2.16s"))
		method_weights[1][PRIOR_16S] = atof(c);
	else if(fieldname == string("pass2.bwa"))
		method_weights[1][NT_BWA] = atof(c);

	else if(fieldname == string("pass2.nt.1nn"))
		method_weights[1][NT_NN] = atof(c);
	else if(fieldname == string("pass2.nt.bayes"))
		method_weights[1][NT_NB] = atof(c);
	else if(fieldname == string("pass2.nt.gmm"))
		method_weights[1][NT_MM] = atof(c);
	else if(fieldname == string("pass2.nt.codelta"))
		method_weights[1][NT_CD] = atof(c);
	else if(fieldname == string("pass2.nt.sparse"))
		method_weights[1][NT_SB] = atof(c);

	else if(fieldname == string("pass2.aa.1nn"))
		method_weights[1][AA_NN] = atof(c);
	else if(fieldname == string("pass2.aa.bayes"))
		method_weights[1][AA_NB] = atof(c);
	else if(fieldname == string("pass2.aa.gmm"))
		method_weights[1][AA_MM] = atof(c);
	else if(fieldname == string("pass2.aa.codelta"))
		method_weights[1][AA_CD] = atof(c);
	else if(fieldname == string("pass2.aa.sparse"))
		method_weights[1][AA_SB] = atof(c);

	else if(fieldname == string("length.nt.1nn"))
		nn_prefix_length = atoi(c);
	else if(fieldname == string("length.nt.bayes"))
		bn_prefix_length = atoi(c);
	else if(fieldname == string("length.nt.gmm"))
		mn_prefix_length = atoi(c);
	else if(fieldname == string("length.nt.codelta"))
		cn_prefix_length = atoi(c);
	else if(fieldname == string("length.nt.sparse"))
		sn_prefix_length = atoi(c);

	else if(fieldname == string("length.aa.1nn"))
		np_prefix_length = atoi(c);
	else if(fieldname == string("length.aa.bayes"))
		bp_prefix_length = atoi(c);
	else if(fieldname == string("length.aa.gmm"))
		mp_prefix_length = atoi(c);
	else if(fieldname == string("length.aa.codelta"))
		cp_prefix_length = atoi(c);
	else if(fieldname == string("length.aa.sparse"))
		sp_prefix_length = atoi(c);

	else
		cerr << "Unrecognized field '" << fieldname << "' in .conf file." << endl;
}

void loadLineage(string filename, size_t class_number) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	foreach(i, MAX_RANK) taxa[class_number][i] = 0;

	bool first_line = true;

	if(f.is_open()) {
		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find(": ");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					string value = s.substr(offset + 2);
					string valname = value.substr(value.find(" ") + 1);
					if(valname[0] == '(' && valname[valname.length() - 1] == ')')
						valname = valname.substr(1, valname.length() - 2);
					size_t sval = atoi(value.substr(0, value.find(" ")).c_str());
					// cerr << "Loading " << filename << " --> " << fieldname << " = " << sval << " (class " << class_number << ")" << endl;

					if(first_line) {
						strain_taxon[class_number] = sval;
					}

					if(fieldname == string("species")) {
						taxa[class_number][SPECIES] = sval;
					} else if(fieldname == string("genus")) {
						taxa[class_number][GENUS] = sval;
					} else if(fieldname == string("family")) {
						taxa[class_number][FAMILY] = sval;
					} else if(fieldname == string("order")) {
						taxa[class_number][ORDER] = sval;
					} else if(fieldname == string("class")) {
						taxa[class_number][CLASS] = sval;
					} else if(fieldname == string("phylum")) {
						taxa[class_number][PHYLUM] = sval;
					} else if(fieldname == string("superkingdom")) {
						taxa[class_number][SUPERKINGDOM] = sval;
					} else if(!first_line) {
						continue;
					}

					first_line = false;

					if(APPEND_TYPE) {
						valname.append("/").append(fieldname);
					}
					taxnames[sval] = valname;
				} else {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}

		f.close();
	} else {
		string* fruit = new string("Couldn't load configuration file " + filename);
		cerr << *fruit << endl;
		throw fruit;
	}

	for(size_t i = 1; i < MAX_RANK; i++)
		if(taxa[class_number][i] == 0) taxa[class_number][i] = taxa[class_number][i-1]; // propagate

	/*cout << "Lineage results for class " << class_number << ":" << endl;
	foreach(i, MAX_RANK) {
		cout << "  level " << i << " is taxid #" << taxa[class_number][i] << endl;
	}*/
}

void loadConf(string filename) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	if(f.is_open()) {
		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find(" ");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					string value = s.substr(offset + 1);
					//cerr << "Loading " << filename << " --> " << fieldname << " (" << value.length() << " bytes)" << endl;
					loadConfField(fieldname, value);
				} else {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}

		f.close();
	} else {
		string* fruit = new string("Couldn't load configuration file " + filename);
		cerr << *fruit << endl;
		throw fruit;
	}

	foreach(i, MAX_RUNLEVEL) {
		method_weights[i][NT_NB_RC] = method_weights[i][NT_NB];
		method_weights[i][NT_NN_RC] = method_weights[i][NT_NN];
		method_weights[i][NT_CD_RC] = method_weights[i][NT_CD];
		method_weights[i][NT_MM_RC] = method_weights[i][NT_MM];
		method_weights[i][NT_SB_RC] = method_weights[i][NT_SB];
	}
}

void load16S(string filename) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	if(class_16s != NULL) delete class_16s;
	class_16s = new number[cc];
	foreach(i, cc) class_16s[i] = 1;

	if(f.is_open()) {
		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find("\t");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					number value = atof(s.substr(offset + 1).c_str());
					//cerr << "Loading " << filename << " --> " << fieldname << " (" << value.length() << " bytes)" << endl;
					for(size_t i = 0; i < cc; i++) {
						string cfn = class_filenames->at(i);
						if(cfn.rfind("/") != string::npos) cfn = cfn.substr(cfn.rfind("/") + 1);
						cfn = cfn.substr(0, cfn.find("."));
						if(cfn == fieldname) {
							class_16s[i] = value;
							break;
						}
					}

				} else {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}

		f.close();

		long long meansy = 0;

		foreach(i, cc) meansy += class_16s[i];
		foreach(i, cc) class_16s[i] = log(class_16s[i] / (number)meansy * ((number)cc + 1.0));
	} else {
		string* fruit = new string("Couldn't load 16S count file " + filename);
		cerr << *fruit << endl;
		throw fruit;
	}
}

void loadAPriori(string filename) {
	ifstream f;
	f.open(filename.c_str(), ios::in);

	if(class_ap != NULL) delete class_ap;
	class_ap = new number[cc];

	if(f.is_open()) {
		while(f.good()) {
			string s;
			getline(f, s);
			if(s.length() > 0) {
				int offset = s.find("\t");
				if(offset != -1) {
					string fieldname = s.substr(0, offset);
					number value = atof(s.substr(offset + 1).c_str());
					//cerr << "Loading " << filename << " --> " << fieldname << " (" << value.length() << " bytes)" << endl;
					for(size_t i = 0; i < cc; i++) {
						string cfn = class_filenames->at(i);
						if(cfn.rfind("/") != string::npos) cfn = cfn.substr(cfn.rfind("/") + 1);
						cfn = cfn.substr(0, cfn.find("."));
						if(cfn == fieldname) {
							class_ap[i] = value;
							break;
						}
					}
				} else {
					ostringstream ox;
					ox << filename << ": " << s << endl;
					cerr << ox.str();
				}
			}
		}

		f.close();
	} else {
		string* fruit = new string("Couldn't load a priori weighting file " + filename);
		cerr << *fruit << endl;
		throw fruit;
	}

	foreach(i, cc) class_ap[i] = log(class_ap[i]);
}

void clear_data_tables(size_t first, size_t stop) {
	foreach(K, max_n_prefix_length + 1) {
		if(n_prefix_enabled[K]) for(size_t j = first; j < stop; j++) {
			delete data_n[K][j];
			data_n[K][j] = NULL;
			delete data_n_rc[K][j];
			data_n_rc[K][j] = NULL;
		}
	}

	foreach(K, max_p_prefix_length + 1) {
		if(p_prefix_enabled[K]) for(size_t j = first; j < stop; j++) {
			delete data_p[K][j];
			data_p[K][j] = NULL;
		}
	}

	if(ENABLE_METHOD[NT_SB]) for(size_t j = first; j < stop; j++) delete sparse_data_n[j];

	if(ENABLE_METHOD[NT_SB_RC]) for(size_t j = first; j < stop; j++) delete sparse_data_n_rc[j];

	if(ENABLE_METHOD[AA_SB])  for(size_t j = first; j < stop; j++) delete sparse_data_p[j];
}

void gen_data_tables(size_t first, size_t stop) {

	time_t clock_start_time = time(NULL);
	string clock_unit = "tables";
	string clock_tag = "generating tables";

	size_t rec_max = stop - first;

	size_t prog_status = first;

	// cout << endl<< "generating data tables "<< first << " to " << stop << endl << endl;

	if(data_n == NULL) {
		data_n = new Table**[max_n_prefix_length+1];
		data_n_rc = new Table**[max_n_prefix_length+1];
		data_p = new Table**[max_p_prefix_length+1];

		sparse_data_n = new SparseTable*[dc];
		sparse_data_n_rc = new SparseTable*[dc];
		sparse_data_p = new SparseTable*[dc];

		foreach(K, max_n_prefix_length + 1) {
			if(n_prefix_enabled[K]) {
				data_n[K] = new Table*[dc];
				data_n_rc[K] = new Table*[dc];
			} else {
				data_n[K] = NULL;
				data_n_rc[K] = NULL;
			}
		}

		foreach(K, max_p_prefix_length + 1) {
			if(p_prefix_enabled[K]) {
				data_p[K] = new Table*[dc];
			} else {
				data_p[K] = NULL;
			}
		}
	}

	for(size_t j = first; j < stop; j++) {

		Table* datum_n = NULL;
		Table* datum_n_rc = NULL;
		Table* datum_p = NULL;

		char* nts = dv[j * 2 + 1];
		char* nrc = reverse_complement(nts);

		foreach(K, max_n_prefix_length + 1) if(data_n[K] != NULL) {
			if(n_prefix_length > 0) {
				derands = 0;

				if(debugdump) {
					cout << "Sequence: " << dv[j * 2 + 1] << endl;
					cout << "Reversed sequence: " << nrc << endl;
				}

				datum_n = new Table(nts, K, NUCLEOTIDE);
				datum_n_rc = new Table(nrc, K, NUCLEOTIDE);

				if(derands > 0) cerr << "\n[head] NOTICE: randomly filled in " << derands << " ambiguous nucleotides." << endl;

				if(datum_n->length == 0) {
					cerr << "\n[head] WARNING: read has length 0: " << string(dv[j * 2]) << endl;
				} else {
					*datum_n /= datum_n->length;
					*datum_n_rc /= datum_n_rc->length;
					datum_n->sqrt_tf();
					datum_n_rc->sqrt_tf();

					if(debugdump) {
						cout << "Normal dump: " << endl;
						datum_n->dump();

						cout << "Reversed dump: " << endl;
						datum_n_rc->dump();
					}
				}


			}

			data_n[K][j] = datum_n;
			data_n_rc[K][j] = datum_n_rc;
		}

		if(ENABLE_METHOD[NT_SB])
			sparse_data_n[j] = new SparseTable(nts, sn_prefix_length, NUCLEOTIDE);

		if(ENABLE_METHOD[NT_SB_RC])
			sparse_data_n_rc[j] = new SparseTable(nrc, sn_prefix_length, NUCLEOTIDE);

		delete[] nrc;

		char* trans;

		if(p_prefix_length > 0) {
			try {
				trans = translate(nts, false);
			} catch (void* e) {
				cerr << "panic: unable to translate " << (string)(nts) << endl;
				exit(1);
			}

			if(strlen(trans) < max_p_prefix_length + 1) {
				datum_p = NULL;
				// cout << "No protein sequence was extracted. Protein methods are disabled for this read." << endl;
				has_protein[j] = false;
			} else {
				has_protein[j] = true;
			}
		} else {
			trans = NULL;
			has_protein[j] = false;
		}

		if(has_protein[j]) {
			foreach(K, max_p_prefix_length + 1) {

				if(data_p[K] != NULL)
				{

					datum_p = new Table(trans, K, PROTEIN);

					*datum_p /= datum_p->length;
					datum_p->sqrt_tf();

					if(debugdump) datum_p->dump();

					data_p[K][j] = datum_p;
				}
			}

			if(ENABLE_METHOD[AA_SB])
				sparse_data_p[j] = new SparseTable(trans, sn_prefix_length, PROTEIN);

			delete[] trans;
			trans = NULL;
		} else {
			foreach(K, max_p_prefix_length + 1)
				if(data_p[K] != NULL)
					data_p[K][j] = NULL;

			if(ENABLE_METHOD[AA_SB]) sparse_data_p[j] = NULL;
		}

		++prog_status;

		runClock(prog_status, rec_max, clock_start_time, clock_unit, clock_tag);
	}
}

void gen_data_tables() {
	gen_data_tables(0, dc);
}

// overall_scoretable[runlevel-1][data][class]
// shortlist[class][data]

int MAXQUOTARATIO[MAX_RUNLEVEL] = {2, 12};

void do_shortlist(size_t include_quota, number include_threshold, TAXON_RANKS tr, number bloom, string & clock_tag) {
	time_t clock_start_time = time(NULL);
	string clock_unit = "reads";
	size_t prog_status = 0;

	if(include_quota > cc / MAXQUOTARATIO[RUNLEVEL-1]) {
		cerr << "[head] WARNING: inclusion quota for runlevel " << RUNLEVEL << " is too high: " << include_quota << " exceeds class count / " << MAXQUOTARATIO[RUNLEVEL-1] << "." << endl;
		include_quota = cc / MAXQUOTARATIO[RUNLEVEL-1]; // otherwise, hilarity ensues. Badly.
		if(include_quota == 0) include_quota = 1;
		cerr << "It has been rounded down to " << include_quota << " for utility's sake." << endl;

	}

	foreach(d, (size_t)dc) {
		number bestscore = 0;
		int bestclass = -1;

		// find best hit:
		foreach(c, cc) {
			if(bestscore < overall_scoretable[RUNLEVEL-1][d][c] || bestclass == -1) {
				bestscore = overall_scoretable[RUNLEVEL-1][d][c];
				bestclass = c;
			}
		}

		// cout << "************** Best hit for read " << d << ": " << bestclass << " (score: " << bestscore << ")" << endl;

		// add hits above threshold:
		size_t hits = 0;
		foreach(c, cc) {
			if(bestscore * include_threshold <= overall_scoretable[RUNLEVEL-1][d][c]) {
				shortlist[c][d] = true;
				hits++;
			}
		}

		// cout << "Hits above threshold (" << include_threshold << "): " << hits << endl;

		// pad to meet quota:
		while(hits < include_quota) {
			number sbestscore = 0;
			int sbestclass = -1;
			foreach(c, cc) if(shortlist[c][d] == false) {
				if(sbestscore < overall_scoretable[RUNLEVEL-1][d][c] || sbestclass == -1) {
					sbestscore = overall_scoretable[RUNLEVEL-1][d][c];
					sbestclass = c;
				}
			}
			if(sbestclass == -1) break; // should only be possible when everything has already been included.
			shortlist[sbestclass][d] = true;

			// lie about meeting quota if we resort to pulling zeroes:
			if(overall_scoretable[RUNLEVEL-1][d][sbestclass] == 0) shortlist[sbestclass][d] = false;
			hits++;
		}

		// cout << "Final hits after quota (" << include_quota << ") padding: " << hits << endl;



	/*
		// New version in calc_taxon_inclusion() decides when this is statistically appropriate.
		// A very sharp, single-strain signal discredits the idea of re-including a large taxon.

		// expand to include whole taxa:

		if(tr != NO_RANK) {
			bool hitlist[cc]; // things we just added; cuts down on complexity
			foreach(c, cc) hitlist[c] = false;

			foreach(c, cc) if(shortlist[c][d] && !hitlist[c]) {
				foreach(e, cc) if(!shortlist[e][d]) {
					if(taxa[c][tr] == taxa[e][tr]) {
						shortlist[e][d] = true;
						hitlist[e] = true;
						hits++;
					}
				}
			}
		}

		cout << "Final hits after taxon inhalation: " << hits << endl;*/

		// expand to include whole taxon:

		calc_taxon_inclusion(NULL, 0, MAX_RANK, 0, d, bloom);

		prog_status++;

		runClock(prog_status, dc, clock_start_time, clock_unit, clock_tag);
	}

	// shortlist histogram:
	shortlist_histo = new size_t[cc];
	foreach(i, cc) {
		shortlist_histo[i] = 0;
		foreach(j, (size_t)dc) if(shortlist[i][j]) shortlist_histo[i]++;
	}
}

void scoreprinter() {
	cout << "runlevel " << RUNLEVEL;
	foreach(i, cc) {
		cout << "\t" << taxnames[strain_taxon[i]];
	}
	cout << endl;
	foreach(j, (size_t)dc) {
		cout << string(dv[j * 2]) << "\t";
		foreach(i, cc) {
			cout << overall_scoretable[RUNLEVEL-1][j][i] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void scoretableprinter(size_t k) {
	cout << "method " << k;
	foreach(i, cc) {
		cout << "\t" << taxnames[strain_taxon[i]];
	}
	cout << endl;
	foreach(j, (size_t)dc) {
		cout << string(dv[j * 2]) << "\t";
		foreach(i, cc) {
			cout << scoretable[k][j][i] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void scoretablerecorder(string* prefix, size_t k) {
	ofstream outfile;
	ostringstream ss; ss << k;
	string outfn = string(*prefix).append(".st").append(ss.str());
	outfile.open(outfn.c_str());

	if(outfile.good()) {
		outfile << "method " << k;
		foreach(i, cc) {
			outfile << "\t" << strain_taxon[i];
		}
		outfile << endl;
		foreach(j, (size_t)dc) {
			outfile << string(dv[j * 2]) << "\t";
			foreach(i, cc) {
				outfile << scoretable[k][j][i] << "\t";
			}
			outfile << endl;
		}
		outfile << endl;
		outfile.close();
	} else {
		cerr << "Couldn't open " << outfn << " for recording!" << endl;
		exit(152);
	}
}

// Perform the running of a data analysis pass, calling classify() and score() with appropriate threading options.
// If -disk is specified, this will cause classes to be loaded.

void core_run() {
	string snm = string("scoring pass ").append(my_to_string(RUNLEVEL));
	string cnm = string("classifying pass ").append(my_to_string(RUNLEVEL));

	if(resume_mode) {

		if(num_threads == 0) {
			//prepclock(dc, "reads", cnm);
			classify(0, dc);
		} else {
			//prepclock(dc, "reads", cnm);
			classifythreads(num_threads);
		}
	} else {
		// string sn = string("classifying pass ").append(my_to_string(RUNLEVEL));

		current_batch = 0;

		if(disk_mode) {
			//prepclock(dc * cc, "calcs", snm);

			size_t batch_start = 0;
			size_t batch_stop = dc;

			do {
				if(batch_size > 0) {
					batch_start = current_batch * batch_size;
					batch_stop = (current_batch + 1) * batch_size;
					if(batch_stop > (size_t)dc) batch_stop = dc;
					if(batch_start >= (size_t)dc) break;

					gen_data_tables(batch_start, batch_stop);

				}

				if(num_threads == 0) {
					foreach(i, cc) {
						// if(shortlist_histo[i] == 0 && RUNLEVEL > 1) continue; // completely skip classes that appear irrelevant to the dataset

						classes[i] = new Class(class_filenames->at(i));
						score(batch_start, batch_stop);

						delete classes[i];
						classes[i] = NULL;
					}
				} else {
					scorethreads(num_threads); // batch range is passed with global variables
				}

				if(batch_size > 0) {
					clear_data_tables(batch_start, batch_stop);
					current_batch++;
				}
			} while (batch_size > 0);

			//prepclock(dc, "reads", cnm);
			if(num_threads == 0) {
				classify(0, dc);
			} else {
				classifythreads(num_threads);
			}
		} else if(num_threads == 0) {
			//prepclock(dc, "reads", snm);

			size_t batch_start = 0;
			size_t batch_stop = dc;

			do {
				if(batch_size > 0) {
					batch_start = current_batch * batch_size;
					batch_stop = (current_batch + 1) * batch_size;
					if(batch_stop > (size_t)dc) batch_stop = dc;
					if(batch_start >= (size_t)dc) break;

					gen_data_tables(batch_start, batch_stop);

				}
				score(batch_start, batch_stop);

				if(batch_size > 0) {
					clear_data_tables(batch_start, batch_stop);
					current_batch++;
				}
			} while (batch_size > 0);

			//prepclock(dc, "reads", cnm);
			classify(0, dc);
		} else {
			//prepclock(dc, "reads", snm);

			size_t batch_start = 0;
			size_t batch_stop = dc;

			do {
				if(batch_size > 0) {
					batch_start = current_batch * batch_size;
					batch_stop = (current_batch + 1) * batch_size;
					if(batch_stop > (size_t)dc) batch_stop = dc;
					if(batch_start >= (size_t)dc) break;

					gen_data_tables(batch_start, batch_stop);

				}

				scorethreads(num_threads);

				if(batch_size > 0) {
					clear_data_tables(batch_start, batch_stop);
					current_batch++;
				}
			} while (batch_size > 0);

			//prepclock(dc, "reads", cnm);
			classifythreads(num_threads);
		}
	}
}

// Reload temporary .st* files as generated by -record

extern string METHOD_NAMES[MAX_VOTE];

void reload_scoretables(string* prefix, size_t k) {

	string md = string(METHOD_NAMES[k]).append(" reload (");
	md.append(my_to_string(k+1)).append("/").append(my_to_string(MAX_VOTE)).append(")");

	time_t clock_start_time = time(NULL);
	string clock_unit = "scores/class";
	size_t prog_status = 0;

	//prepclock(dc, "scores/class", md);

	std::map<string, int> readindex;
	foreach(j, dc) {
		readindex[string(dv[j * 2])] = j;
	}

	std::map<int, int> reverse_taxid;
	std::map<int, vector<int> > reverse_taxid2;
	foreach(i, cc) {
		if(reverse_taxid.count(strain_taxon[i]) != 0) {
			if(reverse_taxid2.count(strain_taxon[i]) == 0) {
				reverse_taxid2[strain_taxon[i]].push_back(reverse_taxid[strain_taxon[i]]);
			}

			reverse_taxid2[strain_taxon[i]].push_back(i);
			reverse_taxid[strain_taxon[i]] = -1; // indicates extended data

		} else {
			reverse_taxid[strain_taxon[i]] = i;
		}
	}

	ifstream f;
	ostringstream ss; ss << k;
	string infn = string(*prefix).append(".st").append(ss.str());
	f.open(infn.c_str(), ios::in);

	std::map<int,int> colmap;
	std::map<int,int> colmap2;

	cerr << "[reload] Starting reload." << endl;

	if(f.is_open()) {
		string ibuff;
		if(f.good()) {
			getline(f, ibuff);
			//cout << " --- " << ibuff << endl;
			ibuff.append("\t");
			size_t x = 0, x2 = 0;
			size_t col = 0;
			size_t il = ibuff.length();
			while(x < il) {
				x++;
				if(ibuff[x] == '\t') {
					if(x2 != 0) {
						string val = ibuff.substr(x2, x - x2);
						size_t vali = atoi(val.c_str());

						if(reverse_taxid.count(vali) == 1) {
							colmap[col] = reverse_taxid[vali];
			//				cout << "Report: column " << col << " (" << vali << ")" << " represents class " << reverse_taxid[vali] << " (" << strain_taxon[reverse_taxid[vali]] << ")" << endl;

							if(colmap[col] == -1) {
								colmap2[col] = vali;
								//cout << "Report: colmap2 " << col << " points to reverse_taxid2 vector with " << reverse_taxid2[vali].size() << " members." << endl;
							}

						} else {
							cerr << "[reload] Ignoring scores for missing class (strain taxon id " << vali << ")" << endl;
						}
					}

					x2 = x + 1;
					col++;
				}
			}
		} else {
			cerr << "[reload] Could not read header row from resume file " << infn << endl;
			exit(21);
		}

		size_t r = 1;
		while(f.good()) {
			getline(f, ibuff);
			//cout << " --- " << ibuff << endl;
			ibuff.append("\t");
			size_t x = 0, x2 = 0;
			size_t col = 0;
			size_t il = ibuff.length();
			size_t j = 0;
			while(x < il) {
				x++;
				if(ibuff[x] == '\t') {
					if(colmap.count(col) == 1 || col == 0) { // skip columns we have no class for
						string val = ibuff.substr(x2, x - x2);

						if(col > 0) {
							number vali = atof(val.c_str());
							int i = colmap[col];

							if(colmap.count(col) == 0) {
								cerr << "[reload] Assertion failed; header row truncated." << endl;
								exit(22);
							}

							if(i == -1) {

								int taxid = colmap2[(int)col];

								std::vector<int> *hi = &(reverse_taxid2[taxid]);
								foreach(ii, hi->size()) {
									size_t pi = hi->at(ii);
									if(scoretable[k][j][pi] == -1) {
										i = pi;
										//cout << "i divined to be " << i << endl;
										break;
									} else {
										//cout << "class " << pi << " already filled..." << endl;
									}
								}

								if(r == 1 && i == -1) {
									cerr << "[reload] Input file has more classes with taxid = " << taxid << " than input class spec." << endl;
									cerr << "         (Found in column " << col << ")" << endl;
									cerr << "         Continuation is probably invalid." << endl;
									cerr << "         DO NOT reuse taxids if changing class spec during a continuation." << endl;
									cerr << "         Scores for classes sharing taxid may be switched." << endl;
								}
							}

							scoretable[k][j][i] = vali;

							if((k == AA_NB || k == AA_NN || k == AA_MM || k == AA_SB || k == AA_CD) && vali > 0.0)
								has_protein[j] = true;

							//cout << "Report: i,j,k,score = " << i << "," << j << "," << k << "," << val << "=" << vali << "." << endl;

						} else {
							if(val.length() == 0) {
								// no label, probably a duplicated header
								//cout << "Report: look at how quiet I am." << endl;
								break;
							} else if(readindex.count(val) == 0) {
								if(val.find("method ") == string::npos)
									cerr << "[reload] Ignoring scores for unexpected read: \"" << val << "\"" << endl;
								break;
							} else {

								prog_status++;

								j = readindex[val];
								has_protein[j] = false;
								//cout << "Report: row " << r << " (" << val << ")" << " represents read " << j << " (" << dv[j * 2] << ")" << endl;
							}

						}
					} else if(colmap.count(col) == 0) {
						// cout << "Report: skipped column " << col << endl;
					}

					x2 = x + 1;
					col++;
				}
			} // end line reader

			r++;
			runClock(prog_status, (size_t)dc, clock_start_time, clock_unit, md);
		}

		// cerr << "[reload] Read " << r << " line(s) from " << infn << endl;

		f.close();
	} else {
		cerr << "[reload] Could not open resume file " << infn << endl;
		exit(23);
	}
	cerr << endl;

	cerr << "[reload] Done reload." << endl;
}

int main(int argc, char** argv) {

 	cerr << "generative inference of sequence taxonomy (gist)" << endl;
 	cerr << "version " << GIST_VERSION << endl;
	cerr << "built " << BUILD_TIME << endl;
	cerr << "2012-2017 rhetorica@cs.toronto.edu" << endl;
	cerr << endl;

	output_mode = OUTPUT_READABLE;

	string fn_16s = "", fn_ap = "";

	class_filenames = new vector<string>();

	string *class_dir = NULL;
	string fgs_profile = "";

	bool showhelp = false, showhmms = false;

	string* confpath = new string(getexepath());
	confpath->append("/gist.conf");

	try {
		loadConf(*confpath);
	} catch(string* e) {
		cerr << *e << endl;
		return 111;
	}

	delete confpath;

	// argument loading
	char** args;
	argc = override_args(args, argc, argv);

	for(int argi = 1; argi < argc; argi++) {
		char* argj;

		if(* (args[argi]) == '-' && * (args[argi] + 1) == '-')
			argj = args[argi] + 1;
		else
			argj = args[argi];

		if(string("-?") == string(argj)) {
			showhelp = true;
		} else if(string("-help") == string(argj)) {
			showhelp = true;
		} else if(string("-sqrt") == string(argj)) {
			no_sqrt = false;
		} else if(string("-class") == string(argj)) {
			argi++;
			class_filenames->push_back(args[argi]);
		} else if(string("-classes") == string(argj)) {
			argi++;
			class_dir = new string(args[argi]);
		} else if(string("-n") == string(argj)) {
			argi++;
			n_prefix_length = atoi(args[argi]);
		} else if(string("-p") == string(argj)) {
			argi++;
			p_prefix_length = atoi(args[argi]);
		} else if(string("-t") == string(argj)) {
			argi++;
			num_threads = atoi(args[argi]);
		} else if(string("-b") == string(argj) ||
				string("-batch") == string(argj)) {
			argi++;
			batch_size = atoi(args[argi]);
		} else if(string("-bwa") == string(argj)) {
			argi++;
			bwa_weight = atof(args[argi]);
		} else if(string("-snw") == string(argj)) {
			argi++;
			sn_weight = atof(args[argi]);
		} else if(string("-spw") == string(argj)) {
			argi++;
			sp_weight = atof(args[argi]);
		} else if(string("-bnw") == string(argj)) {
			argi++;
			bn_weight = atof(args[argi]);
		} else if(string("-bpw") == string(argj)) {
			argi++;
			bp_weight = atof(args[argi]);
		} else if(string("-nnw") == string(argj)) {
			argi++;
			nn_weight = atof(args[argi]);
		} else if(string("-npw") == string(argj)) {
			argi++;
			np_weight = atof(args[argi]);
		} else if(string("-cnw") == string(argj)) {
			argi++;
			cn_weight = atof(args[argi]);
		} else if(string("-cpw") == string(argj)) {
			argi++;
			cp_weight = atof(args[argi]);
		} else if(string("-mnw") == string(argj)) {
			argi++;
			mn_weight = atof(args[argi]);
		} else if(string("-mpw") == string(argj)) {
			argi++;
			mp_weight = atof(args[argi]);
		} else if(string("-sn") == string(argj)) {
			argi++;
			sn_prefix_length = atoi(args[argi]);
		} else if(string("-sp") == string(argj)) {
			argi++;
			sp_prefix_length = atoi(args[argi]);
		} else if(string("-bn") == string(argj)) {
			argi++;
			bn_prefix_length = atoi(args[argi]);
		} else if(string("-bp") == string(argj)) {
			argi++;
			bp_prefix_length = atoi(args[argi]);
		} else if(string("-nn") == string(argj)) {
			argi++;
			nn_prefix_length = atoi(args[argi]);
		} else if(string("-np") == string(argj)) {
			argi++;
			np_prefix_length = atoi(args[argi]);
		} else if(string("-cn") == string(argj)) {
			argi++;
			cn_prefix_length = atoi(args[argi]);
		} else if(string("-cp") == string(argj)) {
			argi++;
			cp_prefix_length = atoi(args[argi]);
		} else if(string("-mn") == string(argj)) {
			argi++;
			mn_prefix_length = atoi(args[argi]);
		} else if(string("-mp") == string(argj)) {
			argi++;
			mp_prefix_length = atoi(args[argi]);
		} else if(string("-dump") == string(argj)) {
			debugdump = true;
		} else if(string("-dgmm") == string(argj)) {
			debuggmm = true;
		} else if(string("-noise") == string(argj)) {
			argi++;
			noise_rate = atof(args[argi]);
		} else if(string("-data") == string(argj)) {
			argi++;
			dataname = new string(args[argi]);
		} else if(string("-conf") == string(argj)) {
			argi++;
			loadConf(string(args[argi]));
		} else if(string("-cc") == string(argj)) {
			argi++;
			autocross_generate = true;
			autocross_weights_fn = string(args[argi]);
		} else if(string("-cross") == string(argj)) {
			argi++;
			autocross_use = true;
			autocross_weights_fn = string(args[argi]);
		} else if(string("-save") == string(argj)) {
			outclasses = true;
		} else if(string("-ts") == string(argj)) {
			trust_source = true;
		} else if(string("-tcs") == string(argj)) {
			trust_class_source = false;
		} else if(string("-hmms") == string(argj)) {
			showhmms = true;
		} else if(string("-fgs") == string(argj)) {
			argi++;
			fgs_profile = string(args[argi]);
		} else if(string("-flen") == string(argj)) {
			argi++;
			FGS_THRESHOLD = atof(args[argi]);
		} else if(string("-ct") == string(argj)) {
			use_class_threads = true;
		} else if(string("-r") == string(argj)) {
			argi++;
			class_read_length = atoi(args[argi]);
		} else if(string("-perf") == string(argj)) {
			performance_assessment = true;
		} else if(string("-o") == string(argj)) {
			argi++;
			char o = args[argi][0];
			switch(o) {
			case 's':
				output_mode = OUTPUT_SCORES;
				break;
			case 'S':
				output_mode = OUTPUT_SCORETABLES;
				break;
			case 'r':
				output_mode = OUTPUT_READABLE;
				break;
			case 'c':
				output_mode = OUTPUT_CSV;
				break;
			default:
				cerr << "[head] Unknown output mode: " << o << endl;
				exit(40);
			}
		} else if(string("-os") == string(argj)) {
			output_mode = OUTPUT_SCORES;
		} else if(string("-oS") == string(argj)) {
			output_mode = OUTPUT_SCORETABLES;
		} else if(string("-or") == string(argj)) {
			output_mode = OUTPUT_READABLE;
		} else if(string("-oc") == string(argj)) {
			output_mode = OUTPUT_CSV;
		} else if(string("-fullread") == string(argj)) {
			class_read_length = 0;
		} else if(string("-s1q") == string(argj)) {
			argi++;
			s1_include_quota = atoi(args[argi]);
		} else if(string("-s1t") == string(argj)) {
			argi++;
			s1_include_threshold = atof(args[argi]);
		} else if(string("-s2q") == string(argj)) {
			argi++;
			s2_include_quota = atoi(args[argi]);
		} else if(string("-s2t") == string(argj)) {
			argi++;
			s2_include_threshold = atof(args[argi]);
		} else if(string("-16s") == string(argj)) {
			argi++;
			fn_16s = string(args[argi]);
		} else if(string("-ap") == string(argj)) {
			argi++;
			fn_ap = string(args[argi]);
		} else if(string("-16sw") == string(argj)) {
			argi++;
			p16s_weight = atof(args[argi]);
		} else if(string("-apw") == string(argj)) {
			argi++;
			ap_weight = atof(args[argi]);
		} else if(string("-crank") == string(argj)) {
			argi++;
			permissable_taxrank = atoi(args[argi]);
		} else if(string("-resume") == string(argj)) {
			resume_mode = true;
		} else if(string("-record") == string(argj)) {
			output_mode = OUTPUT_RECORD;
		} else if(string("-disk") == string(argj)) {
			disk_mode = true;
		} else if(string("-nc") == string(argj)) {
			no_clock_display = true;
		} else {
			cerr << "[head] Unrecognized parameter: " << args[argi] << endl << endl;
			cerr << "       Consult " << args[0] << " -help for more information." << endl;
		}
	}

	if(autocross_use && autocross_generate) {
		cerr << "[head] Can't both generate and use autocross weights in a single session." << endl;
		exit(41);
	}

	if(performance_assessment && autocross_generate) {
		cerr << "[head] WARNING: -cc automatically performs performance assessment. -perf ignored." << endl;
		performance_assessment = false;
	}

	if(num_threads < 0) {
		cerr << "[head] Can't use " << num_threads << " as a number of threads." << endl;
		exit(42);
	}

	// show HMMs
	if(showhmms) {
		cout << "built-in hidden markov profiles:" << endl;
		cout << "    complete        0.0% read error rate (for reference sequences)" << endl;
		cout << "    454_5           454 pyrosequencing with 0.5% read error rate" << endl;
		cout << "    454_10          454 pyrosequencing with 1.0% read error rate" << endl;
		cout << "    454_30          454 pyrosequencing with 3.0% read error rate" << endl;
		cout << "    illumina_1      454 pyrosequencing with 0.1% read error rate" << endl;
		cout << "    illumina_5      454 pyrosequencing with 0.5% read error rate" << endl;
		cout << "    illumina_10     454 pyrosequencing with 1.0% read error rate" << endl;
		cout << "    sanger_5        454 pyrosequencing with 0.5% read error rate" << endl;
		cout << "    sanger_10       454 pyrosequencing with 1.0% read error rate" << endl;
		cout << endl;
		cout << "    (you can add new HMM profiles to " << getexepath() << "/fgs_profiles/" << endl;
		cout << "    but this list is hard-coded and will not show them.)" << endl;
		cout << endl;

		return 0;
	}

	// show help
	if(showhelp) {
		cout << "usage:" << endl;
		cout << "    " << argv[0] << " -data <file> -classes <dir>                  [<options> ... ]" << endl;
		cout << "    " << argv[0] << " -save        -classes <dir>                  [<options> ... ]" << endl;
		cout << "    " << argv[0] << " -data <file> -class <file> -class <file> ... [<options> ... ]" << endl;
		cout << "    " << argv[0] << " -save        -class <file> -class <file> ... [<options> ... ]" << endl;
		cout << endl;
		cout << "information:" << endl;
		cout << "    -?              show help" << endl;
		cout << "    -help           also show help" << endl;
		cout << "    -hmms           list available HMM profiles and quit" << endl;
//		cout << endl;
//		cout << "    -dump           dump extra information during output" << endl;
//		cout << "    -dgmm           dump GMM development information" << endl;
		cout << endl;
		cout << "output format:" << endl;
		cout << "    -o s            output raw final-pass scores (for manual adjustment)" << endl;
		cout << "    -o S            output raw score tables (for all methods)" << endl;
		cout << "    -o c            output CSV" << endl;
		cout << "    -o r            output readable taxon names (default)" << endl;
		cout << "    -perf           report performance, assuming autocross labelled data" << endl;
		cout << endl;
		cout << "data sources:" << endl;
		cout << "    -class <.ffn>   generate a class from a nucleotide FASTA file" << endl;
		cout << "    -class <.gist>  re-load a class from a pre-calculated profile" << endl;
		cout << "    -save           save class profiles to GIST files after calculation" << endl;
		cout << "    -classes <dir>  bulk load FASTA and/or GIST files from a directory" << endl;
		cout << endl;
		cout << "    -ap <.tsv>      load a priori log weights from 'strain, score' TSV file" << endl;
		cout << "    -16s <.tsv>     load 16S count weights from 'strain, count' TSV file" << endl;
		cout << endl;
		cout << "    -data <.fa>     list of reads to classify" << endl;
		cout << endl;
		cout << "cross validation:" << endl;
		cout << "    -cc <.tsv>      generate autocross weightings and store in TSV file" << endl;
		cout << "                    (see the manual for how to prepare autocross data)" << endl;
		cout << "    -cross <.tsv>   load autocross weighting results from TSV file" << endl;
		cout << "    -crank <int>    desired accuracy [default: 2, genus level]" << endl;
		cout << "                    (see manual for details)" << endl;
		cout << endl;
		cout << "delayed/distributed processing:" << endl;
		cout << "    -resume         skip score generation; use old scores from temp files" << endl;
		cout << "    -record         do not classify; just write scores to temp files" << endl;
		cout << endl;
		cout << "running behavior:" << endl;
		cout << "    -t <int>        how many threads to use for classification [current: " << num_threads << "]" << endl;
		cout << "    -ct             enable (maximum) threads for class construction" << endl;
		cout << "    -nc             disable display of progress meter" << endl;
		cout << endl;
		cout << "    -disk           load classes one at a time (use with many classes)" << endl;
		cout << "    -b <int>        how many data tables to load at once [current: " << batch_size << "]" << endl;
		cout << "                    (use with large tuple lengths and large data sets)" << endl;
		cout << "                    (0 = all at once)" << endl;
		cout << endl;
		cout << "parameters:" << endl;
/*		cout << "    -noise <frac>   randomly add bad nucleotide calls everywhere [current: " << noise_rate << "]" << endl;
		cout << "    -sqrt           square root data normalization (not recommended)" << endl;
		cout << endl;*/
		cout << "    -bwa <ratio>    re-weighting of nucleotide BWA [current: " << bwa_weight << "]" << endl;
		cout << endl;
		cout << "    -conf <.conf>   load base weighting config file [default: gist.conf]" << endl;
		cout << "    -bnw <ratio>    re-weighting of nucleotide NB [current: " << bn_weight << "]" << endl;
		cout << "    -bpw <ratio>    re-weighting of protein NB [current: " << bp_weight << "]" << endl;
		cout << "    -snw <ratio>    re-weighting of sparse nucleotide NB [current: " << sn_weight << "]" << endl;
		cout << "    -spw <ratio>    re-weighting of sparse protein NB [current: " << sp_weight << "]" << endl;
		cout << "    -nnw <ratio>    re-weighting of nucleotide 1NN [current: " << nn_weight << "]" << endl;
		cout << "    -npw <ratio>    re-weighting of protein 1NN [current: " << np_weight << "]" << endl;
		cout << "    -cnw <ratio>    re-weighting of nucleotide codelta [current: " << cn_weight << "]" << endl;
		cout << "    -cpw <ratio>    re-weighting of protein codelta [current: " << cp_weight << "]" << endl;
		cout << "    -mnw <ratio>    re-weighting of nucleotide GMM [current: " << mn_weight << "]" << endl;
		cout << "    -mpw <ratio>    re-weighting of protein GMM [current: " << mp_weight << "]" << endl;
		cout << "    -mnk <int>      clusters for nucleotide GMM [current: " << NT_MM_K << "]" << endl;
		cout << "    -mpk <int>      clusters for protein GMM [current: " << AA_MM_K << "]" << endl;
		cout << "    -apw <ratio>    re-weighting of a priori weights [current: " << ap_weight << "]" << endl;
		cout << "    -16sw <ratio>   re-weighting of 16S weights [current: " << p16s_weight << "]" << endl;
		cout << endl;
		cout << "                    (set a weight to 0 to disable a method)" << endl;
		cout << endl;
		cout << "    -n <int>        prefix length for nucleotide K-mers** [current: " << n_prefix_length << "]" << endl;
		cout << "                    (1 will be added to get final K-mer length)" << endl;
		cout << "    -p <int>        prefix length for protein K-mers** [current: " << p_prefix_length << "]" << endl;
		cout << "                    (1 will be added to get final K-mer length)" << endl;
		cout << "    -bn <int>       real prefix length for nucleotide NB [current: " << bn_prefix_length << "]" << endl;
		cout << "    -bp <int>       real prefix length for protein NB [current: " << bp_prefix_length << "]" << endl;
		cout << "    -sn <int>       real prefix length for sparse nucleotide NB [current: " << sn_prefix_length << "]" << endl;
		cout << "    -sp <int>       real prefix length for sparse protein NB [current: " << sp_prefix_length << "]" << endl;
		cout << "    -nn <int>       real prefix length for nucleotide 1NN [current: " << nn_prefix_length << "]" << endl;
		cout << "    -np <int>       real prefix length for protein 1NN [current: " << np_prefix_length << "]" << endl;
		cout << "    -cn <int>       real prefix length for nucleotide codelta [current: " << cn_prefix_length << "]" << endl;
		cout << "    -cp <int>       real prefix length for protein codelta [current: " << cp_prefix_length << "]" << endl;
		cout << "    -mn <int>       real prefix length for nucleotide GMM [current: " << mn_prefix_length << "]" << endl;
		cout << "    -mp <int>       real prefix length for protein GMM [current: " << mp_prefix_length << "]" << endl;
/*		cout << " *  -r <int>        expected read length [default: 100 nt] [UNIMPLEMENTED]" << endl;
		cout << "    -fullread       assume read length is infinite, don't resample [default]" << endl;*/
		cout << endl;
		cout << "    -s1q <int>      pass 1 should return at least <int> hits [current: " << s1_include_quota << "]" << endl;
		cout << "    -s1t <frac>     pass 1 should return all hits with log-p above <frac>" << endl;
		cout << "                    of best hit [current: " << s1_include_threshold << "]" << endl;
		cout << "    -s2q <int>      final pass should return at least <int> hits [current: " << s2_include_quota << "]" << endl;
		cout << "    -s2t <frac>     final pass should return all hits with log-p above <frac>" << endl;
		cout << "                    of best hit [current: " << s1_include_threshold << "]" << endl;
		cout << endl;
		cout << "                    (expected read lengths are used in class construction)" << endl;
		cout << endl;
		cout << "    -ts             assume data reads are in the correct frame" << endl;
		cout << "    -tcs            assume class reads are NOT in the correct frame" << endl;
		cout << "    -fgs <hmm>      HMM profile to use for FragGeneScan ORF detection" << endl;
		cout << "    -flen <int>     minimum read size for ORF detection [default: " << FGS_THRESHOLD << "]" << endl;
		cout << endl;
//		cout << " * = unimplemented" << endl;
		cout << " ** = does not affect sparse methods" << endl;
		cout << endl;

		return 0;
	}

	// ENABLE_METHODs management:
	if(p_prefix_length > 0) {
		ENABLE_METHOD[AA_NB] = (bp_weight > 0);
		ENABLE_METHOD[AA_NN] = (np_weight > 0);
		ENABLE_METHOD[AA_CD] = (cp_weight > 0);
		ENABLE_METHOD[AA_MM] = (mp_weight > 0);
		ENABLE_METHOD[AA_SB] = (sp_weight > 0);
	} else {
		ENABLE_METHOD[AA_NB] = false;
		ENABLE_METHOD[AA_CD] = false;
		ENABLE_METHOD[AA_NN] = false;
		ENABLE_METHOD[AA_MM] = false;
		ENABLE_METHOD[AA_SB] = false;
	}

	ENABLE_METHOD[NT_BWA] = (bwa_weight > 0);

	if(n_prefix_length > 0) {
		ENABLE_METHOD[NT_NB] = (bn_weight > 0);
		ENABLE_METHOD[NT_NN] = (nn_weight > 0);
		ENABLE_METHOD[NT_CD] = (cn_weight > 0);
		ENABLE_METHOD[NT_MM] = (mn_weight > 0);
		ENABLE_METHOD[NT_SB] = (sn_weight > 0);
	} else {
		ENABLE_METHOD[NT_NB] = false;
		ENABLE_METHOD[NT_CD] = false;
		ENABLE_METHOD[NT_NN] = false;
		ENABLE_METHOD[NT_MM] = false;
		ENABLE_METHOD[NT_SB] = false;
	}

	foreach(i, MAX_VOTE) {
		if(i != PRIOR_16S && i != PRIOR_AP)
			UNWEIGHTED_METHOD[i] = false;
		else
			UNWEIGHTED_METHOD[i] = true;

		if(method_weights[0][i] == 0 && method_weights[1][i] == 0) ENABLE_METHOD[i] = false;
	}

	ENABLE_METHOD[NT_NB_RC] = ENABLE_METHOD[NT_NB];
	ENABLE_METHOD[NT_NN_RC] = ENABLE_METHOD[NT_NN];
	ENABLE_METHOD[NT_CD_RC] = ENABLE_METHOD[NT_CD];
	ENABLE_METHOD[NT_MM_RC] = ENABLE_METHOD[NT_MM];
	ENABLE_METHOD[NT_SB_RC] = ENABLE_METHOD[NT_SB];

	foreach(i, MAX_VOTE) USE_METHOD[i] = ENABLE_METHOD[i];

	// begin max prefix setup - sparse methods are excluded since they handle differently

	max_n_prefix_length = 0;
	if(ENABLE_METHOD[NT_NB] && bn_prefix_length > max_n_prefix_length) max_n_prefix_length = bn_prefix_length;
	if(ENABLE_METHOD[NT_CD] && cn_prefix_length > max_n_prefix_length) max_n_prefix_length = cn_prefix_length;
	if(ENABLE_METHOD[NT_NN] && nn_prefix_length > max_n_prefix_length) max_n_prefix_length = nn_prefix_length;
	if(ENABLE_METHOD[NT_MM] && mn_prefix_length > max_n_prefix_length) max_n_prefix_length = mn_prefix_length;

	max_p_prefix_length = 0;
	if(ENABLE_METHOD[AA_NB] && bp_prefix_length > max_p_prefix_length) max_p_prefix_length = bp_prefix_length;
	if(ENABLE_METHOD[AA_CD] && cp_prefix_length > max_p_prefix_length) max_p_prefix_length = cp_prefix_length;
	if(ENABLE_METHOD[AA_NN] && np_prefix_length > max_p_prefix_length) max_p_prefix_length = np_prefix_length;
	if(ENABLE_METHOD[AA_MM] && mp_prefix_length > max_p_prefix_length) max_p_prefix_length = mp_prefix_length;

	n_prefix_enabled = new bool[max_n_prefix_length + 1];
	p_prefix_enabled = new bool[max_p_prefix_length + 1];

	foreach(i, max_n_prefix_length + 1) n_prefix_enabled[i] = false;
	foreach(i, max_p_prefix_length + 1) p_prefix_enabled[i] = false;

	if(ENABLE_METHOD[NT_NB]) n_prefix_enabled[bn_prefix_length] = true;
	if(ENABLE_METHOD[NT_CD]) n_prefix_enabled[cn_prefix_length] = true;
	if(ENABLE_METHOD[NT_NN]) n_prefix_enabled[nn_prefix_length] = true;
	if(ENABLE_METHOD[NT_MM]) n_prefix_enabled[mn_prefix_length] = true;

	if(ENABLE_METHOD[AA_NB]) p_prefix_enabled[bp_prefix_length] = true;
	if(ENABLE_METHOD[AA_CD]) p_prefix_enabled[cp_prefix_length] = true;
	if(ENABLE_METHOD[AA_NN]) p_prefix_enabled[np_prefix_length] = true;
	if(ENABLE_METHOD[AA_MM]) p_prefix_enabled[mp_prefix_length] = true;

	// end max prefix setup

	if(!trust_source && p_prefix_length > 0) {
		string mepath = getexepath();
		mepath.append("/fgs_profiles/");
		if(fgs_profile.length() == 0) {
			fgs_profile = "complete";

			if(!resume_mode)
				cerr << "[head] WARNING: No HMM profile was specified; assuming 0% read error rate." << endl
				 	 << "       If you are working with in-frame genes and not raw reads, use -ts to" << endl
				 	 << "       switch to full-sequence translation instead of ORF prediction." << endl
				 	 << "       Use -fgs <profile> to specify a profile and -hmms to see a list." << endl;
		}
		FragGeneScan_setup(mepath.c_str(), fgs_profile.c_str(), 0);
	}

	// get classes:

	if(class_dir != NULL) {
		get_paths(*class_dir, class_filenames);
		delete class_dir;
		class_dir = NULL;
	}

	cc = class_filenames->size();

	classes = new Class*[cc];

	if(cc == 0) {
		cerr << "[head] No classes specified. Terminating." << endl
		     << "       See gist -help for information on specifying classes." << endl;
		return 1;
	}

	cerr << endl;

	// load class priors:

	if(p16s_weight > 0 && fn_16s != "") {
		try {
			load16S(fn_16s);
		} catch(string* e) {
			cerr << e;
			return 112;
		}
		ENABLE_METHOD[PRIOR_16S] = true;
	} else {
		ENABLE_METHOD[PRIOR_16S] = false;
	}

	if(ap_weight > 0 && fn_ap != "") {
		try {
			loadAPriori(fn_ap);
		} catch(string* e) {
			cerr << e;
			return 113;
		}

		ENABLE_METHOD[PRIOR_AP] = true;
	} else {
		ENABLE_METHOD[PRIOR_AP] = false;
	}

	// load classes:
	// This step is skipped if -disk is enabled, in which case the program
	// will load classes as needed during the scoring process.

	if(!disk_mode) {
		if(use_class_threads)
			classthreads(class_filenames, classes);
		else
			classnothreads(class_filenames, classes);

		if(debugdump) {
			cout << endl << "-------------------------------------" << endl;
			cout << "Class specifications";
			cout << endl << "-------------------------------------" << endl;
			foreach(i, cc) {
				classes[i]->info();
			}
		}

		cerr << "[head] " << cc << " classes loaded in simultaneous mode." << endl;

		if(outclasses) cerr << "[head] Class profiles were saved." << endl;
	} else {
		foreach(i, cc) classes[i] = NULL;
		cerr << "[head] " << cc << " classes found; loading deferred for disk mode." << endl;
	}

	// load data:

	if(dataname == NULL) {
		cerr << "[head] No data file provided." << endl;
		return 2;
	}

	if(ENABLE_METHOD[NT_BWA] && !resume_mode) {
		cerr << "[head] Starting BWA runs." << endl;
		runBWA(class_filenames, num_threads, *dataname);
		cerr << "[head] Finished BWA runs." << endl;
	}

	cerr << endl;

	dc = readFasta(*dataname, dv, resume_mode);

	if(dc == -1) {
		cerr << "[head] No data specified." << endl;
		return 2;
	} else if(dc == -2) {
		cerr << "[head] Data file not in FASTA format." << endl;
		return 3;
	} else {
		cerr << "[head] Loaded data file: " << *dataname << endl;
	}

	// create score tables:

	has_protein = new bool[dc];

	foreach(i, MAX_VOTE) if(ENABLE_METHOD[i]) {
		scoretable[i] = new number*[dc];
		foreach(j, (size_t)dc) {
			scoretable[i][j] = new number[cc];
			foreach(k, cc) {
				scoretable[i][j][k] = -numeric_limits<number>::min();
			}
		}
	} else {
		scoretable[i] = NULL;
	}

	foreach(i, MAX_RUNLEVEL) {
		overall_scoretable[i] = new number*[dc];
		foreach(j, (size_t)dc) {
			overall_scoretable[i][j] = new number[cc];
			foreach(k, cc) {
				overall_scoretable[i][j][k] = 0;
			}
		}
	}

	if(!resume_mode) {
		scoretable_NT_MM = new number**[NT_MM_K];
		scoretable_NT_MM_RC = new number**[NT_MM_K];
		foreach(i, NT_MM_K) {
			scoretable_NT_MM[i] = new number*[dc];
			scoretable_NT_MM_RC[i] = new number*[dc];
			foreach(j, (size_t)dc) {
				scoretable_NT_MM[i][j] = new number[cc];
				scoretable_NT_MM_RC[i][j] = new number[cc];
			}
		}

		scoretable_AA_MM = new number**[AA_MM_K];
		foreach(i, AA_MM_K) {
			scoretable_AA_MM[i] = new number*[dc];
			foreach(j, (size_t)dc) {
				scoretable_AA_MM[i][j] = new number[cc];
			}
		}
	}

	cerr << "[head] Score tables generated." << endl;
	cerr << "[head] Loading lineage data." << endl;

	// load lineage information:

	taxa = new size_t*[cc];
	strain_taxon = new size_t[cc];
	foreach(i, cc) {
		taxa[i] = new size_t[MAX_RANK];
		string lin_name = string(class_filenames->at(i)).append(".lin");
		if(lin_name.substr(lin_name.length() - 9, 9) == string(".gist.lin")) {
			lin_name = lin_name.substr(0, lin_name.length() - 9).append(".lin");
		}

		if(exists(lin_name)) {
			loadLineage(lin_name, i);
		} else {
			bool warn = true;
			if(warn) cerr << "Couldn't find lineage file " << lin_name << "!" << endl;
			foreach(j, MAX_RANK) {
				taxa[i][j] = 0;
			}
		}
	}

	cerr << "[head] Lineage data loaded." << endl;

	if(autocross_generate) {
		autocross_load_labels_from_fasta();
		autocross_prep();
	} else if(autocross_use) {
		autocross_prep();
		autocross_load_weights(autocross_weights_fn);
	}

	// construct shortlist space:

	shortlist = new bool*[cc];
	foreach(i, cc) {
		shortlist[i] = new bool[dc];
		foreach(j, (size_t)dc) shortlist[i][j] = false;
	}

	// analysis:

	if(ENABLE_METHOD[NT_BWA] && !resume_mode) {
		cerr << "[head] Extracting BWA scores." << endl;
		//prepclock(cc * dc, "matches", "extracting aligner scores");
		parseSAMs(*dataname);
		cerr << "[head] Done extracting BWA scores." << endl;
	}

	taxscores = new std::map<int,number>[dc];

	if(resume_mode) {
		cerr << "[head] Reloading score tables." << endl;
		forclassifiers(k) {
			reload_scoretables(dataname, k);
		}
	} else if(batch_size == 0) {
		cerr << "[head] Generating data tables." << endl;
		gen_data_tables();
	}

	cerr << endl;

	if(autocross_generate || output_mode == OUTPUT_SCORETABLES || output_mode == OUTPUT_RECORD) {
		// full runlevel:
		RUNLEVEL = 1;

		if(autocross_generate)
			cerr << "[head] Beginning autocross runlevel." << endl;
		else
			cerr << "[head] Beginning full runlevel." << endl;

		foreach(j, MAX_VOTE) {
			USE_METHOD[j] = ENABLE_METHOD[j];
			ADD_METHOD[j] = USE_METHOD[j];
		}

		core_run();

		if(output_mode == OUTPUT_RECORD) {
			cerr << endl << "[head] Writing scoretables to temp files." << endl;
			forclassifiers(k) {
				scoretablerecorder(dataname, k);
			}
		} else if(output_mode == OUTPUT_SCORETABLES) {
			cerr << endl << "[head] Writing scoretables to stdout." << endl;
			forclassifiers(k) {
				scoretableprinter(k);
			}
		} else {
			cerr << endl << "[head] Beginning autocross learning." << endl;

			autocross_train();

			cerr << endl << "[head] Writing autocross weights to file." << endl;
			autocross_save_weights(autocross_weights_fn);
		}
	} else {
		// runlevel 1:

		RUNLEVEL = 1;

		cerr << "[head] Beginning runlevel 1." << endl;

		foreach(j, MAX_VOTE) {
			if(method_weights[RUNLEVEL-1][j] == 0)
				USE_METHOD[j] = false;
			else
				USE_METHOD[j] = ENABLE_METHOD[j];

			ADD_METHOD[j] = USE_METHOD[j];
		}

		core_run();

		// TODO: antihip here

		// write shortlists:
		cerr << endl << "[head] Beginning shortlist processing." << endl; cerr.flush();

		// generate shortlist:
		string pass_name = string("constructing shortlist for pass 1");
		do_shortlist(s1_include_quota, s1_include_threshold, FAMILY, s1_bloom_threshold, pass_name);

		cerr << endl << "[head] Done shortlist processing." << endl;

		// runlevel 2:

		RUNLEVEL = 2;

		cerr << "[head] Beginning runlevel 2." << endl;

		foreach(j, MAX_VOTE) {
			if(method_weights[RUNLEVEL-1][j] == 0)
				USE_METHOD[j] = false;
			else
				USE_METHOD[j] = ENABLE_METHOD[j];

			if(ADD_METHOD[j])
				ADD_METHOD[j] = false;
			else
				ADD_METHOD[j] = USE_METHOD[j];

			USE_METHOD[j] = ENABLE_METHOD[j]; // last round
		}

		core_run();

		// TODO: antihip here again

		if(output_mode == OUTPUT_SCORES) {
			scoreprinter();

			if(performance_assessment) {
				cerr << endl << "[head] Writing performance summary." << endl;
				autocross_load_labels_from_fasta();
				size_t hitcount = 0;
				foreach(j, dc) {
					int maxi = -1;
					number maxv = 0;
					foreach(i, cc) {
						if(maxi == -1 || overall_scoretable[RUNLEVEL-1][j][i] > maxv) {
							maxv = overall_scoretable[RUNLEVEL-1][j][i];
							maxi = i;
						}
					}
					if(readlabels[j][maxi] == 1.0) hitcount++;
				}

				cout << "Performance: " << hitcount << "/" << dc << " reads classified such that best hit is within tolerances on file " << string(*dataname) << endl;
			}
		} else {
			cerr << endl << "[head] Determining final shortlist." << endl;

			foreach(i, cc) foreach(j, (size_t)dc) shortlist[i][j] = false;

			// generate shortlist:
			pass_name = string("constructing shortlist for pass 2");
			do_shortlist(s2_include_quota, s2_include_threshold, GENUS, s2_bloom_threshold, pass_name);

			cerr << endl << "[head] Writing report." << endl;
			write_output();


			if(performance_assessment) {
				cerr << endl << "[head] Writing performance summary." << endl;
				autocross_load_labels_from_fasta();
				size_t hitcount = 0;
				foreach(j, dc) {
					foreach(i, cc) {
						if(shortlist[i][j] == true && readlabels[j][i] == 1.0) {
							hitcount++;
							break;
						}
					}
				}

				cout << "Performance: " << hitcount << "/" << dc << " reads classified so that some hit is within tolerances on file " << string(*dataname) << endl;
			}
		}
	}

	cerr << "[head] Run complete." << endl;

	// cleanup:

	if(args != argv) {
		foreach(i, (size_t)argc) {
			delete[] args[i];
		}

		delete[] args;
	}

	delete class_filenames;

	delete dataname;

	foreach(i, cc) {
		delete classes[i];
		delete[] taxa[i];

		delete[] shortlist[i];
	}

	delete[] classes;
	delete[] taxa;
	delete[] strain_taxon;

	delete[] shortlist;
	delete[] shortlist_histo;

	foreach(i, (size_t)dc) {
		delete[] dv[i];
		dv[i] = NULL;

		if(!resume_mode) {
			foreach(K, max_n_prefix_length + 1) {
				if(data_n[K] != NULL) {
					delete data_n[K][i]; data_n[K][i] = NULL;
					delete data_n_rc[K][i]; data_n_rc[K][i] = NULL;
				}
			}
			foreach(K, max_p_prefix_length + 1) {
				if(data_p[K] != NULL) {
					delete data_p[K][i]; data_p[K][i] = NULL;
				}
			}
		}

		if(ENABLE_METHOD[NT_SB]) { delete sparse_data_n[i]; sparse_data_n[i] = NULL; }
		if(ENABLE_METHOD[NT_SB_RC]) { delete sparse_data_n_rc[i]; sparse_data_n_rc[i] = NULL; }
		if(ENABLE_METHOD[AA_SB]) { delete sparse_data_p[i]; sparse_data_p[i] = NULL; }
	}

	if(ENABLE_METHOD[NT_SB]) { delete[] sparse_data_n; sparse_data_n = NULL; }
	if(ENABLE_METHOD[NT_SB_RC]) { delete[] sparse_data_n_rc; sparse_data_n_rc = NULL; }
	if(ENABLE_METHOD[AA_SB]) { delete[] sparse_data_p; sparse_data_p = NULL; }

	delete[] dv; dv = NULL;

	if(!resume_mode) {
		foreach(K, max_n_prefix_length + 1) {
			delete[] data_n[K]; data_n[K] = NULL;
			delete[] data_n_rc[K]; data_n_rc[K] = NULL;
		}

		foreach(K, max_p_prefix_length + 1) {
			delete[] data_p[K]; data_p[K] = NULL;
		}

		delete[] data_n; data_n = NULL;
		delete[] data_n_rc; data_n_rc = NULL;
		delete[] data_p; data_p = NULL;
	}

	foreach(i, MAX_RUNLEVEL) {
		foreach(j, (size_t)dc) {
			delete[] overall_scoretable[i][j];
			overall_scoretable[i][j] = NULL;
		}
		delete[] overall_scoretable[i];
	}

	foreach(i, NT_MM_K) {
		foreach(j, (size_t)dc) {
			delete[] scoretable_NT_MM[i][j]; scoretable_NT_MM[i][j] = NULL;
			delete[] scoretable_NT_MM_RC[i][j]; scoretable_NT_MM_RC[i][j] = NULL;
		}

		delete[] scoretable_NT_MM[i]; scoretable_NT_MM[i] = NULL;
		delete[] scoretable_NT_MM_RC[i]; scoretable_NT_MM_RC[i] = NULL;
	}

	delete[] scoretable_NT_MM; scoretable_NT_MM = NULL;
	delete[] scoretable_NT_MM_RC; scoretable_NT_MM_RC = NULL;

	foreach(i, AA_MM_K) {
		foreach(j, (size_t)dc) {
			delete[] scoretable_AA_MM[i][j]; scoretable_AA_MM[i][j] = NULL;
		}

		delete[] scoretable_AA_MM[i]; scoretable_AA_MM[i] = NULL;
	}

	delete[] scoretable_AA_MM; scoretable_AA_MM = NULL;

	delete[] has_protein; has_protein = NULL;

	delete[] class_16s; class_16s = NULL;
	delete[] class_ap; class_ap = NULL;

	foreach(i, MAX_VOTE) if(ENABLE_METHOD[i]) {
		foreach(j, (size_t)dc) {
			delete[] scoretable[i][j];
			scoretable[i][j] = NULL;
		}

		delete[] scoretable[i];
	}

	if(readlabels != NULL) {
		foreach(j, dc) {
			delete[] readlabels[j];
		}
		delete[] readlabels;
	}

	delete[] taxscores; taxscores = NULL;

	delete[] n_prefix_enabled;
	delete[] p_prefix_enabled;

	if(autocross_use || autocross_generate) autocross_cleanup();

	cerr << "[head] Done." << endl;
}
