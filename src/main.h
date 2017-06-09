/*
 	Gist main.h

 		a giant bucket o' constants and compile-time parameters

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

#ifndef MAIN_H_
#define MAIN_H_

#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>
#include <cstddef>

#include <unistd.h>
#include <pthread.h>

using namespace std;

#define number double
// #define number float

// see also METHOD_NAMES in autocross.cpp
enum VOTETYPE {
	NT_NB, // 0
	AA_NB, // 1
	NT_NB_RC, // 2
	NT_NN, // 3
	AA_NN, // 4
	NT_NN_RC, // 5
	NT_CD, // 6
	AA_CD, // 7
	NT_CD_RC, // 8
	NT_MM, // 9
	AA_MM, // 10
	NT_MM_RC, // 11
	NT_SB, // 12
	AA_SB, // 13
	NT_SB_RC, // 14
	NT_BWA, // 15
	PRIOR_16S, // 16
	PRIOR_AP, // 17
	PRIOR_BIAS, // 18
	MAX_VOTE
};

#define BUILD_TIME __DATE__ " " __TIME__

#define GIST_VERSION "0.8.02 asterisks + fgs 1.16 + alglib 3.8.0 + autocross 3c (unweighted priors)"
// #define GIST_VERSION "0.8.01 you are everything + fgs 1.16 + alglib 3.8.0 + autocross 3c (unweighted priors)"
// #define GIST_VERSION "0.8.00 don't panic + fgs 1.16 + alglib 3.8.0 + autocross 3c (unweighted priors)"
// #define GIST_VERSION "0.7.17 shipwrecked and comatose + fgs 1.16 + alglib 3.8.0 + autocross 3b (no count regularization)"
// #define GIST_VERSION "0.7.16 thundersnail + fgs 1.16 + alglib 3.8.0 + autocross 3b (no count regularization)"
// #define GIST_VERSION "0.7.15 tenebre MCMLXXXIV + fgs 1.16 + alglib 3.8.0 + autocross 4"
// #define GIST_VERSION "0.7.14 statistics + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.13 wisdom to elisande + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.12 the witchmother commands it + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.11 chemical X + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.10 paprika + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.09 Experiencing a Gravitas Shortfall + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.08 the return of the thin white duke + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.07 trust the source + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.06 back to the valgrind + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.05 eyes just like mine + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.04 interceptor + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.03 ah, slippers? + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.02 paper packages tied up with string + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.01 pop-music theorems + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.7.00 akasha's euphemistic scepter + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04g2 shaggier prometheus + fgs 1.16 + alglib 3.8.0 + autocross 3a (no dog-ear)"
// #define GIST_VERSION "0.6.04g shaggy prometheus + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04f flawless theory + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04e electric feel + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04d the final countdown! + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04c the final countdown? + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04b prototype + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04a watch out + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.04 edge of disgrace + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.03a carry on + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.03 man overboard + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.02 roadkill is the path to awe + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.01 deus ex machina + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.6.00 supreme tango koala + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.99b bits bucketed + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.99a walk like an egyptian + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.99 we have accidentally borrowed your votedisk + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.98 five faces + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.97 goin' down the fast way + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.13 and do it more and more + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12f fake electronic lightshow + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12e one more time + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12d glass moon + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12c dog-earred disaster + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12b sugar rush + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12a waking dreams + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.12 anyone, anything + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.11 blistering sky + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.10 let's roll + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9g hilariously inefficient + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9f syncopation + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9e personal responsibility + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9d the light of day + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9c sublimation + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9b everything you ever did is coming back around + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9a the queen bee keeper + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.9 no wonder + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.8b space melodies + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.8a wumbaloo + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.8 the action of the clown + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.7a superb owl sunday + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.7 oregano rollercoaster + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.6 skyfall + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.5a actually a good idea + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.5 #ifndef what_the_hell_was_i_thinking + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.4 spin to win + fgs 1.16 + alglib 3.8.0 + autocross 3a"
// #define GIST_VERSION "0.5.3a ain't got time to take a fast train + fgs 1.16 + alglib 3.8.0 + autocross 3"
// #define GIST_VERSION "0.5.3 klaatu barada nikto + fgs 1.16 + alglib 3.8.0 + autocross 3"
// #define GIST_VERSION "0.5.2 disintegration of the persistence of network attached storage + fgs 1.16 + alglib 3.8.0 + autocross 3"
// #define GIST_VERSION "0.5.1 under your spell + fgs 1.16 + alglib 3.8.0 + autocross 3"
// #define GIST_VERSION "0.5.0 final flight of the osiris + fgs 1.16 + alglib 3.8.0 + autocross 3"
// #define GIST_VERSION "0.4.9 final flight of the osiris + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.8c just about half past ten + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.8b revenge of the grapevine + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.8a animus vox + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.8 genesis + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.7c bouton de rouge + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.7b far too loud + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.7a a really huge chicken with teeth + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.7 let's go bowling + fgs 1.16 + alglib 3.8.0 + autocross 3 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.6 through the grapevine + fgs 1.16 + alglib 3.8.0 + autocross 2a (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.5 skyline + fgs 1.16 + alglib 3.8.0 + autocross 2 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.5 skyline + fgs 1.16 + alglib 3.8.0 + autocross 1 (" __DATE__ " " __TIME__ ")"
// #define GIST_VERSION "0.4.4 the ones they'll blame + fgs 1.16 + alglib 3.8.0 + autocross 1"
// #define GIST_VERSION "0.4.3a they're back + fgs 1.16 + alglib 3.8.0"
// #define GIST_VERSION "0.4.3 this is the sound of C + fgs 1.16 + alglib 3.8.0"
// #define GIST_VERSION "0.4.2b don't stop the crock + fgs 1.16 + alglib 3.8.0"
// #define GIST_VERSION "0.4.2a don't stop the rock + fgs 1.16 + alglib 3.8.0"
// #define GIST_VERSION "0.4.2 find your light + fgs 1.16 + alglib 3.8.0"
// #define GIST_VERSION "0.4.1 running in the nineties + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.4.0d beginning today + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.4.0c more bells and whistles + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.4.0b must be the feeling + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.4.0a spaghetti + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.4.0 invisible pink unit testing + fgs 1.16 + alglib 3.8.0 + autocross 0"
// #define GIST_VERSION "0.3.2 octave multiplexer + fraggenescan 1.16"
// #define GIST_VERSION "0.3.1 electric mistress + fraggenescan 1.16"
// #define GIST_VERSION "0.3.0 remain indoors + fraggenescan 1.16"
// #define GIST_VERSION "0.2.6 with the lemons + fraggenescan 1.16"
// #define GIST_VERSION "0.2.5 could be a replicant + fraggenescan 1.16"
// #define GIST_VERSION "0.2.4 ameliorated artichoke + fraggenescan 1.16"
// #define GIST_VERSION "0.2.3 fluffy little clouds + fraggenescan 1.16"
// #define GIST_VERSION "0.2.2 celebrity fame zeppelin + fraggenescan 1.16"
// #define GIST_VERSION "0.2.1 pink fluffy unicorns dancing on rainbows + fraggenescan 1.16"
// #define GIST_VERSION "0.2.0 bloodshot streets + fraggenescan 1.16"
// #define GIST_VERSION "0.1.8 a whim away + fraggenescan 1.16"
// #define GIST_VERSION "0.1.7 baptismal records of death + fraggenescan 1.16"
// #define GIST_VERSION "0.1.6 eternal canned peaches + fraggenescan 1.16"
// #define GIST_VERSION "0.1.5 calm down dear + fraggenescan 1.16"
// #define GIST_VERSION "0.1.4 crazy ant conspiracy + fraggenescan 1.16"
// #define GIST_VERSION "0.1.3 making bacon"
// #define GIST_VERSION "0.1.2 poundcake"
// #define GIST_VERSION "0.1.1 artillery (100%)"
// #define GIST_VERSION "0.1.0 silly sibyl"
// #define GIST_VERSION "0.0.6 enlightenment"
// #define GIST_VERSION "0.0.5 god of hellfire"
// #define GIST_VERSION "0.0.4 conte hit and run"
// #define GIST_VERSION "0.0.3 infamous parallax"
// #define GIST_VERSION "0.0.2 linoleum floor"
// #define GIST_VERSION "0.0.1 blooming feeling"

// never report a distance metric value smaller than this after inversion:

#define MINV 0.00001

// enable dog-ear codelta optimization? (cuts memory usage by codelta in half)

#define DOG_EAR

// autocross option: return the weights that got the best classification performance, NOT the best error margin

#define SCORES_BEFORE_LOWS

// use D_2^* correlation for 1NN instead of euclidean

#define NN_D2_STAR

// maximum number of runlevels (hardcoded in some places; don't just change this!)

#define MAX_RUNLEVEL 2

// current program state:

extern size_t RUNLEVEL;

class Table;
class SparseTable;


// density for FGS:

extern int FGS_THRESHOLD;

// Component counts for mixture models:

extern size_t NT_MM_K;
extern size_t AA_MM_K;

// Iterations for mixture model training:
#define MM_ITER 120

char twobits(char e, int i);
char fivebits(char e);
int readFasta(string filename, char**& records, bool skim_mode);

char* translate(char* nts, bool class_mode);
char* reverse_complement(char* nrs);

void wait();

int sgn(number val);
number logsumexp(number* xs, size_t count);

typedef enum {
	NUCLEOTIDE = 0,
	PROTEIN = 1
} tableType;

extern bool use_class_threads;
extern int num_threads;

bool exists(string& name);

extern unsigned int class_read_length;

extern const char* codons;
extern const char* nucleotides;
extern const char* peptides;

extern bool outclasses;

extern bool disk_mode;
extern size_t batch_size;
extern size_t current_batch;

extern bool no_sqrt;

extern bool debugdump;
extern bool debuggmm;
extern bool debuglse;

extern size_t n_prefix_length;
extern size_t p_prefix_length;

// re-weights (from console input):

extern number bwa_weight;

extern number bn_weight;
extern number bp_weight;

extern number nn_weight;
extern number np_weight;

extern number cn_weight;
extern number cp_weight;

extern number mn_weight;
extern number mp_weight;

extern number sn_weight;
extern number sp_weight;

extern number p16s_weight;
extern number ap_weight;

// prefix lengths (from config file; changeable by console options):

extern size_t bn_prefix_length;
extern size_t bp_prefix_length;

extern size_t nn_prefix_length;
extern size_t np_prefix_length;

extern size_t cn_prefix_length;
extern size_t cp_prefix_length;

extern size_t mn_prefix_length;
extern size_t mp_prefix_length;

extern size_t sn_prefix_length;
extern size_t sp_prefix_length;

extern size_t max_n_prefix_length;
extern size_t max_p_prefix_length;

extern bool* n_prefix_enabled;
extern bool* p_prefix_enabled;

extern size_t derands; // how many bad nucleotides have we randomly guessed values for?

// computed method enabledness for this pass:

extern bool ADD_METHOD[MAX_VOTE]; // determine scores for these methods in this pass (for classify.cpp)
extern bool USE_METHOD[MAX_VOTE]; // sum scores for these methods in this pass (for score.cpp)
extern bool ENABLE_METHOD[MAX_VOTE]; // all methods ever considered by this session (global)
extern bool UNWEIGHTED_METHOD[MAX_VOTE]; // ignore autocross weightings for this method

extern size_t cc; // class count - how many classes are we working with?

extern std::string *dataname; // data file name - what data are we working with?

extern int dc; // data record count - how many reads are we working with?
extern char** dv; // data - every even-numbered row is a label, every odd-numbered row is a sequence

// base weights (from configuration file):

extern number method_weights[MAX_RUNLEVEL][MAX_VOTE];

// phylogenetic reconstruction parameters (from configuration file, can be overridden at command line):

extern size_t s1_include_quota;
extern number s1_include_threshold;
extern number s1_bloom_threshold;

extern size_t s2_include_quota;
extern number s2_include_threshold;
extern number s2_bloom_threshold;

// prior scores:

extern number *class_16s;
extern number *class_ap;

// calculated sequence markov tables:

extern Table*** data_n;
extern Table*** data_n_rc;
extern Table*** data_p;

// markov tables for sparse bayes:

extern SparseTable** sparse_data_n;
extern SparseTable** sparse_data_n_rc;
extern SparseTable** sparse_data_p;

// phylogenetic + output.cpp stuff

extern bool** shortlist; // [cc][dc] (who's a hit for what reads?)
extern size_t* shortlist_histo; // [cc] (how often are classes considered important?)

extern size_t **taxa; // [cc][MAX_RANK] total taxon_id database from .lin files
extern size_t *strain_taxon; // [cc] strain taxa for each class -- this information is not in taxa[][]

enum output_t {
	OUTPUT_READABLE, // output normal FASTA-like human-readable reports
	OUTPUT_CSV, // output neat CSV reports
	OUTPUT_SCORES, // output raw aggregate scores with no final phylogenetic processing
	OUTPUT_SCORETABLES, // output raw scoretables with no aggregation or phylogeny
	OUTPUT_RECORD // execute one pass for all methods, and dump scoretables in intermediate format
};

extern output_t output_mode;

// taxon levels in taxa[][]:

enum TAXON_RANKS {
	SPECIES,
	GENUS,
	FAMILY,
	ORDER,
	CLASS,
	PHYLUM,
	SUPERKINGDOM,
	MAX_RANK,
	NO_RANK
};

// autocross stuff:

extern number* autocross_CW; // [cc] class weights - deprecated
extern number* autocross_MW; // [MAX_VOTE] method weights - deprecated

extern int permissable_taxrank; // include how many taxonomic ranks as correct label? (-1 = strain, 0 = species, 1 = genus, etc.)

// data table maintenance and manipulation:

void clear_data_tables(size_t first, size_t stop); // from first to stop - 1, inclusive
void gen_data_tables(size_t first, size_t stop); // from first to stop - 1, inclusive
void gen_data_tables(); // do all

// convenience macros:

#define foreach(i, max) for(size_t i = 0; i < (size_t)(max); i++)
#define forclassifiers(i) for(size_t i = 0; i < MAX_VOTE; i++) if(USE_METHOD[i])
#define displayTime(t) (t / 3600) << ":" << (((t / 60) % 60) < 10 ? "0" : "") << ((t / 60) % 60) << ":" << ((t % 60) < 10 ? "0" : "") << (t % 60)

#endif
