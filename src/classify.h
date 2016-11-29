/*
 	Gist classify.h

 		headers for classification manager

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

#ifndef CLASSIFY_H_
#define CLASSIFY_H_

void classifythreads(size_t num_threads);
void classify(size_t start_i, size_t end_i);

void runClock(size_t prog_status, size_t prog_target, time_t& starttime, string& prog_unit, string& prog_tag);

void prepclock(int process_target, string processing_units);
void prepclock(int process_target, string processing_units, string phase_tag);

void backupclock();
void restoreclock();

extern bool* has_protein; // [data]
extern number** scoretable[MAX_VOTE]; // scoretable[method][data][class]
extern number** overall_scoretable[MAX_RUNLEVEL]; // [runlevel-1][data][class]
extern number*** scoretable_NT_MM; // [cluster][data][class]
extern number*** scoretable_NT_MM_RC; // [cluster][data][class]
extern number*** scoretable_AA_MM; // [cluster][data][class]

#endif /* CLASSIFY_H_ */
