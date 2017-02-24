/*
 * progress.h
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#ifndef PROGRESS_H_
#define PROGRESS_H_

#include <string>
#include "main.h"

extern string prog_unit;
extern int prog_target;
extern int prog_status;
extern bool clock_running;

extern time_t starttime;

extern string prog_tag;

void prepclock(int pt, string pu, string ptt);

void prepclock(int pt, string pu);

void runClock();


#endif /* PROGRESS_H_ */
