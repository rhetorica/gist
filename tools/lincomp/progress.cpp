/*
 * progress.cpp
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */


#include <sstream>

#include "progress.h"

#define displayTime(t) (t / 3600) << ":" << (((t / 60) % 60) < 10 ? "0" : "") << ((t / 60) % 60) << ":" << ((t % 60) < 10 ? "0" : "") << (t % 60)

string prog_unit;
int prog_target = 0;
int prog_status = 0;
bool clock_running = false;

time_t starttime;

string prog_tag;

void prepclock(int pt, string pu, string ptt) {
	prog_unit = pu;
	prog_target = pt;
	prog_status = 0;
	prog_tag = ptt;
	clock_running = false;
	starttime = time(NULL);
}

void prepclock(int pt, string pu) {
	prepclock(pt, pu, string(""));
}

void runClock() {
	if(clock_running) return;
	clock_running = true;

	ostringstream rlog;

	try {
		double progfrac = (double)prog_status / (double)prog_target;
		double oldprogfrac = (double)(prog_status - 1) / (double)prog_target;

		if((int)(progfrac * 1000) > (int)(oldprogfrac * 1000)) {
			time_t elapsed = time(NULL) - starttime;
			double speed = 0;
			/*static time_t lasttime;
			if((int)(lasttime) == (int)(elapsed)) return;
			lasttime = elapsed;*/

			if(elapsed > 0) speed = ((double)prog_status / (double)elapsed); // * (double)((double)dc / (double)prog_target);
			if(speed > 10) speed = int(speed);
			else if(speed > 1) speed = (double)((int)(speed * 10.0)) / 10.0;
			else if(speed > 0.1) speed = (double)((int)(speed * 100.0)) / 100.0;
			else speed = (double)((int)(speed * 1000.0)) / 1000.0;

			time_t total = 0;
			if(prog_status > 0) total = elapsed / progfrac;

			time_t remaining = total * (1 - progfrac);

			std::streamsize ss = rlog.precision();

			rlog << "\x1b[2K\r";
			rlog.precision(1);

			if(remaining < 0 || prog_status < 1) {
				rlog << "(calculating ETA)...";
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
