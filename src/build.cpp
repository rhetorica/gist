/*
 	Gist build.cpp

 		execution management for class construction

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
#include "build.h"
#include <cstring>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <pthread.h>

struct threadparams {
	string filename;
	Class** container;
	size_t i;
};

// force resaving for debug purposes
#define FORCE_RESAVE false

void* runclassthread(void* t) {
	threadparams* ts = (threadparams*)t;
	try {
		ts->container[ts->i] = new Class(ts->filename);
	} catch(string* e) {
		cerr << *e << endl;
		return (void*)1;
	}
	try {
		if(outclasses) {
			if(ts->filename.substr(ts->filename.length() - 5, 5) != string(".gist") || FORCE_RESAVE) {
				string sfn = ts->filename;
				sfn.append(".gist");
				ts->container[ts->i]->saveClass(sfn);
				if(dataname == NULL) {
					std::cerr << "[head] Running generation only; de-allocating class " << sfn << " to save memory." << std::endl;
					delete ts->container[ts->i];
					ts->container[ts->i] = NULL;
				}
			}
		}
	} catch(string* e) {
		cerr << *e << endl;
		delete ts;
		return (void*)2;
	}
	// std::cerr << endl << "[thread] Class construction thread finished." << std::endl;
	delete ts;
	return NULL;
}

void oneclass(vector<string>* class_filenames, Class** classes, size_t i) {
	// void* status;

	std::cerr << "[head] Multithreading for class construction disabled." << std::endl;

	threadparams* thrps = new threadparams();
	thrps->filename = class_filenames->at(i);
	thrps->i = i;
	thrps->container = classes;
	void* thrpsy = (void*)thrps;

	long rc = (long)runclassthread(thrpsy);
	if(rc) {
		std::cerr << "[head] Unable to run class " << thrps->filename << ". Condition " << rc << std::endl;
		exit(124);
	}
}

void classnothreads(vector<string>* class_filenames, Class** classes) {
	// void* status;

	std::cerr << "[head] Multithreading for class construction disabled." << std::endl;

	size_t threads = class_filenames->size();

	for(size_t i = 0; i < threads; i++) {
		threadparams* thrps = new threadparams();
		thrps->filename = class_filenames->at(i);
		thrps->i = i;
		thrps->container = classes;
		void* thrpsy = (void*)thrps;

		long rc = (long)runclassthread(thrpsy);
		if(rc) {
			std::cerr << "[head] Unable to run class " << thrps->filename << ". Condition " << rc << std::endl;
			exit(123);
		}
	}
}

void classthreads(vector<string>* class_filenames, Class** classes) {
	void* status;

	std::cerr << "[head] Multithreading for class construction enabled." << std::endl;

	size_t threads = class_filenames->size();

	// Threads must be "joinable" in order to wait on them to finish:
	pthread_t threadids[threads];
	pthread_attr_t attribs;
	pthread_attr_init(&attribs);
	pthread_attr_setdetachstate(&attribs, PTHREAD_CREATE_JOINABLE);

	for(size_t i = 0; i < threads; i++) {
		threadparams* thrps = new threadparams();
		thrps->filename = class_filenames->at(i);
		thrps->i = i;
		thrps->container = classes;
		void* thrpsy = (void*)thrps;

		int rc = pthread_create(&threadids[i], NULL, runclassthread, thrpsy);
		if(rc) {
			std::cerr << "[head] Unable to create thread " << thrps->filename << ". Condition " << rc << std::endl;
			exit(120);
		}
	}

	pthread_attr_destroy(&attribs);

	// wait on each thread to finish:
	for(size_t i = 0; i < threads; i++) {
		int rc = pthread_join(threadids[i], &status);
		if(rc) {
			std::cerr << "[head] Unable to monitor thread " << i << ". Condition " << rc << std::endl;
			exit(121);
		}
		if((unsigned long)status > 0) {
			std::cerr << "[head] Thread " << i << " exiting with status " << (unsigned long)status << std::endl;
			exit(122);
		}
	}
}
