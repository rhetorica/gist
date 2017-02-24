/*
 * main.h
 *
 *  Created on: 2013-08-08
 *      Author: rhetorica
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <cstring>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <algorithm>

#define foreach(i, max) for(size_t i = 0; i < max; i++)

using namespace std;

extern map<int, string> nrank;
extern map<int, int> parent;
extern map<int, string> nodename;

extern std::vector<std::string> class_filenames;

#endif /* MAIN_H_ */
