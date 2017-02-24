/*
 * dircrawl.cpp
 *
 *  Created on: 2014-01-09
 *      Author: rhetorica
 */

#include "dircrawl.h"
#include "parselin.h"

vector<string> class_filenames;

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

bool is_dir(string& name) {
	struct stat buffer;
	stat(name.c_str(), &buffer);
	return S_ISDIR(buffer.st_mode);
//	return (buffer.st_mode & S_IFREG) > 0;
}

void crawl_tree(string& dirname, vector<string>* filenames) {
	cerr << "[crawl] Loading .lin files..." << endl;
	DIR *dir;
	struct dirent *ent = NULL;
	dir = opendir(dirname.c_str());
//	struct stat buf, bufg;

	unsigned long ccnt = 0;

	if(dir) {
		while((ent = readdir(dir))) {
			// if(ent->d_type == 4) continue;

			if(ent->d_name == string("..") || ent->d_name == string(".")) continue;

			string fname = string(dirname).append("/").append(string(ent->d_name));

			string ext4 = fname.substr(fname.length() - 5, 5);
			string ext3 = fname.substr(fname.length() - 4, 4);
			string ext2 = fname.substr(fname.length() - 3, 3);

			if(ext2 == string(".gz")) continue;
			if(ext2 == string(".sa")) continue;
			if(ext3 == string(".txt")) continue;
			if(ext3 == string(".amb")) continue;
			if(ext3 == string(".ann")) continue;
			if(ext3 == string(".bwt")) continue;
			if(ext3 == string(".pac")) continue;
			if(ext3 == string(".sam")) continue;
			if(ext3 == string(".bam")) continue;
			if(ext3 == string(".sai")) continue;
			if(ext3 == string(".fa")) continue;
			if(ext3 == string(".ffn")) continue;
			if(ext3 == string(".lin")) tax_id_from_lin(fname);
			if(ext4 == string(".gist")) continue;

			// cerr << "-> " << fname << " (" << (int)(ent->d_type) << ")" << endl;

			if(exists(fname) && is_file(fname)) {
//				cerr << "    " << ent->d_name << endl;
				class_filenames.push_back(fname);
				ccnt++;
			} else if(is_dir(fname)) {
				crawl_tree(fname, filenames);
			}

		}
	} else {
		cerr << "Could not open " << dirname << "." << endl;
		exit(85);
	}

	closedir(dir);
}
