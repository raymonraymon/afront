
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/


#include <gtb/gtb.hpp>
#ifndef WIN32
#include <gtb/io/io.hpp>
#include <gtb/error/error.hpp>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/io/io.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


void get_file_base_name(const char *file_name, char *base_name, unsigned size)
{
	assert(NULL != file_name);
	assert(NULL != base_name);
	assert(size > 0);
	unsigned n = strlen(file_name);

	// degenerate case
	if (0 == n) {
		base_name[0] = '\0';
		return;
	}

	// find right end
	unsigned r = n - 1;
	while ((r > 0) && ('/' == file_name[r])) {
		r--;
	}

	// find left end
	unsigned l = r;
	while ((l > 0) && ('/' != file_name[l])) {
		l--;
	}
	if (('/' == file_name[l]) && (l < r)) {
		l++;
	}

	// check size
	if (r - l + 2 > size) {
		fprintf(stderr, "%s:%d: buffer for base name is too small\n",
			__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	// copy base name
	unsigned k = 0;
	for (unsigned i = l; i <= r; i++, k++) {
		assert(k < size - 1);
		base_name[k] = file_name[i];
	}
	base_name[k] = '\0';
}


void get_file_base_name(const char *file_name,
			const char *suffix,
			char *base_name,
			unsigned size)
{
	assert(NULL != file_name);
	assert(NULL != suffix);
	assert(NULL != base_name);
	assert(size > 0);
	int n = strlen(file_name);

	// degenerate case
	if (0 == n) {
		base_name[0] = '\0';
		return;
	}

	// find right end
	int r = n - 1;
	while ((r > 0) && ('/' == file_name[r])) {
		r--;
	}

	// check for suffix
	int m = strlen(suffix);
	if (r - (m - 1) >= 0) {
		const char *p = file_name + r - (m - 1);
		if (strcmp(p, suffix) == 0) {
			if (m == n) {
				base_name[0] = '\0';
				return;
			}
			assert(r >= m);
			r -= m;
		}
	}

	// find left end
	int l = r;
	while ((l > 0) && ('/' != file_name[l])) {
		l--;
	}
	if (('/' == file_name[l]) && (l < r)) {
		l++;
	}

	// check size
	if (r - l + 2 > (int) size) {
		fprintf(stderr, "%s:%d: buffer for base name is too small\n",
			__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	// copy base name
	assert(l >= 0);
	assert(r >= l);
	int k = 0;
	for (int i = l; i <= r; i++, k++) {
		assert(k < (int) size - 1);
		base_name[k] = file_name[i];
	}
	base_name[k] = '\0';
}


void get_file_extension(const char *file_name, char *extension, unsigned size)
{
	assert(NULL != file_name);
	assert(NULL != extension);
	assert(size > 0);
	unsigned n = strlen(file_name);

	// degenerate case
	if (0 == n) {
		extension[0] = '\0';
		return;
	}

	// check for directory
	if ('/' == file_name[n - 1]) {
		fprintf(stderr, "%s is a directory\n", file_name);
		exit(EXIT_FAILURE);
	}

	// find '.'
	unsigned r = n - 1;
	unsigned l = n - 1;
	while ((l > 0) && ('.' != file_name[l])) {
		l--;
	}

	// no extension case
	if ('.' != file_name[l]) {
		extension[0] = '\0';
		return;
	}

	// skip '.'
	l++;

	// check size
	if (r - l + 2 > size) {
		fprintf(stderr, "%s:%d: buffer for extension is too small\n",
			__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	// copy extension
	unsigned k = 0;
	for (unsigned i = l; i <= r; i++, k++) {
		assert(k < size - 1);
		extension[k] = file_name[i];
	}
	extension[k] = '\0';
}


FILE *xfopen(const char *file_name, const char *mode)
{
	assert(NULL != file_name);
	assert(NULL != mode);
	FILE *fp = fopen(file_name, mode);
	if (NULL == fp) {
		perror(file_name);
#ifndef NO_EXCEPTIONS
		throw CErr("xfopen: failed");
#else
		exit(EXIT_FAILURE);
#endif
	}
	return fp;
}

/*
 * Open a file for write
 * if the file already exist, save a version of the old file
 *
 * max_backup - maximal number of backup files to save
 * 
 * file backup format name is name.ext;version
 */
#if 0
FILE* xfopenv(const char* file_name, int max_backup)
{
#ifdef _WIN32
    if (file_exists(file_name))
    {
        std::vector<std::string> fileslist;

        char searchfile[1000];
        WIN32_FIND_DATA finddata;
        sprintf(searchfile, "%s;*", file_name);
        HANDLE h = FindFirstFile(searchfile, &finddata);
        if (h != INVALID_HANDLE_VALUE)
        {
            while (1)
            {
                fileslist.push_back(finddata.cFileName);
                if (!FindNextFile(h, &finddata)) break;
            }
            FindClose(h);

            const char* p = strrchr(fileslist.back().c_str(), ';');
            assert(p != 0);
            int last_backup_idx = atoi(p+1)+1;
            char newfilename[1024];
            sprintf(newfilename, "%s;%05d", file_name, last_backup_idx);
            if (!MoveFile(file_name, newfilename))
            {
                printf("Failed to move file to its backup aborintg open");
                return 0;
            }

            // delete old versions
            if ((max_backup > 0) && (fileslist.size() >= (unsigned)max_backup))
            {
                std::vector<std::string>::const_iterator f = fileslist.begin();
                std::vector<std::string>::const_iterator l = fileslist.end() - max_backup;
                for (; f != l; ++f)
                {
                    DeleteFile(f->c_str());
                }
            }
        }
    }

    return fopen(file_name, "wb");

#else
#error "xfopenv Not supported"
#endif
}
#endif

bool file_exists(const char *file_name)
{
	struct stat buf;
	if ((stat(file_name, &buf) == 0)
	    && S_ISREG(buf.st_mode)) {
		return true;
	}
	return false;
}


bool dir_exists(const char *dir_name)
{
	struct stat buf;
	if ((stat(dir_name, &buf) == 0)
	    && S_ISDIR(buf.st_mode)) {
		return true;
	}
	return false;
}


GTB_END_NAMESPACE
