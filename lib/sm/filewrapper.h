
#ifndef FILEWRAPPER_H
#define FILEWRAPPER_H

#ifdef _WIN32

#include <stdio.h>

inline FILE* myfopen(const char *fname, const char *mode) {
	return fopen(fname, mode);
}

inline FILE* myfopen(FILE *f, const char *mode) {
	return f;
}


#else

#include <stdio.h>
#include <stdarg.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>


typedef struct {
	bool gz;
	void *f;
} myfile;


inline myfile* myfopen(const char *file, const char *mode) {
	myfile *ret = new myfile;

	if (!strstr(file, ".gz")) {
		// regular file
		ret->gz = false;
		ret->f = fopen(file, mode);
	} else {
		// gzip file
		ret->gz = true;
		ret->f = gzopen(file, mode);
	}

	if (!ret->f) {
		delete ret;
		ret = NULL;
	}
		
	return ret;
}

inline myfile* myfopen(FILE *f, const char *mode) {
	myfile *ret = new myfile;
	ret->gz = false;
	ret->f = f;
	return ret;
}

inline void fclose(myfile *f) {
	if (f->gz)
		gzclose(f->f);
	else
		fclose((FILE*)(f->f));
	delete f;
}

inline void fprintf(myfile *f, const char *format, ...) {
	va_list ap;
	va_start(ap, format);
	char buffer[1024];
	vsnprintf(buffer, 1024, format, ap);
	va_end(ap);

	if (f->gz)
		gzprintf(f->f, "%s", buffer);
	else
		fprintf((FILE*)(f->f), "%s", buffer);
}

inline char* fgets(char *buf, int len, myfile *f) {
	if (f->gz)
		return gzgets(f->f, buf, len);
	else
		return fgets(buf, len, (FILE*)(f->f));
}

inline int getc(myfile *f) {
	if (f->gz)
		return gzgetc(f->f);
	else
		return getc((FILE*)(f->f));
}

inline int fputc(int c, myfile *f) {
	if (f->gz)
		return gzputc(f->f, c);
	else
		return fputc(c, (FILE*)(f->f));
}

inline int fgetc(myfile *f) {
	if (f->gz)
		return gzgetc(f->f);
	else
		return fgetc((FILE*)(f->f));
}


inline int fread(void *buf, int size, int num, myfile *f) {
	if (f->gz)
		return gzread(f->f, buf, size*num)/size;
	else
		return fread(buf, size, num, (FILE*)(f->f));
}


inline int fwrite(const void *buf, int size, int num, myfile *f) {
	if (f->gz)
		return gzwrite(f->f, buf, size*num);
	else
		return fwrite(buf, size, num, (FILE*)(f->f));
}

#undef FILE
#define FILE myfile




#endif



#endif



