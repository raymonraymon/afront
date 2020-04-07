#include <gtb/replace/replace.h>
#include <ctype.h>


int stricmp(const char *s1, const char *s2)
{
 	int n = 0;
	while ((n == 0) && (*s1 || *s2)) {
		n += tolower(*s1++) - tolower(*s2++);
	}
        return n;
}
