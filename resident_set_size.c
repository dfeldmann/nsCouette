#include <stdio.h>
#include <unistd.h>

long ftimings_resident_set_size() {
	long rss = 0L;
	FILE* fp = NULL;
	if ((fp = fopen( "/proc/self/statm", "r" )) == NULL ) {
		return 0L;
	}
	if (fscanf(fp, "%*s%ld", &rss) != 1) {
		fclose(fp);
		return (size_t)0L;	  /* Can't read? */
	}
	fclose(fp);
	return rss * sysconf( _SC_PAGESIZE);
}
