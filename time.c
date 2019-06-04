#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config-f90.h"
#endif

/* Return number of microseconds since 1.1.1970, in a 64 bit integer.
 * (with 2^64 us ~ 6 * 10^5 years, this should be sufficiently overflow safe)
 */
int64_t ftimings_microseconds_since_epoch(void) {
	struct timeval tv;
	if (gettimeofday(&tv, NULL) != 0) {
		perror("gettimeofday");
		exit(1);
	}
	return (int64_t) (tv.tv_sec) * ((int64_t) 1000000) + (int64_t)(tv.tv_usec);
}
