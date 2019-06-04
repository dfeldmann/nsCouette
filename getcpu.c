// obtain mapping of threads to cores (by F. Thomas, RZG/IBM) 
   #define _GNU_SOURCE
   #include <stdlib.h>
   # include <stdio.h>
   # include <math.h>
   # include <float.h>
   # include <omp.h>
   # include <limits.h>
   # include <sys/time.h>
   #include <sys/types.h>
   #include <sys/syscall.h>
   #include <unistd.h>
   #include <sched.h>
   #include <time.h>
   #include <pthread.h>
   #define gettid() syscall(SYS_gettid)
   static int
   getProcessorID(cpu_set_t* cpu_set)
   {
       int processorId;
       for (processorId=0;processorId<256;processorId++)
       {
           if (CPU_ISSET(processorId,cpu_set))
           {
               break;
           }
       }
       return processorId;
   }
   int  threadGetProcessorId()
   {
       cpu_set_t  cpu_set;
       CPU_ZERO(&cpu_set);
       sched_getaffinity(gettid(),sizeof(cpu_set_t), &cpu_set);
       return getProcessorID(&cpu_set);
   }
   int where_am_i_running()
   {
#ifdef _OPENMP 
   #pragma omp parallel
     {
       printf ("Thread %d running on processor %d ....\n",omp_get_thread_num(),threadGetProcessorId());
     }
     return threadGetProcessorId();
#else
     return 0;
#endif
   }

