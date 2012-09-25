/*
 * DiskSim Storage Subsystem Simulation Environment (Version 4.0)
 * Revision Authors: John Bucy, Greg Ganger
 * Contributors: John Griffin, Jiri Schindler, Steve Schlosser
 *
 * Copyright (c) of Carnegie Mellon University, 2001-2008.
 *
 * This software is being provided by the copyright holders under the
 * following license. By obtaining, using and/or copying this software,
 * you agree that you have read, understood, and will comply with the
 * following terms and conditions:
 *
 * Permission to reproduce, use, and prepare derivative works of this
 * software is granted provided the copyright and "No Warranty" statements
 * are included with all reproductions and derivative works and associated
 * documentation. This software may also be redistributed without charge
 * provided that the copyright and "No Warranty" statements are included
 * in all redistributions.
 *
 * NO WARRANTY. THIS SOFTWARE IS FURNISHED ON AN "AS IS" BASIS.
 * CARNEGIE MELLON UNIVERSITY MAKES NO WARRANTIES OF ANY KIND, EITHER
 * EXPRESSED OR IMPLIED AS TO THE MATTER INCLUDING, BUT NOT LIMITED
 * TO: WARRANTY OF FITNESS FOR PURPOSE OR MERCHANTABILITY, EXCLUSIVITY
 * OF RESULTS OR RESULTS OBTAINED FROM USE OF THIS SOFTWARE. CARNEGIE
 * MELLON UNIVERSITY DOES NOT MAKE ANY WARRANTY OF ANY KIND WITH RESPECT
 * TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.
 * COPYRIGHT HOLDERS WILL BEAR NO LIABILITY FOR ANY USE OF THIS SOFTWARE
 * OR DOCUMENTATION.
 *
 */



/*
 * DiskSim Storage Subsystem Simulation Environment (Version 2.0)
 * Revision Authors: Greg Ganger
 * Contributors: Ross Cohen, John Griffin, Steve Schlosser
 *
 * Copyright (c) of Carnegie Mellon University, 1999.
 *
 * Permission to reproduce, use, and prepare derivative works of
 * this software for internal use is granted provided the copyright
 * and "No Warranty" statements are included with all reproductions
 * and derivative works. This software may also be redistributed
 * without charge provided that the copyright and "No Warranty"
 * statements are included in all redistributions.
 *
 * NO WARRANTY. THIS SOFTWARE IS FURNISHED ON AN "AS IS" BASIS.
 * CARNEGIE MELLON UNIVERSITY MAKES NO WARRANTIES OF ANY KIND, EITHER
 * EXPRESSED OR IMPLIED AS TO THE MATTER INCLUDING, BUT NOT LIMITED
 * TO: WARRANTY OF FITNESS FOR PURPOSE OR MERCHANTABILITY, EXCLUSIVITY
 * OF RESULTS OR RESULTS OBTAINED FROM USE OF THIS SOFTWARE. CARNEGIE
 * MELLON UNIVERSITY DOES NOT MAKE ANY WARRANTY OF ANY KIND WITH RESPECT
 * TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.
 */

/*
 * DiskSim Storage Subsystem Simulation Environment
 * Authors: Greg Ganger, Bruce Worthington, Yale Patt
 *
 * Copyright (C) 1993, 1995, 1997 The Regents of the University of Michigan 
 *
 * This software is being provided by the copyright holders under the
 * following license. By obtaining, using and/or copying this software,
 * you agree that you have read, understood, and will comply with the
 * following terms and conditions:
 *
 * Permission to use, copy, modify, distribute, and sell this software
 * and its documentation for any purpose and without fee or royalty is
 * hereby granted, provided that the full text of this NOTICE appears on
 * ALL copies of the software and documentation or portions thereof,
 * including modifications, that you make.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS," AND COPYRIGHT HOLDERS MAKE NO
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE,
 * BUT NOT LIMITATION, COPYRIGHT HOLDERS MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
 * THAT THE USE OF THE SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY
 * THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. COPYRIGHT
 * HOLDERS WILL BEAR NO LIABILITY FOR ANY USE OF THIS SOFTWARE OR
 * DOCUMENTATION.
 *
 *  This software is provided AS IS, WITHOUT REPRESENTATION FROM THE
 * UNIVERSITY OF MICHIGAN AS TO ITS FITNESS FOR ANY PURPOSE, AND
 * WITHOUT WARRANTY BY THE UNIVERSITY OF MICHIGAN OF ANY KIND, EITHER
 * EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE REGENTS
 * OF THE UNIVERSITY OF MICHIGAN SHALL NOT BE LIABLE FOR ANY DAMAGES,
 * INCLUDING SPECIAL , INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 * WITH RESPECT TO ANY CLAIM ARISING OUT OF OR IN CONNECTION WITH THE
 * USE OF OR IN CONNECTION WITH THE USE OF THE SOFTWARE, EVEN IF IT HAS
 * BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY OF SUCH DAMAGES
 *
 * The names and trademarks of copyright holders or authors may NOT be
 * used in advertising or publicity pertaining to the software without
 * specific, written prior permission. Title to copyright in this software
 * and any associated documentation will at all times remain with copyright
 * holders.
 */

#include "config.h"

#include "disksim_global.h"
#include "disksim_hptrace.h"
#include "disksim_iotrace.h"

#include "disksim_iodriver.h"
#include "disksim_logorg.h"
#include "../ssdmodel/ssd.h"
#include "../ssdmodel/ssd_timing.h"
#include<sys/stat.h> //for mkdir
#ifdef DEBUG
#include<time.h>
#endif
static void realign_request(ioreq_event *new_ioreq_event);
int time_to_checkpoint(double currtime);

static void iotrace_initialize_iotrace_info ()
{
   disksim->iotrace_info = DISKSIM_malloc (sizeof(iotrace_info_t));
   bzero ((char *)disksim->iotrace_info, sizeof(iotrace_info_t));

   tracebasetime = 0.0;
   firstio = TRUE;
   lasttime = 0.0;
   basebigtime = -1;
   basesmalltime = -1;
   basesimtime = 0.0;
   validate_lastserv = 0.0;
   validate_nextinter = 0.0;
   accumulated_event_time = 0.0;
   lastaccesstime = 0.0;
}


void iotrace_set_format (char *formatname)
{
   if (disksim->iotrace_info == NULL) {
      iotrace_initialize_iotrace_info();
   }

   disksim->traceendian = _LITTLE_ENDIAN;
   disksim->traceheader = TRUE;
   if (strcmp(formatname, "0") == 0) {
      disksim->traceformat = DEFAULT;
   } else if (strcmp(formatname, "ascii") == 0) {
	/* default ascii trace format */
      disksim->traceformat = ASCII;
   } else if (strcmp(formatname, "raw") == 0) {
	/* format of traces collected at Michigan */
      disksim->traceformat = RAW;
   } else if (strcmp(formatname, "validate") == 0) {
	/* format of disk validation traces */
      disksim->traceformat = VALIDATE;
   } else if (strcmp(formatname, "hpl") == 0) {
	/* format of traces provided by HPLabs for research purposes */
      disksim->traceformat = HPL;
      disksim->traceendian = _BIG_ENDIAN;
   } else if (strcmp(formatname, "hpl2") == 0) {
	/* format of traces provided by HPLabs for research purposes,     */
        /* after they have been modified/combined by the `hplcomb' program */
      disksim->traceformat = HPL;
      disksim->traceendian = _BIG_ENDIAN;
      disksim->traceheader = FALSE;
   } else if (strcmp(formatname, "dec") == 0) {
	/* format of some traces provided by dec for research purposes */
      disksim->traceformat = DEC;
   } else if (strcmp(formatname, "emcsymm") == 0) {
	/* format of Symmetrix traces provided by EMC for research purposes */
      disksim->traceformat = EMCSYMM;
   } else if (strcmp(formatname, "emcbackend") == 0) {
        /* format of Symmetrix backend traces that have been converted to map
	   logical devices to physical devices  */
      disksim->traceformat = EMCBACKEND;
   } else if (strcmp(formatname, "batch") == 0) {
        /* ascii traces with added batch information */
      disksim->traceformat = BATCH;
   } else {
      fprintf(stderr, "Unknown trace format - %s\n", formatname);
      exit(1);
   }
}


static int iotrace_read_space (FILE *tracefile, char *spaceptr, int spacesize)
{
   if (fread(spaceptr, spacesize, 1, tracefile) != 1) {
      return(-1);
   }
   return(0);
}


static int iotrace_read_char (FILE *tracefile, char *charptr)
{
   StaticAssert (sizeof(char) == 1);
   if (fread(charptr, sizeof(char), 1, tracefile) != 1) {
      return(-1);
   }
   return(0);
}


static int iotrace_read_short (FILE *tracefile, short *shortptr)
{
   StaticAssert (sizeof(short) == 2);
   if (fread(shortptr, sizeof(short), 1, tracefile) != 1) {
      return(-1);
   }
   if (disksim->endian != disksim->traceendian) {
      *shortptr = ((*shortptr) << 8) + (((*shortptr) >> 8) & 0xFF);
   }
   return(0);
}


static int iotrace_read_int32 (FILE *tracefile, int32_t *intP)
{
   int i;
   intchar swapval;
   intchar intcharval;

   StaticAssert (sizeof(int) == 4);
   if (fread(&intcharval.value, sizeof(int), 1, tracefile) != 1) {
      return(-1);
   }
   if (disksim->endian != disksim->traceendian) {
      for (i=0; i<sizeof(int); i++) {
         swapval.byte[i] = intcharval.byte[(sizeof(int) - i - 1)];
      }
/*
      fprintf (outputfile, "intptr.value %x, swapval.value %x\n", intcharval.value, swapval.value);
*/
      intcharval.value = swapval.value;
   }
   *intP = intcharval.value;
   return(0);
}


#define iotrace_read_float(a, b) iotrace_read_int32(a, b)


ioreq_event * iotrace_validate_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   char line[201];
   char rw;
   double servtime;

   if (fgets(line, 200, tracefile) == NULL) {
      addtoextraq((event *) new_ioreq_event);
      return(NULL);
   }
   new_ioreq_event->time = simtime + (validate_nextinter / (double) 1000);
   if (sscanf(line, "%c %s %ld %d %lf %lf\n", 
	      &rw, 
	      validate_buffaction, 
	      &new_ioreq_event->blkno, 
	      &new_ioreq_event->bcount, 
	      &servtime, 
	      &validate_nextinter) != 6) 
   {
      fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
      ddbg_assert(0);
   }
   validate_lastserv = servtime / (double) 1000;
   if (rw == 'R') {
      new_ioreq_event->flags = READ;
   } else if (rw == 'W') {
      new_ioreq_event->flags = WRITE;
   } else {
      fprintf(stderr, "Invalid R/W value: %c\n", rw);
      exit(1);
   }
   new_ioreq_event->devno = 0;
   new_ioreq_event->buf = 0;
   new_ioreq_event->opid = 0;
   new_ioreq_event->cause = 0;
   new_ioreq_event->busno = 0;
   new_ioreq_event->tempint2 = 0;
   new_ioreq_event->tempint1 = 0;
   validate_lastblkno = new_ioreq_event->blkno;
   validate_lastbcount = new_ioreq_event->bcount;
   validate_lastread = new_ioreq_event->flags & READ;
   return(new_ioreq_event);
}


static ioreq_event * iotrace_dec_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   assert ("removed for distribution" == 0);
   iotrace_read_space (tracefile, NULL, 0);
   return (NULL);
}


static void iotrace_hpl_srt_convert_flags (ioreq_event *curr)
{
   int flags;

   flags = curr->flags;
   curr->flags = 0;
   if (flags & HPL_READ) {
      curr->flags |= READ;
      hpreads++;
   } else {
      hpwrites++;
   }
   if (!(flags & HPL_ASYNC)) {
      curr->flags |= TIME_CRITICAL;
      if (curr->flags & READ) {
         syncreads++;
      } else {
         syncwrites++;
      }
   }
   if (flags & HPL_ASYNC) {
      if (curr->flags & READ) {
         curr->flags |= TIME_LIMITED;
         asyncreads++;
      } else {
         asyncwrites++;
      }
   }
}


static ioreq_event * iotrace_hpl_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   int32_t size;
   int32_t id=0;
   int32_t sec=0;
   int32_t usec;
   int32_t val=0;
   int32_t junkint;
   unsigned int failure = 0;

   while (TRUE) {
      failure |= iotrace_read_int32(tracefile, &size);
      failure |= iotrace_read_int32(tracefile, &id);
      failure |= iotrace_read_int32(tracefile, &sec);
      failure |= iotrace_read_int32(tracefile, &usec);
      if (failure) {
         addtoextraq((event *) new_ioreq_event);
         return(NULL);
      }
      if (((id >> 16) < 1) || ((id >> 16) > 4)) {
         fprintf(stderr, "Error in trace format - id %x\n", id);
         exit(1);
      }
      if (((id & 0xFFFF) != HPL_SHORTIO) && ((id & 0xFFFF) != HPL_SUSPECTIO)) {
         fprintf(stderr, "Unexpected record type - %x\n", id);
         exit(1);
      }
      new_ioreq_event->time = (double) sec * (double) MILLI;
      new_ioreq_event->time += (double) usec / (double) MILLI;

      if ((disksim->traceheader == FALSE) && (new_ioreq_event->time == 0.0)) {
         tracebasetime = simtime;
      }

      failure |= iotrace_read_int32(tracefile, &val);    /* traced request start time */
      new_ioreq_event->tempint1 = val;
      failure |= iotrace_read_int32(tracefile, &val);    /* traced request stop time */
      new_ioreq_event->tempint2 = val;
      new_ioreq_event->tempint2 -= new_ioreq_event->tempint1;
      failure |= iotrace_read_int32(tracefile, &val);
      new_ioreq_event->bcount = val;
      if (new_ioreq_event->bcount & 0x000001FF) {
         fprintf(stderr, "HPL request for non-512B multiple size: %d\n", new_ioreq_event->bcount);
         exit(1);
      }
      new_ioreq_event->bcount = new_ioreq_event->bcount >> 9;
      failure |= iotrace_read_int32(tracefile, &val);
      new_ioreq_event->blkno = val;
      failure |= iotrace_read_int32(tracefile, &val);
      new_ioreq_event->devno = (val >> 8) & 0xFF;
      failure |= iotrace_read_int32(tracefile, &val);       /* drivertype */
      failure |= iotrace_read_int32(tracefile, &val);       /* cylno */
      /* for convenience and historical reasons, this cast is being allowed */
      /* (the value is certain to be less than 32 sig bits, and will not be */
      /* used as a pointer).                                                */
      new_ioreq_event->buf = (void *)(long) val;
      failure |= iotrace_read_int32(tracefile, &val);
      new_ioreq_event->flags = val;
      iotrace_hpl_srt_convert_flags(new_ioreq_event);
      failure |= iotrace_read_int32(tracefile, &junkint);           /* info */
      size -= 13 * sizeof(int32_t);
      if ((id >> 16) == 4) {
         failure |= iotrace_read_int32(tracefile, &val);  /* queuelen */
         new_ioreq_event->slotno = val;
         size -= sizeof(int32_t);
      }
      if ((id & 0xFFFF) == HPL_SUSPECTIO) {
         failure |= iotrace_read_int32(tracefile, &junkint);    /* susflags */
         size -= sizeof(int32_t);
      }
      if (failure) {
         addtoextraq((event *) new_ioreq_event);
         return(NULL);
      }
      if (size) {
         fprintf(stderr, "Unmatched size for record - %d\n", size);
         exit(1);
      }
      new_ioreq_event->cause = 0;
      new_ioreq_event->opid = 0;
      new_ioreq_event->busno = 0;
      if ((id & 0xFFFF) == HPL_SHORTIO) {
         return(new_ioreq_event);
      }
   }
}


static int iotrace_month_convert (char *monthstr, int year)
{
   if (strcmp(monthstr, "Jan") == 0) {
      return(0);
   } else if (strcmp(monthstr, "Feb") == 0) {
      return(31);
   } else if (strcmp(monthstr, "Mar") == 0) {
      return((year % 4) ? 59 : 60);
   } else if (strcmp(monthstr, "Apr") == 0) {
      return((year % 4) ? 90 : 91);
   } else if (strcmp(monthstr, "May") == 0) {
      return((year % 4) ? 120 : 121);
   } else if (strcmp(monthstr, "Jun") == 0) {
      return((year % 4) ? 151 : 152);
   } else if (strcmp(monthstr, "Jul") == 0) {
      return((year % 4) ? 181 : 182);
   } else if (strcmp(monthstr, "Aug") == 0) {
      return((year % 4) ? 212 : 213);
   } else if (strcmp(monthstr, "Sep") == 0) {
      return((year % 4) ? 243 : 244);
   } else if (strcmp(monthstr, "Oct") == 0) {
      return((year % 4) ? 273 : 274);
   } else if (strcmp(monthstr, "Nov") == 0) {
      return((year % 4) ? 304 : 305);
   } else if (strcmp(monthstr, "Dec") == 0) {
      return((year % 4) ? 334 : 335);
   }
   assert(0);
   return(-1);
}


static double iotrace_raw_get_hirestime (int bigtime, int smalltime)
{
   unsigned int loresticks;
   int small, turnovers;
   int smallticks;

   if (basebigtime == -1) {
      basebigtime = bigtime;
      basesmalltime = smalltime;
      basesimtime = 0.0;
   } else {
      small = (basesmalltime - smalltime) & 0xFFFF;
      loresticks = (bigtime - basebigtime) * 11932 - small;
      turnovers = (int) (((double) loresticks / (double) 65536) + (double) 0.5);
      smallticks = turnovers * 65536 + small;
      basebigtime = bigtime;
      basesmalltime = smalltime;
      basesimtime += (double) smallticks * (double) 0.000838574;
   }
   return(basesimtime);
}


/* kept mainly as an example */

static ioreq_event * iotrace_raw_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   int bigtime;
   short small;
   int smalltime;
   int failure = 0;
   char order, crit;
   double schedtime, donetime;
   int32_t val=0;

   failure |= iotrace_read_int32(tracefile, &val);
   bigtime = val;
   failure |= iotrace_read_short(tracefile, &small);
   smalltime = ((int) small) & 0xFFFF;
   new_ioreq_event->time = iotrace_raw_get_hirestime(bigtime, smalltime);
   failure |= iotrace_read_short(tracefile, &small);
   failure |= iotrace_read_int32(tracefile, &val);
   bigtime = val;
   smalltime = ((int) small) & 0xFFFF;
   schedtime = iotrace_raw_get_hirestime(bigtime, smalltime);
   failure |= iotrace_read_int32(tracefile, &val);
   bigtime = val;
   failure |= iotrace_read_short(tracefile, &small);
   smalltime = ((int) small) & 0xFFFF;
   donetime = iotrace_raw_get_hirestime(bigtime, smalltime);
   failure |= iotrace_read_char(tracefile, &order);
   failure |= iotrace_read_char(tracefile, &crit);
   if (crit) {
      new_ioreq_event->flags |= TIME_CRITICAL;
   }
   failure |= iotrace_read_int32(tracefile, &val);
   new_ioreq_event->bcount = val >> 9;
   failure |= iotrace_read_int32(tracefile, &val);
   new_ioreq_event->blkno = val;
   failure |= iotrace_read_int32(tracefile, &val);
   new_ioreq_event->devno = val;
   failure |= iotrace_read_int32(tracefile, &val);
   new_ioreq_event->flags = val & READ;
   new_ioreq_event->cause = 0;
   new_ioreq_event->buf = 0;
   new_ioreq_event->opid = 0;
   new_ioreq_event->busno = 0;
   new_ioreq_event->tempint1 = (int)((schedtime - new_ioreq_event->time) * (double) 1000);
   new_ioreq_event->tempint2 = (int)((donetime - schedtime) * (double) 1000);
   if (failure) {
      addtoextraq((event *) new_ioreq_event);
      new_ioreq_event = NULL;
   }
   return(new_ioreq_event);
}


static ioreq_event * iotrace_emcsymm_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   char line[201];
   char operation[15];
   unsigned int director;

   if (fgets(line, 200, tracefile) == NULL) {
      addtoextraq((event *) new_ioreq_event);
      return(NULL);
   }
   if (sscanf(line, "%lf %s %x %x %ld %d\n", &new_ioreq_event->time, operation, &director, &new_ioreq_event->devno, &new_ioreq_event->blkno, &new_ioreq_event->bcount) != 6) {
      fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
      fprintf(stderr, "line: %s", line);
      ddbg_assert(0);
   }
   if (!strcasecmp(operation,"Read")) {
      new_ioreq_event->flags = READ;
   } else if (!strcasecmp(operation,"Write")) {
      new_ioreq_event->flags = WRITE;
   } else {
      fprintf(stderr, "Unknown operation: %s in iotrace event\n",operation);
      fprintf(stderr, "line: %s", line);
      exit(1);
   }
   new_ioreq_event->buf = 0;
   new_ioreq_event->opid = 0;
   new_ioreq_event->busno = 0;
   new_ioreq_event->cause = 0;
   return(new_ioreq_event);
}

static ioreq_event * iotrace_emcbackend_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   char line[201];
   char operation[15];
   unsigned int director;
   char bus[2];
   unsigned int disk, hyper;

   if (fgets(line, 200, tracefile) == NULL) {
      addtoextraq((event *) new_ioreq_event);
      return(NULL);
   }

   if (sscanf(line, "%lf %s %x %x %ld %d %s %d %d\n", 
              &new_ioreq_event->time, operation, &director, &hyper, &new_ioreq_event->blkno, &new_ioreq_event->bcount, bus, &disk, &new_ioreq_event->devno) != 9) {
      fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
      fprintf(stderr, "line: %s", line);
      exit(0);
   }
   if (!strcasecmp(operation,"Read")) {
      new_ioreq_event->flags = READ;
   } else if (!strcasecmp(operation,"Write")) {
      new_ioreq_event->flags = WRITE;
   } else {
      fprintf(stderr, "Unknown operation: %s in iotrace event\n",operation);
      fprintf(stderr, "line: %s", line);
      exit(0);
   }

   new_ioreq_event->time *= 1000.0;  /* emc trace times are in seconds!!  */

   new_ioreq_event->buf = 0;
   new_ioreq_event->opid = 0;
   new_ioreq_event->busno = 0;
   new_ioreq_event->cause = 0;
   return(new_ioreq_event);
}

static int tot_req[8]= {0,0,0,0,0,0,0,0};
static int element_id = 3;

ioreq_event * iotrace_ascii_get_ioreq_event_detailed_simulation(FILE* tracefile, ioreq_event *new_ioreq_event)
{
	   char line[201];
	   static double last_req_time = 0;
	   static int is_first_request = 1;
	   int i=0;
	   int sub_status = 1;
	   //static current_event = trace_event_list;

	   //by abhishek
	   static int max_count;
	   max_count = disksim->repeat_trace;
	   static int min_count = 0;
	   static double start_time = 0;
           // KJ: The following three params are set to get the duration of a detailed simulation round;
           // Used in temporal granularity analysis
           static double prev_round_end = 0;
           static double curr_round_end = 0;
           static int record_prev = 1;
           static int record_curr = 1;
           assert(tracefile != NULL);

                  
           //printf("Reaches here!\n");
	start_again:
		if (fgets(line, 200, tracefile) == NULL) {
		      
                  
		  //printf("Reaches here!\n");
                //by abhishek
		      min_count++;
                disksim->cur_detailed_trace_count++;
                assert(disksim->cur_detailed_trace_count == min_count);
	        fprintf(stderr,"Pass %d completed\n",min_count);
                
	#ifdef DEBUG
	        logorg **temp = sysorgs;
	        					 fprintf(stderr,"Curr:%d,Max:%d,Qlen:%d\n",temp[0]->stat.outstanding,temp[0]->stat.maxoutstanding,temp[0]->outstandqlen);
	        fflush(stderr);
	#endif        
		    start_time = last_req_time;
                   
                  
		//  printf("Reaches here!\n");
                  // KJ: compute the duration of a simulation round
                   if(!record_prev && record_curr)
                   {
                     curr_round_end = last_req_time;
                     disksim->duration = curr_round_end - prev_round_end;
                     fprintf(stderr, "The Detailed Simulation Duration of current workload is %lf\n", disksim->duration);
                     record_curr = 0;
                   } 

                   if(record_prev)
                   {
                     prev_round_end = last_req_time;
                     record_prev = 0;
                   }               

                  // KJ: dumping the stress distribution matrix of last round to a separate file and reinitialize all the elements of matrix to 0;
                  if(disksim->cur_detailed_trace_count > 4)
                  {
                   // char *output_matrix_file = (char*) malloc(sizeof(char) * 500);
                  //  assert(output_file_name);
                    // set file
                    assert(disksim->output_file_prefix);
                    char *output_matrix_file_name = (char*) malloc(sizeof(char) * 500);
                    assert(output_matrix_file_name);
                    sprintf(output_matrix_file_name, "%s_matrix%d.outv", disksim->output_file_prefix, (disksim->cur_detailed_trace_count - 1));
                    FILE *curr_matrix = fopen(output_matrix_file_name, "w");
                    assert(curr_matrix != NULL);
		    // dumping stats to file
                    int row, col;
		    for(row=0; row<disksim->num_row; row++)
                    { 
                      for(col=0; col<disksim->num_col; col++)
                      {
                        fprintf(curr_matrix, "%d ", disksim->stress_dist_matrix[row][col]);
                      }
                      fprintf(curr_matrix, "\n");
                    }		    

                    // reinitialize matrix
                    for(row=0; row<disksim->num_row; row++)
                    {
                      memset(disksim->stress_dist_matrix[row], 0, sizeof(int) * disksim->num_col);
                    }
                    // check if all the elements are 0;
                    for(row=0; row<disksim->num_row; row++)
                      for(col=0; col<disksim->num_col; col++)
                        assert(disksim->stress_dist_matrix[row][col] == 0);
                  }  
                  
		  //printf("Reaches here!\n");

		    if(min_count > max_count) {
	          //device_notify_trace_done();
	      	  addtoextraq((event *) new_ioreq_event);
	      	  return(NULL);
		    }
		      //reset the file pointer to the beginning...
	        /*if(min_count == max_count) {//for the last final pass, change the checkpointing interval to 1hr..
	          fprintf(stderr,"Writing Snapshot at %lf\n",last_req_time);
	          set_snapshot_file();
	          device_notify_trace_done();//dump the state just before the final hourly snapshot.
	          disksim->ssdinfo->ssds[0]->params.checkpoint_time = 1 * 60 * 60 *1000; //1hr.
	          strcat(disksim->snapshot_output_path,"ds_last_snapshot_dir/"); 
	          mode_t temp = umask(0);
	          mkdir(disksim->snapshot_output_path,S_IRWXU | S_IRWXG | S_IRWXO);
	          umask(temp);
	          fprintf(stderr,"Final pass. Snapshot saved in %s\n",disksim->snapshot_output_path);
	          checkpoint_starttime = last_req_time; //Reset the checkpoint start time to current time.
	        }*/
		      rewind(tracefile);
		      goto start_again;
	     }


		if (sscanf(line, "%lf %d %ld %d %x\n", &new_ioreq_event->time, &new_ioreq_event->devno, &new_ioreq_event->blkno, &new_ioreq_event->bcount, &new_ioreq_event->flags) != 5)
		{
			fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
	  		fprintf(stderr, "line: %s", line);
	  		ddbg_assert(0);
	  	}

	   /*#ifdef DEBUG
	    static int request_count = 1;
	    time_t now;
	    now = time(NULL);
	    fprintf(stderr,"Request id:%d,Time:%s\n",request_count++,ctime(&now));
	   #endif*/
	   // by abhishek
	   new_ioreq_event->time = new_ioreq_event->time + start_time;
	   if (new_ioreq_event->flags & ASYNCHRONOUS) {
	      new_ioreq_event->flags |= (new_ioreq_event->flags & READ) ? TIME_LIMITED : 0;
	   } else if (new_ioreq_event->flags & SYNCHRONOUS) {
	      new_ioreq_event->flags |= TIME_CRITICAL;
	   }

	   if(is_first_request == 1) {
	   	checkpoint_starttime = new_ioreq_event->time;
	    disksim_cleanstats();
	   	is_first_request = 0;
	   }

	   new_ioreq_event->buf = 0;
	   new_ioreq_event->opid = 0;
	   new_ioreq_event->busno = 0;
	   new_ioreq_event->cause = 0;


	   realign_request(new_ioreq_event);
	   if(time_to_checkpoint(new_ioreq_event->time)) {
	      fprintf(stderr,"Writing Snapshot at %lf\n",new_ioreq_event->time);
	      set_snapshot_file();
	      device_notify_trace_done();
	      checkpoint_starttime =new_ioreq_event->time;//for the next start
	      disksim_printstats2(sub_status); //dump stats.
	   }

	   last_req_time = new_ioreq_event->time;
	   return(new_ioreq_event);
}


// KJ: This is the single_threshold mode
ioreq_event * iotrace_ascii_get_ioreq_event_fast_forward_1 (FILE* tracefile, ioreq_event *new_ioreq_event)
{
	// static var's lifetime exsits in the whole execution;
	// This function will get called multiple times during execution;
	// static var will get initialized only at compile time;
	static double last_req_time = 0;
	static int is_first_request = 1;
//	static int is_past_threshold = 0;
	static int cur_repeat_trace_count = 0;
	static double start_time_of_next_round = 0;
//	static cur_detailed_trace_count = 0;
	static int past_threshold = 0;	
	double cur_distance = 0;
	char line[201];
	int total_repeat_trace = disksim->repeat_trace;
	int sub_status = 1;
	int rewinded;
	
	
	rewinded = 0;
	// if already end of trace file
	if(fgets(line, 200, tracefile) == NULL)
	{
		cur_repeat_trace_count++;
                disksim->cur_repeat = cur_repeat_trace_count;
		start_time_of_next_round = last_req_time;
		// if total_repeat_trace completes
		if(cur_repeat_trace_count >= total_repeat_trace)
		{
			//end of entire simulation
			addtoextraq((event *) new_ioreq_event);
			return (NULL);
		}
		
		rewind(tracefile);
		rewinded = 1;
                
		// KJ: two consecutive low points will trigger fastforward mode; avoid outliers
		if(disksim->is_past_threshold)
		{
			disksim->cur_detailed_trace_count++;
			fprintf(stderr, "Detailed Simulation Pass %d in Current Round completed Total Pass %d completed!\n", disksim->cur_detailed_trace_count, cur_repeat_trace_count);
			// if detailed duration completes
			if(disksim->cur_detailed_trace_count >= disksim->detailed_simulation_duration)
			{
				// move into fast-forward mode; return the number of trace_repeat rounds in fast-forward mode
				cur_repeat_trace_count += device_fast_forward_simulation(cur_repeat_trace_count);
  				disksim->cur_repeat = cur_repeat_trace_count;
                                fprintf(stderr, "Total rounds of %d has completed\n", cur_repeat_trace_count);
				// if already done; no need to continue
				if(cur_repeat_trace_count >= total_repeat_trace)
				{
					fprintf(stderr, "Total rounds of %d has completed\n", cur_repeat_trace_count);
					//end of entire simulation
					addtoextraq((event *) new_ioreq_event);
					return (NULL);
				}
				// reset the prev_trace_file after fast_forward simulation
				set_prev_trace_file(disksim->curr_record);
				
				// after fast-forward mode, reset the detailed_trace_count to 0 to start another round of detailed
				// simulation
				disksim->cur_detailed_trace_count = 0;
				// after fast forward mode
				start_time_of_next_round = simtime;
			}
		}else
		{
			fprintf(stderr, "Detailed Simulation Pass %d completed\n", cur_repeat_trace_count);
		}

		// avoid outliers
		if(past_threshold == 3)
		{
			disksim->is_past_threshold = 1;
		}
			
	}
	if(rewinded)
	{
		assert(fgets(line, 200, tracefile) != NULL);
                rewinded = 0;
	}
	
	// otherwise, read in the current request from tracefile
	if (sscanf(line, "%lf %d %ld %d %x\n", &new_ioreq_event->time, &new_ioreq_event->devno, &new_ioreq_event->blkno, &new_ioreq_event->bcount, &new_ioreq_event->flags) != 5)
	{
		fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
  		fprintf(stderr, "line: %s", line);
  		ddbg_assert(0);
  	}
	
	// calc the request time of the event
	new_ioreq_event->time = new_ioreq_event->time + start_time_of_next_round;
	last_req_time = new_ioreq_event->time;
	
	if(new_ioreq_event->flags & ASYNCHRONOUS) 
	{
		new_ioreq_event->flags |= (new_ioreq_event->flags & READ) ? TIME_LIMITED : 0;
	}else if(new_ioreq_event->flags & SYNCHRONOUS) {
		new_ioreq_event->flags |= TIME_CRITICAL;
	}
	
	if(is_first_request)
	{
		checkpoint_starttime = new_ioreq_event->time;
		disksim_cleanstats();
		is_first_request = 0;
	}
	
	new_ioreq_event->buf = 0;
	new_ioreq_event->opid = 0;
	new_ioreq_event->busno = 0;
	new_ioreq_event->cause = 0;
	
	// TODO: not entire sure why need to do realign_request
	realign_request(new_ioreq_event);
	
	// if time to checkpoint has come: dump both stat output and snapshot
	if(time_to_checkpoint(new_ioreq_event->time))
	{
		fprintf(stderr, "Detailed Simulation Writing Snapshots at %lf\n", new_ioreq_event->time);
		set_snapshot_file();
		device_notify_trace_done();
		checkpoint_starttime = new_ioreq_event->time;
		// TODO: not sure if it's right, just elminate the dumping of output stats due to unclear seg fault
		disksim_printstats2(sub_status);
		
		// calc the earth mover distance between the most recent 2 snapshots and compare with the threshold to decide whether
		// to enter fast-forward mode or not;
		// TODO: after changing dumping snapshots come back
		if(disksim->snapshot_currindex > 2 && (disksim->is_past_threshold == 0))
		{
			char *curr_snapshot_name = (char*)malloc((strlen(disksim->snapshot_output_path)+150) * sizeof(char));
			char *prev_snapshot_name = (char*)malloc((strlen(disksim->snapshot_output_path)+150) * sizeof(char));
			sprintf(curr_snapshot_name, "%ssnapshot%d.outv", disksim->snapshot_output_path, disksim->snapshot_currindex-1);
		 	sprintf(prev_snapshot_name, "%ssnapshot%d.outv", disksim->snapshot_output_path, disksim->snapshot_currindex-2);
			FILE *curr_snapshot = fopen(curr_snapshot_name, "r");
			FILE *prev_snapshot = fopen(prev_snapshot_name, "r");
			assert(curr_snapshot);
			assert(prev_snapshot);
			// KJ: use opencv to calculate earth mover dist
			cur_distance = device_calculate_distance(prev_snapshot, curr_snapshot);
			// if distance of neighbouring snapshots, enter fast-forward mode; order will be detail->ff->detail->ff ...
                        free(curr_snapshot_name);
                        free(prev_snapshot_name);
			if(cur_distance < disksim->threshold)
			{
				past_threshold ++;
			//	fprintf(stderr, "Acceleration Started!\n");
			//	fprintf(stderr, "##############################################################\n");
			}
		}			
	}
	
	return (new_ioreq_event);	
}



ioreq_event * iotrace_ascii_get_ioreq_event (FILE* tracefile, ioreq_event *new_ioreq_event)
{
	ioreq_event *q = NULL;
        switch(disksim->simulation_strategy)
	{
		case 0:
			q = iotrace_ascii_get_ioreq_event_detailed_simulation(tracefile, new_ioreq_event);
                        break;
			
		case 1:
			q = iotrace_ascii_get_ioreq_event_fast_forward_1(tracefile, new_ioreq_event);
                        break;
	} 
        return q; 
	
}

inline int time_to_checkpoint(double currtime)
{
  return (currtime > (checkpoint_starttime + disksim->ssdinfo->ssds[0]->params.checkpoint_time));
}


void realign_request(ioreq_event *new_ioreq_event)
{

	//align requests to page size;
	int diskblocks_per_page = disksim->ssdinfo->ssds[0]->params.page_size;
	if(new_ioreq_event->blkno%diskblocks_per_page !=0)
		new_ioreq_event->blkno -= (new_ioreq_event->blkno%diskblocks_per_page);
	if(new_ioreq_event->bcount%diskblocks_per_page!=0)
		new_ioreq_event->bcount+= (diskblocks_per_page - new_ioreq_event->bcount%diskblocks_per_page);
	assert(new_ioreq_event->blkno%diskblocks_per_page == 0 && new_ioreq_event->bcount%diskblocks_per_page == 0);

}

static ioreq_event * iotrace_batch_get_ioreq_event (FILE *tracefile, ioreq_event *new_ioreq_event)
{
   char line[201];

   if (fgets(line, 200, tracefile) == NULL) {
      addtoextraq((event *) new_ioreq_event);
      return(NULL);
   }
   if (sscanf(line, "%lf %d %ld %d %x %d\n", &new_ioreq_event->time, &new_ioreq_event->devno, &new_ioreq_event->blkno, &new_ioreq_event->bcount, &new_ioreq_event->flags, &new_ioreq_event->batchno) != 6) {
      fprintf(stderr, "Wrong number of arguments for I/O trace event type\n");
      fprintf(stderr, "line: %s", line);
      ddbg_assert(0);
   }
   if (new_ioreq_event->flags & ASYNCHRONOUS) {
      new_ioreq_event->flags |= (new_ioreq_event->flags & READ) ? TIME_LIMITED : 0;
   } else if (new_ioreq_event->flags & SYNCHRONOUS) {
      new_ioreq_event->flags |= TIME_CRITICAL;
   }
   if (new_ioreq_event->flags & BATCH_COMPLETE) {
     new_ioreq_event->batch_complete = 1;
   } else {
     new_ioreq_event->batch_complete = 0;
   }

   new_ioreq_event->buf = 0;
   new_ioreq_event->opid = 0;
   new_ioreq_event->busno = 0;
   new_ioreq_event->cause = 0;
   return(new_ioreq_event);
}


ioreq_event * iotrace_get_ioreq_event (FILE *tracefile, int traceformat, ioreq_event *temp)
{
   switch (traceformat) {
      
   case ASCII:
      temp = iotrace_ascii_get_ioreq_event(tracefile, temp);
      break;
      
   case RAW:
      temp = iotrace_raw_get_ioreq_event(tracefile, temp);
      break;
      
   case HPL:
      temp = iotrace_hpl_get_ioreq_event(tracefile, temp);
      break;

   case DEC:
      temp = iotrace_dec_get_ioreq_event(tracefile, temp);
      break;

   case VALIDATE:
      temp = iotrace_validate_get_ioreq_event(tracefile, temp);
      break;

   case EMCSYMM:
      temp = iotrace_emcsymm_get_ioreq_event(tracefile, temp);
      break;

   case EMCBACKEND:
      temp = iotrace_emcbackend_get_ioreq_event(tracefile, temp);
      break;

   case BATCH:
      temp = iotrace_batch_get_ioreq_event(tracefile, temp);
      break;

   default:
      fprintf(stderr, "Unknown traceformat in iotrace_get_ioreq_event - %d\n", traceformat);
      exit(1);
   }

   return ((ioreq_event *)temp);
}


static void iotrace_hpl_srt_tracefile_start (char *tracedate)
{
   char crap[40];
   char monthstr[40];
   int day;
   int hour;
   int minute;
   int second;
   int year;

   if (sscanf(tracedate, "%s\t= \"%s %s %d %d:%d:%d %d\";\n", crap, crap, monthstr, &day, &hour, &minute, &second, &year) != 8) {
      fprintf(stderr, "Format problem with 'tracedate' line in HPL trace - %s\n", tracedate);
      exit(1);
   }
   if (baseyear == 0) {
      baseyear = year;
   }
   day = day + iotrace_month_convert(monthstr, year);
   if (year != baseyear) {
      day += (baseyear % 4) ? 365 : 366;
   }
   if (baseday == 0) {
      baseday = day;
   }
   second = second + (60 * minute) + (3600 * hour) + (86400 * (day - baseday));
   if (basesecond == 0) {
      basesecond = second;
   }
   second -= basesecond;
   tracebasetime += (double) 1000 * (double) second;
}


static void iotrace_hpl_initialize_file (FILE *tracefile, int print_tracefile_header)
{
   char letter = '0';
   char line[201];
   char linetype[40];

   if (disksim->traceheader == FALSE) {
      return;
   }
   while (1) {
      if (fgets(line, 200, tracefile) == NULL) {
         fprintf(stderr, "No 'tracedate' line in HPL trace\n");
         exit(1);
      }
      sscanf(line, "%s", linetype);
      if (strcmp(linetype, "tracedate") == 0) {
         break;
      }
   }
   iotrace_hpl_srt_tracefile_start(line);
   while (letter != 0x0C) {
      if (fscanf(tracefile, "%c", &letter) != 1) {
         fprintf(stderr, "End of header information never found - end of file\n");
         exit(1);
      }
      if ((print_tracefile_header) && (letter != 0x0C)) {
         printf("%c", letter);
      }
   }
}


void iotrace_initialize_file (FILE *tracefile, int traceformat, int print_tracefile_header)
{
   if (traceformat == HPL) {
      iotrace_hpl_initialize_file(tracefile, print_tracefile_header);
   }
}


void iotrace_printstats (FILE *outfile)
{
   if (disksim->iotrace_info == NULL) {
      return;
   }
   if (hpreads | hpwrites) {
      fprintf (outfile, "\n");
      fprintf(outfile, "Total reads:    \t%d\t%5.2f\n", hpreads, ((double) hpreads / (double) (hpreads + hpwrites)));
      fprintf(outfile, "Total writes:   \t%d\t%5.2f\n", hpwrites, ((double) hpwrites / (double) (hpreads + hpwrites)));
      fprintf(outfile, "Sync Reads:  \t%d\t%5.2f\t%5.2f\n", syncreads, ((double) syncreads / (double) (hpreads + hpwrites)), ((double) syncreads / (double) hpreads));
      fprintf(outfile, "Sync Writes: \t%d\t%5.2f\t%5.2f\n", syncwrites, ((double) syncwrites / (double) (hpreads + hpwrites)), ((double) syncwrites / (double) hpwrites));
      fprintf(outfile, "Async Reads: \t%d\t%5.2f\t%5.2f\n", asyncreads, ((double) asyncreads / (double) (hpreads + hpwrites)), ((double) asyncreads / (double) hpreads));
      fprintf(outfile, "Async Writes:\t%d\t%5.2f\t%5.2f\n", asyncwrites, ((double) asyncwrites / (double) (hpreads + hpwrites)), ((double) asyncwrites / (double) hpwrites));
   }
}

