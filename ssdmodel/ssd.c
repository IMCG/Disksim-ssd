// DiskSim SSD support
// �2008 Microsoft Corporation. All Rights Reserved

#include "ssd.h"
#include "ssd_timing.h"
#include "ssd_clean.h"
#include "ssd_gang.h"
#include "ssd_init.h"
#include "ssd_utils.h"
#include "modules/ssdmodel_ssd_param.h"
#include "ssd_refresh.h"
// KJ: to compute distance and histogram
#include "cv.h"
#include "highgui.h"
#include <math.h>

#ifndef sprintf_s
#define sprintf_s3(x,y,z) sprintf(x,z)
#define sprintf_s4(x,y,z,w) sprintf(x,z,w)
#define sprintf_s5(x,y,z,w,s) sprintf(x,z,w,s)
#else
#define sprintf_s3(x,y,z) sprintf_s(x,z)
#define sprintf_s4(x,y,z,w) sprintf_s(x,z,w)
#define sprintf_s5(x,y,z,w,s) sprintf_s(x,z,w,s)
#endif

#ifndef _strdup
#define _strdup strdup
#endif

static void ssd_request_complete(ioreq_event *curr);
static void ssd_media_access_request(ioreq_event *curr);
// KJ: extern these 3 methods here; from disksim_iotrace.c and disksim.c and disksim.c
int time_to_checkpoint(double currtime);
void set_snapshot_file();
void disksim_printstats2(int sub_status);

struct ssd *getssd (int devno)
{
   struct ssd *s;
   ASSERT1((devno >= 0) && (devno < MAXDEVICES), "devno", devno);

   s = disksim->ssdinfo->ssds[devno];
   return (disksim->ssdinfo->ssds[devno]);
}

int ssd_set_depth (int devno, int inbusno, int depth, int slotno)
{
   ssd_t *currdisk;
   int cnt;

   currdisk = getssd (devno);
   assert(currdisk);
   cnt = currdisk->numinbuses;
   currdisk->numinbuses++;
   if ((cnt + 1) > MAXINBUSES) {
      fprintf(stderr, "Too many inbuses specified for ssd %d - %d\n", devno, (cnt+1));
      exit(1);
   }
   currdisk->inbuses[cnt] = inbusno;
   currdisk->depth[cnt] = depth;
   currdisk->slotno[cnt] = slotno;
   return(0);
}

int ssd_get_depth (int devno)
{
   ssd_t *currdisk;
   currdisk = getssd (devno);
   return(currdisk->depth[0]);
}

int ssd_get_slotno (int devno)
{
   ssd_t *currdisk;
   currdisk = getssd (devno);
   return(currdisk->slotno[0]);
}

int ssd_get_inbus (int devno)
{
   ssd_t *currdisk;
   currdisk = getssd (devno);
   return(currdisk->inbuses[0]);
}

int ssd_get_maxoutstanding (int devno)
{
   ssd_t *currdisk;
   currdisk = getssd (devno);
   return(currdisk->maxqlen);
}

double ssd_get_blktranstime (ioreq_event *curr)
{
   ssd_t *currdisk;
   double tmptime;

   currdisk = getssd (curr->devno);
   tmptime = bus_get_transfer_time(ssd_get_busno(curr), 1, (curr->flags & READ));
   if (tmptime < currdisk->blktranstime) {
      tmptime = currdisk->blktranstime;
   }
   return(tmptime);
}

int ssd_get_busno (ioreq_event *curr)
{
   ssd_t *currdisk;
   intchar busno;
   int depth;

   currdisk = getssd (curr->devno);
   busno.value = curr->busno;
   depth = currdisk->depth[0];
   return(busno.byte[depth]);
}

static void ssd_assert_current_activity(ssd_t *currdisk, ioreq_event *curr)
{
    assert(currdisk->channel_activity != NULL &&
        currdisk->channel_activity->devno == curr->devno &&
        currdisk->channel_activity->blkno == curr->blkno &&
        currdisk->channel_activity->bcount == curr->bcount);
}

/*
 * ssd_send_event_up_path()
 *
 * Acquires the bus (if not already acquired), then uses bus_delay to
 * send the event up the path.
 *
 * If the bus is already owned by this device or can be acquired
 * immediately (interleaved bus), the event is sent immediately.
 * Otherwise, ssd_bus_ownership_grant will later send the event.
 */
static void ssd_send_event_up_path (ioreq_event *curr, double delay)
{
   ssd_t *currdisk;
   int busno;
   int slotno;

   // fprintf (outputfile, "ssd_send_event_up_path - devno %d, type %d, cause %d, blkno %d\n", curr->devno, curr->type, curr->cause, curr->blkno);

   currdisk = getssd (curr->devno);

   ssd_assert_current_activity(currdisk, curr);

   busno = ssd_get_busno(curr);
   slotno = currdisk->slotno[0];

   /* Put new request at head of buswait queue */
   curr->next = currdisk->buswait;
   currdisk->buswait = curr;

   curr->tempint1 = busno;
   curr->time = delay;
   if (currdisk->busowned == -1) {

      // fprintf (outputfile, "Must get ownership of the bus first\n");

      if (curr->next) {
         //fprintf(stderr,"Multiple bus requestors detected in ssd_send_event_up_path\n");
         /* This should be ok -- counting on the bus module to sequence 'em */
      }
      if (bus_ownership_get(busno, slotno, curr) == FALSE) {
         /* Remember when we started waiting (only place this is written) */
         currdisk->stat.requestedbus = simtime;
      } else {
         currdisk->busowned = busno;
         bus_delay(busno, DEVICE, curr->devno, delay, curr); /* Never for SCSI */
      }
   } else if (currdisk->busowned == busno) {

      //fprintf (outputfile, "Already own bus - so send it on up\n");

      bus_delay(busno, DEVICE, curr->devno, delay, curr);
   } else {
      fprintf(stderr, "Wrong bus owned for transfer desired\n");
      exit(1);
   }
}

/* The idea here is that only one request can "possess" the channel back to the
   controller at a time. All others are enqueued on queue of pending activities.
   "Completions" ... those operations that need only be signaled as done to the
   controller ... are given on this queue.  The "channel_activity" field indicates
   whether any operation currently possesses the channel.

   It is our hope that new requests cannot enter the system when the channel is
   possessed by another operation.  This would not model reality!!  However, this
   code (and that in ssd_request_arrive) will handle this case "properly" by enqueuing
   the incoming request.  */

static void ssd_check_channel_activity (ssd_t *currdisk)
{
   while (1) {
       ioreq_event *curr = currdisk->completion_queue;
       currdisk->channel_activity = curr;
       if (curr != NULL) {
           currdisk->completion_queue = curr->next;
           if (currdisk->neverdisconnect) {
               /* already connected */
               if (curr->flags & READ) {
                   /* transfer data up the line: curr->bcount, which is still set to */
                   /* original requested value, indicates how many blks to transfer. */
                   curr->type = DEVICE_DATA_TRANSFER_COMPLETE;
                   ssd_send_event_up_path(curr, (double) 0.0);
               } else {
                   ssd_request_complete (curr);
               }
           } else {
               /* reconnect to controller */
               curr->type = IO_INTERRUPT_ARRIVE;
               curr->cause = RECONNECT;
               ssd_send_event_up_path (curr, currdisk->bus_transaction_latency);
               currdisk->reconnect_reason = DEVICE_ACCESS_COMPLETE;
           }
       } else {
           curr = ioqueue_get_next_request(currdisk->queue);
           currdisk->channel_activity = curr;
           if (curr != NULL) {
               if (curr->flags & READ) {
                   if (!currdisk->neverdisconnect) {
                       ssd_media_access_request(ioreq_copy(curr));
                       curr->type = IO_INTERRUPT_ARRIVE;
                       curr->cause = DISCONNECT;
                       ssd_send_event_up_path (curr, currdisk->bus_transaction_latency);
                   } else {
                       ssd_media_access_request(curr);
                       continue;
                   }
               } else {
                   curr->cause = RECONNECT;
                   curr->type = IO_INTERRUPT_ARRIVE;
                   currdisk->reconnect_reason = IO_INTERRUPT_ARRIVE;
                   ssd_send_event_up_path (curr, currdisk->bus_transaction_latency);
               }
           }
       }
       break;
   }
}

/*
 * ssd_bus_ownership_grant
 * Calls bus_delay to handle the event that the disk has been granted the bus.  I believe
 * this is always initiated by a call to ssd_send_even_up_path.
 */
void ssd_bus_ownership_grant (int devno, ioreq_event *curr, int busno, double arbdelay)
{
   ssd_t *currdisk;
   ioreq_event *tmp;

   currdisk = getssd (devno);

   ssd_assert_current_activity(currdisk, curr);
   tmp = currdisk->buswait;
   while ((tmp != NULL) && (tmp != curr)) {
      tmp = tmp->next;
   }
   if (tmp == NULL) {
      fprintf(stderr, "Bus ownership granted to unknown ssd request - devno %d, busno %d\n", devno, busno);
      exit(1);
   }
   currdisk->busowned = busno;
   currdisk->stat.waitingforbus += arbdelay;
   //ASSERT (arbdelay == (simtime - currdisk->stat.requestedbus));
   currdisk->stat.numbuswaits++;
   bus_delay(busno, DEVICE, devno, tmp->time, tmp);
}

void ssd_bus_delay_complete (int devno, ioreq_event *curr, int sentbusno)
{
   ssd_t *currdisk;
   intchar slotno;
   intchar busno;
   int depth;

   currdisk = getssd (devno);
   ssd_assert_current_activity(currdisk, curr);

   // fprintf (outputfile, "Entered ssd_bus_delay_complete\n");

   // EPW: I think the buswait logic doesn't do anything, is confusing, and risks
   // overusing the "next" field, although an item shouldn't currently be a queue.
   if (curr == currdisk->buswait) {
      currdisk->buswait = curr->next;
   } else {
      ioreq_event *tmp = currdisk->buswait;
      while ((tmp->next != NULL) && (tmp->next != curr)) {
         tmp = tmp->next;
      }
      if (tmp->next != curr) {
          // fixed a warning here
          //fprintf(stderr, "Bus delay complete for unknown ssd request - devno %d, busno %d\n", devno, busno.value);
          fprintf(stderr, "Bus delay complete for unknown ssd request - devno %d, busno %d\n", devno, curr->busno);
         exit(1);
      }
      tmp->next = curr->next;
   }
   busno.value = curr->busno;
   slotno.value = curr->slotno;
   depth = currdisk->depth[0];
   slotno.byte[depth] = slotno.byte[depth] >> 4;
   curr->time = 0.0;
   if (depth == 0) {
      intr_request ((event *)curr);
   } else {
      bus_deliver_event(busno.byte[depth], slotno.byte[depth], curr);
   }
}


/*
 * send completion up the line
 */
static void ssd_request_complete(ioreq_event *curr)
{
   ssd_t *currdisk;
   ioreq_event *x;

   // fprintf (outputfile, "Entering ssd_request_complete: %12.6f\n", simtime);

   currdisk = getssd (curr->devno);
   ssd_assert_current_activity(currdisk, curr);

   if ((x = ioqueue_physical_access_done(currdisk->queue,curr)) == NULL) {
      fprintf(stderr, "ssd_request_complete:  ioreq_event not found by ioqueue_physical_access_done call\n");
      exit(1);
   }

   /* send completion interrupt */
   curr->type = IO_INTERRUPT_ARRIVE;
   curr->cause = COMPLETION;
   ssd_send_event_up_path(curr, currdisk->bus_transaction_latency);
}

static void ssd_bustransfer_complete (ioreq_event *curr)
{
   // fprintf (outputfile, "Entering ssd_bustransfer_complete for disk %d: %12.6f\n", curr->devno, simtime);

   if (curr->flags & READ) {
      ssd_request_complete (curr);
   } else {
      ssd_t *currdisk = getssd (curr->devno);
      ssd_assert_current_activity(currdisk, curr);
      if (currdisk->neverdisconnect == FALSE) {
          /* disconnect from bus */
          ioreq_event *tmp = ioreq_copy (curr);
          tmp->type = IO_INTERRUPT_ARRIVE;
          tmp->cause = DISCONNECT;
          ssd_send_event_up_path (tmp, currdisk->bus_transaction_latency);
          ssd_media_access_request (curr);
      } else {
          ssd_media_access_request (curr);
          ssd_check_channel_activity (currdisk);
      }
   }
}

/*
 * returns the logical page number within an element given a block number as
 * issued by the file system
 */
int ssd_logical_pageno(int blkno, ssd_t *s)
{
    int apn;
    int lpn;

    // absolute page number is the block number as written by the above layer
    apn = blkno/s->params.page_size;

    // find the logical page number within the ssd element. we maintain the
    // mapping between the logical page number and the actual physical page
    // number. an alternative is that we could maintain the mapping between
    // apn we calculated above and the physical page number. but the range
    // of apn is several times bigger and so we chose to go with the mapping
    // b/w lpn --> physical page number
    lpn = ((apn - (apn % (s->params.element_stride_pages * s->params.nelements)))/
                      s->params.nelements) + (apn % s->params.element_stride_pages);

    return lpn;
}

int ssd_already_present(ssd_req **reqs, int total, ioreq_event *req)
{
    int i;
    int found = 0;

    for (i = 0; i < total; i ++) {
        if ((req->blkno == reqs[i]->org_req->blkno) &&
            (req->flags == reqs[i]->org_req->flags)) {
            found = 1;
            break;
        }
    }

    return found;
}

double _ssd_invoke_element_cleaning(int elem_num, ssd_t *s)
{
    double clean_cost = ssd_clean_element(s, elem_num);
    return clean_cost;
}

static int ssd_invoke_element_cleaning(int elem_num, ssd_t *s)
{
    double max_cost = 0;
    int cleaning_invoked = 0;
    ssd_element *elem = &s->elements[elem_num];

    // element must be free
    ASSERT(elem->media_busy == FALSE);

    max_cost = _ssd_invoke_element_cleaning(elem_num, s);
    //printf("mxcost %d\n",max_cost);  
    // cleaning was invoked on this element. we can start
    // the next operation on this elem only after the cleaning
    // gets over.
    if (max_cost > 0) {
        ioreq_event *tmp;

        elem->media_busy = TRUE;
        cleaning_invoked = 1;

        // we use the 'blkno' field to store the element number
        tmp = (ioreq_event *)getfromextraq();
        tmp->devno = s->devno;
        tmp->time = simtime + max_cost;
        tmp->blkno = elem_num;
        tmp->ssd_elem_num = elem_num;
        tmp->type = SSD_CLEAN_ELEMENT;
        tmp->flags = SSD_CLEAN_ELEMENT;
        tmp->busno = -1;
        tmp->bcount = -1;
        stat_update (&s->stat.acctimestats, max_cost);
        addtointq ((event *)tmp);
#ifdef DEBUG
        //fprintf(stderr,"Cleaning invoked @%f on element %d\n",simtime,elem_num);
        //fflush(stderr);
#endif
        // stat
        elem->stat.tot_clean_time += max_cost;
    }

    return cleaning_invoked;
}

static void ssd_activate_elem(ssd_t *currdisk, int elem_num)
{
    ioreq_event *req;
    ssd_req **read_reqs;
    ssd_req **write_reqs;
    int i;
    int read_total = 0;
    int write_total = 0;
    double schtime = 0;
    int max_reqs;
    int tot_reqs_issued;
    double max_time_taken = 0;


    ssd_element *elem = &currdisk->elements[elem_num];

#ifdef DEBUG
    static int debug_print = 0;
    if (debug_print == 1) {
      if (elem_num == 7) {
        fprintf(stderr,"Media busy:%d,simtime:%lf\n",elem->media_busy,simtime);
      }
    }
#endif
 
    // if the media is busy, we can't do anything, so return
    if (elem->media_busy == TRUE) {
        return;
    }

    ASSERT(ioqueue_get_reqoutstanding(elem->queue) == 0);

    // we can invoke cleaning in the background whether there
    // is request waiting or not
    if (currdisk->params.cleaning_in_background) {
        // if cleaning was invoked, wait until
        // it is over ...
        if (ssd_invoke_element_cleaning(elem_num, currdisk)) {
            return;
        }
    }

    ASSERT(elem->metadata.reqs_waiting == ioqueue_get_number_in_queue(elem->queue));

    if (elem->metadata.reqs_waiting > 0) {

        // invoke cleaning in foreground when there are requests waiting
        if (!currdisk->params.cleaning_in_background) {
            // if cleaning was invoked, wait until
            // it is over ...
            if (ssd_invoke_element_cleaning(elem_num, currdisk)) {
                return;
            }
        }

        // how many reqs can we issue at once
        if (currdisk->params.copy_back == SSD_COPY_BACK_DISABLE) {
            max_reqs = 1;
        } else {
            if (currdisk->params.num_parunits == 1) {
                max_reqs = 1;
            } else {
                max_reqs = MAX_REQS_ELEM_QUEUE;
            }
        }

        // ideally, we should issue one req per plane, overlapping them all.
        // in order to simplify the overlapping strategy, let's issue
        // requests of the same type together.

        read_reqs = (ssd_req **) malloc(max_reqs * sizeof(ssd_req *));
        write_reqs = (ssd_req **) malloc(max_reqs * sizeof(ssd_req *));

        // collect the requests
        while ((req = ioqueue_get_next_request(elem->queue)) != NULL) {
            int found = 0;

            elem->metadata.reqs_waiting --;

            // see if we already have the same request in the list.
            // this usually doesn't happen -- but on synthetic traces
            // this weird case can occur.







            if (req->flags & READ) {
                found = ssd_already_present(read_reqs, read_total, req);
            } else {
                found = ssd_already_present(write_reqs, write_total, req);
            }

            if (!found) {
                // this is a valid request
                ssd_req *r = (ssd_req*) malloc(sizeof(ssd_req));
                r->blk = req->blkno;
                r->count = req->bcount;
                r->is_read = req->flags & READ;
                r->org_req = req;
                r->plane_num = -1; // we don't know to which plane this req will be directed at

                if (req->flags & READ) {
                    read_reqs[read_total] = r;
                    read_total ++;
                } else {
                    write_reqs[write_total] = r;
                    write_total ++;
                }

                // if we have more reqs than we can handle, quit
                if ((read_total >= max_reqs) ||
                    (write_total >= max_reqs)) {
                    break;
                }
            } else {
                // throw this request -- it doesn't make sense
                stat_update (&currdisk->stat.acctimestats, 0);
                req->time = simtime;
                req->ssd_elem_num = elem_num;
                req->type = DEVICE_ACCESS_COMPLETE;
                addtointq ((event *)req);
            }
        }

        if (read_total > 0) {
            // first issue all the read requests (it doesn't matter what we
            // issue first). i chose read because reads are mostly synchronous.
            // find the time taken to serve these requests.
            ssd_compute_access_time(currdisk, elem_num, read_reqs, read_total);

            // add an event for each request completion
            for (i = 0; i < read_total; i ++) {
              elem->media_busy = TRUE;

              // find the maximum time taken by a request
              if (schtime < read_reqs[i]->schtime) {
                  schtime = read_reqs[i]->schtime;
              }

              stat_update (&currdisk->stat.acctimestats, read_reqs[i]->acctime);
              read_reqs[i]->org_req->time = simtime + read_reqs[i]->schtime;
              read_reqs[i]->org_req->ssd_elem_num = elem_num;
              read_reqs[i]->org_req->type = DEVICE_ACCESS_COMPLETE;

              //printf("R: blk %d elem %d acctime %f simtime %f\n", read_reqs[i]->blk,
                //  elem_num, read_reqs[i]->acctime, read_reqs[i]->org_req->time);

              addtointq ((event *)read_reqs[i]->org_req);
              free(read_reqs[i]);
            }
        }

        free(read_reqs);

        max_time_taken = schtime;

        if (write_total > 0) {
            // next issue the write requests
            ssd_compute_access_time(currdisk, elem_num, write_reqs, write_total);

            // add an event for each request completion.
            // note that we can issue the writes only after all the reads above are
            // over. so, include the maximum read time when creating the event.
            for (i = 0; i < write_total; i ++) {
              elem->media_busy = TRUE;

              stat_update (&currdisk->stat.acctimestats, write_reqs[i]->acctime);
              write_reqs[i]->org_req->time = simtime + schtime + write_reqs[i]->schtime;
              //printf("blk %d elem %d acc time %f\n", write_reqs[i]->blk, elem_num, write_reqs[i]->acctime);

              if (max_time_taken < (schtime+write_reqs[i]->schtime)) {
                  max_time_taken = (schtime+write_reqs[i]->schtime);
              }

              write_reqs[i]->org_req->ssd_elem_num = elem_num;
              write_reqs[i]->org_req->type = DEVICE_ACCESS_COMPLETE;
              //printf("W: blk %d elem %d acctime %f simtime %f\n", write_reqs[i]->blk,
                //  elem_num, write_reqs[i]->acctime, write_reqs[i]->org_req->time);

              addtointq ((event *)write_reqs[i]->org_req);
              free(write_reqs[i]);
            }
        }

        free(write_reqs);

        // statistics
        tot_reqs_issued = read_total + write_total;
        ASSERT(tot_reqs_issued > 0);
        currdisk->elements[elem_num].stat.tot_reqs_issued += tot_reqs_issued;
        currdisk->elements[elem_num].stat.tot_time_taken += max_time_taken;
    }
}


static void ssd_do_refresh(ssd_t *currdisk,double now)
{
  if(is_queue_empty(currdisk->refresh_queue))
    return; //Queue is empty. No refresh to do.

  
  switch(currdisk->params.refresh_policy) {

    case SSD_FCFS_REFRESH:
    case SSD_EDF_REFRESH:
      fprintf(stderr, "SSD refresh is involved\n");
      ssd_perform_refresh(currdisk,now);
      break;
    case SSD_REFRESH_NONE:
      fprintf(stderr, "SSD_REFRESH_NONE\n");
      return;
    default:
      fprintf(stderr,"Invalid refresh policy specified:%d\n",currdisk->params.refresh_policy);
      break;
  }
}

//Is it time to do a refresh? Return yes or no.
int should_i_refresh_now(double now, double next_refresh)
{
 return (now > next_refresh);
}

static void ssd_media_access_request_element(ioreq_event *curr)
{
   ssd_t *currdisk = getssd(curr->devno);
   // TODO: deactivate refresh
   
   if(should_i_refresh_now(simtime, currdisk->next_refresh_time)) {  //Is it time to do refresh yet?
   //   fprintf(stderr, "refresh started\n");
      currdisk->next_refresh_time = simtime + currdisk->params.refresh_interval;
      ssd_do_refresh(currdisk,curr->time);
   }
   
   int blkno = curr->blkno;
   int count = curr->bcount;
#ifdef DEBUG
   int prev_block = -1;
#endif
   /* **** CAREFUL ... HIJACKING tempint2 and tempptr2 fields here **** */
   curr->tempint2 = count;
   while (count != 0) {

       // find the element (package) to direct the request
       int elem_num = currdisk->timing_t->choose_element(currdisk->timing_t, blkno);
       ssd_element *elem = &currdisk->elements[elem_num];

       // create a new sub-request for the element
       ioreq_event *tmp = (ioreq_event *)getfromextraq();
       memcpy((void*)tmp,(void*)curr,sizeof(curr));
       tmp->devno = curr->devno;
       tmp->busno = curr->busno;
       tmp->flags = curr->flags;
       tmp->blkno = blkno;
       tmp->opid = curr->opid;
       tmp->bcount = ssd_choose_aligned_count(currdisk->params.page_size, blkno, count);
       ASSERT(tmp->bcount == currdisk->params.page_size);

       tmp->tempptr2 = curr;
#ifdef DEBUG
       if(prev_block != -1) {
          ASSERT(tmp->blkno - prev_block == currdisk->params.page_size);
       }
       prev_block = tmp->blkno;
#endif
       blkno += tmp->bcount;
       count -= tmp->bcount;

       elem->metadata.reqs_waiting ++;

       // add the request to the corresponding element's queue
       ioqueue_add_new_request(elem->queue, (ioreq_event *)tmp);
       ssd_activate_elem(currdisk, elem_num);
   }
}

static void ssd_media_access_request (ioreq_event *curr)
{
    ssd_t *currdisk = getssd(curr->devno);

    switch(currdisk->params.alloc_pool_logic) {
        case SSD_ALLOC_POOL_PLANE:
        case SSD_ALLOC_POOL_CHIP:
            ssd_media_access_request_element(curr);
        break;

        case SSD_ALLOC_POOL_GANG:
#if SYNC_GANG
            ssd_media_access_request_gang_sync(curr);
#else
            ssd_media_access_request_gang(curr);
#endif
        break;

        default:
            printf("Unknown alloc pool logic %d\n", currdisk->params.alloc_pool_logic);
            ASSERT(0);
            break;
    }
}

static void ssd_reconnect_done (ioreq_event *curr)
{
   ssd_t *currdisk;

   // fprintf (outputfile, "Entering ssd_reconnect_done for disk %d: %12.6f\n", curr->devno, simtime);

   currdisk = getssd (curr->devno);
   ssd_assert_current_activity(currdisk, curr);

   if (curr->flags & READ) {
      if (currdisk->neverdisconnect) {
         /* Just holding on to bus; data transfer will be initiated when */
         /* media access is complete.                                    */
         addtoextraq((event *) curr);
         ssd_check_channel_activity (currdisk);
      } else {
         /* data transfer: curr->bcount, which is still set to original */
         /* requested value, indicates how many blks to transfer.       */
         curr->type = DEVICE_DATA_TRANSFER_COMPLETE;
         ssd_send_event_up_path(curr, (double) 0.0);
      }

   } else {
      if (currdisk->reconnect_reason == DEVICE_ACCESS_COMPLETE) {
         ssd_request_complete (curr);

      } else {
         /* data transfer: curr->bcount, which is still set to original */
         /* requested value, indicates how many blks to transfer.       */
         curr->type = DEVICE_DATA_TRANSFER_COMPLETE;
         ssd_send_event_up_path(curr, (double) 0.0);
      }
   }
}

static void ssd_request_arrive (ioreq_event *curr)
{
   ssd_t *currdisk;

   // fprintf (outputfile, "Entering ssd_request_arrive: %12.6f\n", simtime);
   // fprintf (outputfile, "ssd = %d, blkno = %d, bcount = %d, read = %d\n",curr->devno, curr->blkno, curr->bcount, (READ & curr->flags));

   currdisk = getssd(curr->devno);

   /* verify that request is valid. */
   if ((curr->blkno < 0) || (curr->bcount <= 0) ||
       ((curr->blkno + curr->bcount) > currdisk->numblocks)) {
      fprintf(stderr, "Invalid set of blocks requested from ssd - blkno %ld, bcount %d, numblocks %d\n", curr->blkno, curr->bcount, currdisk->numblocks);
      exit(1);
   }

   /* create a new request, set it up for initial interrupt */
   ioqueue_add_new_request(currdisk->queue, curr);
   if (currdisk->channel_activity == NULL) {

      curr = ioqueue_get_next_request(currdisk->queue);
      currdisk->busowned = ssd_get_busno(curr);
      currdisk->channel_activity = curr;
      currdisk->reconnect_reason = IO_INTERRUPT_ARRIVE;

      if (curr->flags & READ) {
          if (!currdisk->neverdisconnect) {
              ssd_media_access_request (ioreq_copy(curr));
              curr->cause = DISCONNECT;
              curr->type = IO_INTERRUPT_ARRIVE;
              ssd_send_event_up_path(curr, currdisk->bus_transaction_latency);
          } else {
              ssd_media_access_request (curr);
              ssd_check_channel_activity(currdisk);
          }
      } else {
         curr->cause = READY_TO_TRANSFER;
         curr->type = IO_INTERRUPT_ARRIVE;
         ssd_send_event_up_path(curr, currdisk->bus_transaction_latency);
      }
   }
}

/*
 * cleaning in an element is over.
 */
static void ssd_clean_element_complete(ioreq_event *curr)
{
   ssd_t *currdisk;
   int elem_num;

   currdisk = getssd (curr->devno);
   elem_num = curr->ssd_elem_num;
   ASSERT(currdisk->elements[elem_num].media_busy == TRUE);

   // release this event
   addtoextraq((event *) curr);

   // activate the gang to serve the next set of requests
   currdisk->elements[elem_num].media_busy = FALSE;
   ssd_activate_elem(currdisk, elem_num);
}

static void ssd_refresh_element_complete(ioreq_event *curr)
{
   ssd_t *currdisk;
   int elem_num;

   currdisk = getssd (curr->devno);
   elem_num = curr->ssd_elem_num;
   ASSERT(currdisk->elements[elem_num].media_busy == TRUE);

   // release this event
   addtoextraq((event *) curr);

   // activate the gang to serve the next set of requests
   currdisk->elements[elem_num].media_busy = FALSE;
   ssd_activate_elem(currdisk, elem_num);
}


void ssd_complete_parent(ioreq_event *curr, ssd_t *currdisk)
{
    ioreq_event *parent;

    /* **** CAREFUL ... HIJACKING tempint2 and tempptr2 fields here **** */
    parent = (ioreq_event*) curr->tempptr2;
    parent->tempint2 -= curr->bcount;

    if (parent->tempint2 == 0) {
      ioreq_event *prev;

      assert(parent != currdisk->channel_activity);
      prev = currdisk->completion_queue;
      if (prev == NULL) {
         currdisk->completion_queue = parent;
         parent->next = prev;
      } else {
         while (prev->next != NULL)
            prev = prev->next;
            parent->next = prev->next;
            prev->next = parent;
      }
      if (currdisk->channel_activity == NULL) {
         ssd_check_channel_activity (currdisk);
      }
    }
}

static void ssd_access_complete_element(ioreq_event *curr)
{
   ssd_t *currdisk;
   int elem_num;
   ssd_element  *elem;
   ioreq_event *x;


   #ifdef DEBUG
   static int debug_print = 0;
   if(debug_print == 1) {
        fprintf(stderr,"%lf, %ld,%ld\n",curr->time,curr->blkno,curr->opid);
   }
   #endif //DEBUG
   currdisk = getssd (curr->devno);
   elem_num = currdisk->timing_t->choose_element(currdisk->timing_t, curr->blkno);
   ASSERT(elem_num == curr->ssd_elem_num);
   elem = &currdisk->elements[elem_num];

   if ((x = ioqueue_physical_access_done(elem->queue,curr)) == NULL) {
      fprintf(stderr, "ssd_access_complete:  ioreq_event not found by ioqueue_physical_access_done call\n");
      exit(1);
   }

   // all the reqs are over
   if (ioqueue_get_reqoutstanding(elem->queue) == 0) {
    elem->media_busy = FALSE;
   }

   ssd_complete_parent(curr, currdisk);
   addtoextraq((event *) curr);
   ssd_activate_elem(currdisk, elem_num);
}

static void ssd_access_complete(ioreq_event *curr)
{
    ssd_t *currdisk = getssd (curr->devno);;

    switch(currdisk->params.alloc_pool_logic) {
        case SSD_ALLOC_POOL_PLANE:
        case SSD_ALLOC_POOL_CHIP:
            ssd_access_complete_element(curr);
        break;

        case SSD_ALLOC_POOL_GANG:
#if SYNC_GANG
            ssd_access_complete_gang_sync(curr);
#else
            ssd_access_complete_gang(curr);
#endif
        break;

        default:
            printf("Unknown alloc pool logic %d\n", currdisk->params.alloc_pool_logic);
            ASSERT(0);
            break;
    }
}

/* intermediate disconnect done */
static void ssd_disconnect_done (ioreq_event *curr)
{
   ssd_t *currdisk;

   currdisk = getssd (curr->devno);
   ssd_assert_current_activity(currdisk, curr);

   // fprintf (outputfile, "Entering ssd_disconnect for disk %d: %12.6f\n", currdisk->devno, simtime);

   addtoextraq((event *) curr);

   if (currdisk->busowned != -1) {
      bus_ownership_release(currdisk->busowned);
      currdisk->busowned = -1;
   }
   ssd_check_channel_activity (currdisk);
}

/* completion disconnect done */
static void ssd_completion_done (ioreq_event *curr)
{
   ssd_t *currdisk = getssd (curr->devno);
   ssd_assert_current_activity(currdisk, curr);

   // fprintf (outputfile, "Entering ssd_completion for disk %d: %12.6f\n", currdisk->devno, simtime);

   addtoextraq((event *) curr);
#ifdef DEBUG
   static int num_request = 0;
   num_request++;
#endif
   if (currdisk->busowned != -1) {
      bus_ownership_release(currdisk->busowned);
      currdisk->busowned = -1;
   }

   ssd_check_channel_activity (currdisk);
}

static void ssd_interrupt_complete (ioreq_event *curr)
{
   // fprintf (outputfile, "Entered ssd_interrupt_complete - cause %d\n", curr->cause);

   switch (curr->cause) {

      case RECONNECT:
         ssd_reconnect_done(curr);
     break;

      case DISCONNECT:
     ssd_disconnect_done(curr);
     break;

      case COMPLETION:
     ssd_completion_done(curr);
     break;

      default:
          assert(0);
          break;
         //ddbg_assert2(0, (const char*)"bad event type");
   }
}


void ssd_event_arrive (ioreq_event *curr)
{
   ssd_t *currdisk;

   // fprintf (outputfile, "Entered ssd_event_arrive: time %f (simtime %f)\n", curr->time, simtime);
   // fprintf (outputfile, " - devno %d, blkno %d, type %d, cause %d, read = %d\n", curr->devno, curr->blkno, curr->type, curr->cause, curr->flags & READ);

   currdisk = getssd (curr->devno);

   switch (curr->type) {

      case IO_ACCESS_ARRIVE:
         curr->time = simtime + currdisk->overhead;
         curr->type = DEVICE_OVERHEAD_COMPLETE;
         addtointq((event *) curr);
         break;

      case DEVICE_OVERHEAD_COMPLETE:
         ssd_request_arrive(curr);
         break;

      case DEVICE_ACCESS_COMPLETE:
         ssd_access_complete (curr);
         break;

      case DEVICE_DATA_TRANSFER_COMPLETE:
         ssd_bustransfer_complete(curr);
         break;

      case IO_INTERRUPT_COMPLETE:
         ssd_interrupt_complete(curr);
         break;

      case IO_QLEN_MAXCHECK:
         /* Used only at initialization time to set up queue stuff */
         curr->tempint1 = -1;
         curr->tempint2 = ssd_get_maxoutstanding(curr->devno);
         curr->bcount = 0;
         break;

      case SSD_CLEAN_GANG:
          ssd_clean_gang_complete(curr);
          break;

      case SSD_CLEAN_ELEMENT:
          ssd_clean_element_complete(curr);
          break;

	  case SSD_REFRESH_ELEMENT:
		ssd_refresh_element_complete(curr);
		break;

	  case SSD_REFRESH_GANG:
		assert(0);
		break;

      default:
        fprintf(stderr, "Unrecognized event type at ssd_event_arrive\n");
        exit(1);
   }

   // fprintf (outputfile, "Exiting ssd_event_arrive\n");
}


int ssd_get_number_of_blocks (int devno)
{
   ssd_t *currdisk = getssd (devno);
   return (currdisk->numblocks);
}


int ssd_get_numcyls (int devno)
{
   ssd_t *currdisk = getssd (devno);
   return (currdisk->numblocks);
}


void ssd_get_mapping (int maptype, int devno, int blkno, int *cylptr, int *surfaceptr, int *blkptr)
{
   ssd_t *currdisk = getssd (devno);

   if ((blkno < 0) || (blkno >= currdisk->numblocks)) {
      fprintf(stderr, "Invalid blkno at ssd_get_mapping: %d\n", blkno);
      exit(1);
   }

   if (cylptr) {
      *cylptr = blkno;
   }
   if (surfaceptr) {
      *surfaceptr = 0;
   }
   if (blkptr) {
      *blkptr = 0;
   }
}


int ssd_get_avg_sectpercyl (int devno)
{
   return (1);
}


int ssd_get_distance (int devno, ioreq_event *req, int exact, int direction)
{
   /* just return an arbitrary constant, since acctime is constant */
   return 1;
}


// returning 0 to remove warning
double  ssd_get_servtime (int devno, ioreq_event *req, int checkcache, double maxtime)
{
   fprintf(stderr, "device_get_seektime not supported for ssd devno %d\n",  devno);
   assert(0);
   return 0;
}


// returning 0 to remove warning
double  ssd_get_acctime (int devno, ioreq_event *req, double maxtime)
{
   fprintf(stderr, "device_get_seektime not supported for ssd devno %d\n",  devno);
   assert(0);
   return 0;
}


int ssd_get_numdisks (void)
{
   return(numssds);
}


void ssd_cleanstats (void)
{
   int i, j;

   for (i=0; i<MAXDEVICES; i++) {
      ssd_t *currdisk = getssd (i);
      if (currdisk) {
          ioqueue_cleanstats(currdisk->queue);
          for (j=0; j<currdisk->params.nelements; j++)
              ioqueue_cleanstats(currdisk->elements[j].queue);
      }
   }
}

void ssd_setcallbacks ()
{
   ioqueue_setcallbacks();
}

int ssd_add(struct ssd *d) {
  int c;

  if(!disksim->ssdinfo) ssd_initialize_diskinfo();

  for(c = 0; c < disksim->ssdinfo->ssds_len; c++) {
    if(!disksim->ssdinfo->ssds[c]) {
      disksim->ssdinfo->ssds[c] = d;
      numssds++;
      return c;
    }
  }

  /* note that numdisks must be equal to diskinfo->disks_len */
  disksim->ssdinfo->ssds =
    (struct ssd**) realloc(disksim->ssdinfo->ssds,
        2 * c * sizeof(struct ssd *));

  bzero(disksim->ssdinfo->ssds + numssds,
    numssds);

  disksim->ssdinfo->ssds[c] = d;
  numssds++;
  disksim->ssdinfo->ssds_len *= 2;
  return c;
}


struct ssd *ssdmodel_ssd_loadparams(struct lp_block *b, int *num)
{
  /* temp vars for parameters */
  int n;
  struct ssd *result;

  if(!disksim->ssdinfo) ssd_initialize_diskinfo();

  result = (struct ssd*) malloc(sizeof(struct ssd));
  if(!result) return 0;
  bzero(result, sizeof(struct ssd));

  n = ssd_add(result);

  result->hdr = ssd_hdr_initializer;
  if(b->name)
    result->hdr.device_name = _strdup(b->name);

  lp_loadparams(result, b, &ssdmodel_ssd_mod);

  device_add((struct device_header *)result, n);
  if (num != NULL)
	  *num = n;
  return result;
}


struct ssd *ssd_copy(struct ssd *orig) {
  int i;
  struct ssd *result = (struct ssd*) malloc(sizeof(struct ssd));
  bzero(result, sizeof(struct ssd));
  memcpy(result, orig, sizeof(struct ssd));
  result->queue = ioqueue_copy(orig->queue);
  for (i=0;i<orig->params.nelements;i++)
      result->elements[i].queue = ioqueue_copy(orig->elements[i].queue);
  return result;
}

void ssd_set_syncset (int setstart, int setend)
{
}


static void ssd_acctime_printstats (int *set, int setsize, char *prefix)
{
   int i;
   statgen * statset[MAXDEVICES];

   if (device_printacctimestats) {
      for (i=0; i<setsize; i++) {
         ssd_t *currdisk = getssd (set[i]);
         statset[i] = &currdisk->stat.acctimestats;
      }
      stat_print_set(statset, setsize, prefix);
   }
}


static void ssd_other_printstats (int *set, int setsize, char *prefix)
{
   int i;
   int numbuswaits = 0;
   double waitingforbus = 0.0;

   for (i=0; i<setsize; i++) {
      ssd_t *currdisk = getssd (set[i]);
      numbuswaits += currdisk->stat.numbuswaits;
      waitingforbus += currdisk->stat.waitingforbus;
   }

   fprintf(outputfile, "%sTotal bus wait time: %f\n", prefix, waitingforbus);
   fprintf(outputfile, "%sNumber of bus waits: %d\n", prefix, numbuswaits);
}

void ssd_print_block_lifetime_distribution(int elem_num, ssd_t *s, int ssdno, double avg_lifetime, char *sourcestr)
{
    
    const int bucket_size = 20;
    int no_buckets = (100/bucket_size + 1);
    int i;
    int *hist;
    int *idle_hist;
    int dead_blocks = 0;
    int n;
    double sum,sum_idle;
    double sum_sqr,sum_sqr_idle;
    double mean;
    double variance;
    ssd_element_metadata *metadata = &(s->elements[elem_num].metadata);

    // allocate the buckets
    hist = (int *) malloc(no_buckets * sizeof(int));
    if(hist == NULL)
    {
    	fprintf(stderr, "cannot allocate memory for hist\n");
    }
    memset(hist, 0, no_buckets * sizeof(int));

    // allocate the buckets
    idle_hist = (int *) malloc(no_buckets * sizeof(int));
    if(idle_hist == NULL)
    {
    	fprintf(stderr, "cannot allocate memory for idle_hist\n");
    }
    memset(idle_hist, 0, no_buckets * sizeof(int));

    // to calc the variance
    n = s->params.blocks_per_element;
    sum = 0;
    sum_sqr = 0;
    
    sum_idle = 0;
    sum_sqr_idle = 0;
    long int written_pages_per_element = 0;
    long int total_stress_cycles = 0,total_actual_stress_cycles =0;
    long int sum_used_lifetime = 0;
		double sum_retention_period = 0,sum_sqr_retention_period = 0,l_r_p=100; //100 as in 100 years.
		double sum_recovery_period=0;
		int l_r_p_id = -1;
		double block_retention_period =0;
    for (i = 0; i < s->params.blocks_per_element; i ++) {
        int bucket,idle_bucket=0;

        int rem_lifetime = metadata->block_usage[i].rem_lifetime;
        double perc = (rem_lifetime * 100.0) / avg_lifetime;
        // find out how many blocks have completely been erased.
        if (ssd_block_dead(&metadata->block_usage[i],s)) {
            dead_blocks ++;
        }

        if (perc >= 100) {
            // this can happen if a particular block was not
            // cleaned at all and so its remaining life time
            // is greater than the average life time. put these
            // blocks in the last bucket.
            bucket = no_buckets - 1;
        } else {
            bucket = (int) perc / bucket_size;
        }

        hist[bucket] ++;
        // calculate the variance
        sum = sum + rem_lifetime;
        sum_sqr = sum_sqr + (rem_lifetime*rem_lifetime);
        sum_used_lifetime+= (SSD_MAX_ERASURES - rem_lifetime);

			block_retention_period = metadata->block_usage[i].least_retention_page->retention_period/YEAR;
			sum_retention_period += block_retention_period;
			sum_sqr_retention_period += (block_retention_period*block_retention_period);
			if(l_r_p > block_retention_period) { 
				l_r_p = block_retention_period; 
				l_r_p_id = i;
			}
      written_pages_per_element += metadata->block_usage[i].total_pages_written;
    }

    fprintf(outputfile, "%s #%d elem #%d   ", sourcestr, ssdno, elem_num);
    fprintf(outputfile, "Block Lifetime Distribution\n");
    // print the bucket size
    fprintf(outputfile, "%s #%d elem #%d   ", sourcestr, ssdno, elem_num);
    for (i = bucket_size; i <= 100; i += bucket_size) {
        fprintf(outputfile, "< %d\t", i);
    }
    fprintf(outputfile, ">= 100\t\n");

    // print the histogram bar lengths
    fprintf(outputfile, "%s #%d elem #%d   ", sourcestr, ssdno, elem_num);
    for (i = bucket_size; i <= 100; i += bucket_size) {
        fprintf(outputfile, "%d\t", hist[i/bucket_size - 1]);
    }
    fprintf(outputfile, "%d\t\n", hist[no_buckets - 1]);
    mean = sum/n;
    variance = (sum_sqr - sum*mean)/(n - 1);
    fprintf(outputfile, "%s #%d elem #%d   Total used life time:\t%ld\n",
        sourcestr, ssdno, elem_num, sum_used_lifetime);
    fprintf(outputfile, "%s #%d elem #%d   Average of life time:\t%f\n",
        sourcestr, ssdno, elem_num, mean);

	//Print retention period statistics.
	mean = sum_retention_period/n;
	variance = (sum_sqr_retention_period - sum_retention_period*mean)/(n-1);
	fprintf(outputfile, "%s #%d elem #%d   Average Retention Period:\t%lf\n",
        sourcestr, ssdno, elem_num, mean);
  fprintf(outputfile, "%s #%d elem #%d   Variance of Retention Period:\t%lf\n",
        sourcestr, ssdno, elem_num, variance);
	//Also print the blocks with the least retention period
  fprintf(outputfile, "%s #%d elem #%d   Least Retention Period:\t%lf (block #%d)\n",
        sourcestr, ssdno, elem_num, l_r_p, l_r_p_id);

  ssd_element_stat *stat = &(s->elements[elem_num].stat);
	mean = stat->tot_recovery_period/stat->tot_logical_stresses;
	fprintf(outputfile, "%s #%d elem #%d   Average Recovery Period :\t%lf (in seconds)\n",
        sourcestr, ssdno, elem_num, mean);


    fprintf(outputfile, "%s #%d elem #%d   Written Pages:\t%ld\n",
        sourcestr, ssdno, elem_num, written_pages_per_element-INITIAL_STRESSES*s->params.blocks_per_element*s->params.pages_per_block);
    fprintf(outputfile, "%s #%d elem #%d   Stress cycles:\t%ld\n",
        sourcestr, ssdno, elem_num, stat->tot_logical_stresses);
    fprintf(outputfile, "%s #%d elem #%d   Actual cycles:\t%ld\n",
        sourcestr, ssdno, elem_num, stat->tot_physical_stresses);
    fprintf(outputfile, "%s #%d elem #%d   Variance of life time:\t%f\n",
        sourcestr, ssdno, elem_num, variance);
    fprintf(outputfile, "%s #%d elem #%d   Total dead blocks:\t%d\n",
        sourcestr, ssdno, elem_num, dead_blocks);
}

//prints the cleaning algo statistics
void ssd_printcleanstats(int *set, int setsize, char *sourcestr)
{
    int i;
    int tot_ssd = 0;
    int elts_count = 0;
    double iops = 0;

    fprintf(outputfile, "\n\nSSD CLEANING STATISTICS\n");
    fprintf(outputfile, "---------------------------------------------\n\n");
    for (i = 0; i < setsize; i ++) {
        int j;
        int tot_elts = 0;
        ssd_t *s = getssd(set[i]);

        if (s->params.write_policy == DISKSIM_SSD_WRITE_POLICY_OSR) {

            elts_count += s->params.nelements;

            for (j = 0; j < s->params.nelements; j ++) {
                int plane_num;
                double avg_lifetime;
                double elem_iops = 0;
                double elem_clean_iops = 0;

                ssd_element_stat *stat = &(s->elements[j].stat);

                avg_lifetime = ssd_compute_avg_lifetime(-1, j, s);
                fprintf(outputfile, "%s #%d elem #%d   Total reqs issued:\t%d\n",
                    sourcestr, set[i], j, s->elements[j].stat.tot_reqs_issued);
                fprintf(outputfile, "%s #%d elem #%d   Total time taken:\t%f\n",
                    sourcestr, set[i], j, s->elements[j].stat.tot_time_taken);
                if (s->elements[j].stat.tot_time_taken > 0) {
                    elem_iops = ((s->elements[j].stat.tot_reqs_issued*1000.0)/s->elements[j].stat.tot_time_taken);
                    fprintf(outputfile, "%s #%d elem #%d   IOPS:\t%f\n",
                        sourcestr, set[i], j, elem_iops);
                }

                fprintf(outputfile, "%s #%d elem #%d   Total cleaning reqs issued:\t%d\n",
                    sourcestr, set[i], j, s->elements[j].stat.num_clean);
                fprintf(outputfile, "%s #%d elem #%d   Total cleaning time taken:\t%f\n",
                    sourcestr, set[i], j, s->elements[j].stat.tot_clean_time);
                fprintf(outputfile, "%s #%d elem #%d   Total migrations:\t%d\n",
                    sourcestr, set[i], j, s->elements[j].metadata.tot_migrations);
                fprintf(outputfile, "%s #%d elem #%d   Total pages migrated:\t%d\n",
                    sourcestr, set[i], j, s->elements[j].metadata.tot_pgs_migrated);
                fprintf(outputfile, "%s #%d elem #%d   Total migrations cost:\t%f\n",
                    sourcestr, set[i], j, s->elements[j].metadata.mig_cost);


                if (s->elements[j].stat.tot_clean_time > 0) {
                    elem_clean_iops = ((s->elements[j].stat.num_clean*1000.0)/s->elements[j].stat.tot_clean_time);
                    fprintf(outputfile, "%s #%d elem #%d   clean IOPS:\t%f\n",
                        sourcestr, set[i], j, elem_clean_iops);
                }

                fprintf(outputfile, "%s #%d elem #%d   Overall IOPS:\t%f\n",
                    sourcestr, set[i], j, ((s->elements[j].stat.num_clean+s->elements[j].stat.tot_reqs_issued)*1000.0)/(s->elements[j].stat.tot_clean_time+s->elements[j].stat.tot_time_taken));

                iops += elem_iops;

                fprintf(outputfile, "%s #%d elem #%d   Number of free blocks:\t%d\n",
                    sourcestr, set[i], j, s->elements[j].metadata.tot_free_blocks);
                fprintf(outputfile, "%s #%d elem #%d   Number of cleans:\t%d\n",
                    sourcestr, set[i], j, stat->num_clean);
                fprintf(outputfile, "%s #%d elem #%d   Pages moved:\t%d\n",
                    sourcestr, set[i], j, stat->pages_moved);
                fprintf(outputfile, "%s #%d elem #%d   Total xfer time:\t%f\n",
                    sourcestr, set[i], j, stat->tot_xfer_cost);
                if (stat->tot_xfer_cost > 0) {
                    fprintf(outputfile, "%s #%d elem #%d   Xfer time per page:\t%f\n",
                        sourcestr, set[i], j, stat->tot_xfer_cost/(1.0*stat->pages_moved));
                } else {
                    fprintf(outputfile, "%s #%d elem #%d   Xfer time per page:\t0\n",
                        sourcestr, set[i], j);
                }
                fprintf(outputfile, "%s #%d elem #%d   Average lifetime:\t%f\n",
                    sourcestr, set[i], j, avg_lifetime);
                fprintf(outputfile, "%s #%d elem #%d   Plane Level Statistics\n",
                    sourcestr, set[i], j);
                fprintf(outputfile, "%s #%d elem #%d   ", sourcestr, set[i], j);
                for (plane_num = 0; plane_num < s->params.planes_per_pkg; plane_num ++) {
                    fprintf(outputfile, "%d:(%d)  ",
                        plane_num, s->elements[j].metadata.plane_meta[plane_num].num_cleans);
                }
                fprintf(outputfile, "\n");

                // KJ: for now didn't printout the retention stats    
     //           ssd_print_block_lifetime_distribution(j, s, set[i], avg_lifetime, sourcestr);
                fprintf(outputfile, "\n");

                tot_elts += stat->pages_moved;
            }

            //fprintf(outputfile, "%s SSD %d average # of pages moved per element %d\n",
            //  sourcestr, set[i], tot_elts / s->params.nelements);

            tot_ssd += tot_elts;
            fprintf(outputfile, "\n");
        }
    }

    if (elts_count > 0) {
        fprintf(outputfile, "%s   Total SSD IOPS:\t%f\n",
            sourcestr, iops);
        fprintf(outputfile, "%s   Average SSD element IOPS:\t%f\n",
            sourcestr, iops/elts_count);
    }

    //fprintf(outputfile, "%s SSD average # of pages moved per ssd %d\n\n",
    //  sourcestr, tot_ssd / setsize);
}


void ssd_printrefreshstats(int *set, int setsize, char *sourcestr)
{
  int i;
  int tot_ssd = 0;
  int elts_count = 0;
  double iops = 0;

  fprintf(outputfile, "\n\nSSD REFRESH STATISTICS\n");
  fprintf(outputfile, "---------------------------------------------\n\n");
  for (i = 0; i < setsize; i ++) {
    int j;
    int tot_elts = 0;
    ssd_t *s = getssd(set[i]);

    if (s->params.refresh_policy == SSD_FCFS_REFRESH || s->params.refresh_policy == SSD_EDF_REFRESH) {
      for (j = 0; j < s->params.nelements; j ++) {
        ssd_element_stat *stat = &(s->elements[j].stat);
        fprintf(outputfile,"%s SSD #%d elem #%d Number of refresh operations:%d\n",
                sourcestr, set[i], j, stat->num_refresh);
        fprintf(outputfile,"%s SSD #%d elem #%d Total refresh time:%lf\n",
                sourcestr, set[i], j, stat->tot_refresh_time);              
      }
      fprintf(outputfile, "\n");
    }
  }
}

void ssd_printsetstats (int *set, int setsize, char *sourcestr)
{
   int i;
   struct ioq * queueset[MAXDEVICES*SSD_MAX_ELEMENTS];
   int queuecnt = 0;
   int reqcnt = 0;
   char prefix[80];

   //using more secure functions
   sprintf_s4(prefix, 80, "%sssd ", sourcestr);
   for (i=0; i<setsize; i++) {
      ssd_t *currdisk = getssd (set[i]);
      struct ioq *q = currdisk->queue;
      queueset[queuecnt] = q;
      queuecnt++;
      reqcnt += ioqueue_get_number_of_requests(q);
   }
   if (reqcnt == 0) {
      fprintf (outputfile, "\nNo ssd requests for members of this set\n\n");
      return;
   }
   ioqueue_printstats(queueset, queuecnt, prefix);

   ssd_acctime_printstats(set, setsize, prefix);
   ssd_other_printstats(set, setsize, prefix);
}


void ssd_printstats (void)
{
   struct ioq * queueset[MAXDEVICES*SSD_MAX_ELEMENTS];
   int set[MAXDEVICES];
   int i,j;
   int reqcnt = 0;
   char prefix[80];
   int diskcnt;
   int queuecnt;

   fprintf(outputfile, "\nSSD STATISTICS\n");
   fprintf(outputfile, "---------------------\n\n");

   sprintf_s3(prefix, 80, "ssd ");

   diskcnt = 0;
   queuecnt = 0;
   for (i=0; i<MAXDEVICES; i++) {
      ssd_t *currdisk = getssd (i);
      if (currdisk) {
         struct ioq *q = currdisk->queue;
         queueset[queuecnt] = q;
         queuecnt++;
         reqcnt += ioqueue_get_number_of_requests(q);
         diskcnt++;
      }
   }
   assert (diskcnt == numssds);

   if (reqcnt == 0) {
      fprintf(outputfile, "No ssd requests encountered\n");
      return;
   }

   //ioqueue_printstats(queueset, queuecnt, prefix);

   diskcnt = 0;
   for (i=0; i<MAXDEVICES; i++) {
      ssd_t *currdisk = getssd (i);
      if (currdisk) {
         set[diskcnt] = i;
         diskcnt++;
      }
   }
   assert (diskcnt == numssds);

   ssd_acctime_printstats(set, numssds, prefix);
   ssd_other_printstats(set, numssds, prefix);

   ssd_printcleanstats(set, numssds, prefix);
   // KJ: eleminate refresh stats as well
 //  ssd_printrefreshstats(set,numssds,prefix);

   fprintf (outputfile, "\n\n");

   for (i=0; i<numssds; i++) {
      ssd_t *currdisk = getssd (set[i]);
      if (currdisk->printstats == FALSE) {
          continue;
      }
      reqcnt = 0;
      {
          struct ioq *q = currdisk->queue;
          reqcnt += ioqueue_get_number_of_requests(q);
      }
      if (reqcnt == 0) {
          fprintf(outputfile, "No requests for ssd #%d\n\n\n", set[i]);
          continue;
      }
      fprintf(outputfile, "ssd #%d:\n\n", set[i]);
      sprintf_s4(prefix, 80, "ssd #%d ", set[i]);
      {
          struct ioq *q;
          q = currdisk->queue;
          ioqueue_printstats(&q, 1, prefix);
      }
      for (j=0;j<currdisk->params.nelements;j++) {
          char pprefix[100];
          struct ioq *q;
          sprintf_s5(pprefix, 100, "%s elem #%d ", prefix, j);
          q = currdisk->elements[j].queue;
          ioqueue_printstats(&q, 1, pprefix);
      }
      ssd_acctime_printstats(&set[i], 1, prefix);
      ssd_other_printstats(&set[i], 1, prefix);
      fprintf (outputfile, "\n\n");

			//ll_release(currdisk->refresh_queue);
			//currdisk->refresh_queue = 0;
   }
}

void ssd_dump_block_status();

void ssd_notify_trace_done()
{
  ssd_dump_block_status();
}

void ssd_dump_block_status()
{
	int ssd=0,el=0,bl=0;
	assert(snapshotfile);
	const char *deadoralive[] = {"A", "D"};
  fprintf(snapshotfile,"#Retention period,Pages written,#Stresses,#Actual Stresses,delvth\n");
	fprintf(snapshotfile,"#---------------------------\n");
	for(ssd=0;ssd<MAXDEVICES;ssd++)
	{
		ssd_t* currdisk = getssd(ssd);
		if(currdisk) {
			for(el=0;el<currdisk->params.nelements;el++) {
				ssd_element_metadata *metadata;
				metadata = &(currdisk->elements[el].metadata);
				for(bl=0;bl<currdisk->params.blocks_per_element;bl++) {
          block_metadata *currblock = &metadata->block_usage[bl];
          ssd_page_metadata *lrp = currblock->least_retention_page;
          //NOTE: If you update any of the following fprintf's, make sure
          //to change traceanalyzer.py's parserow() and updateblock()
					fprintf(snapshotfile,"D #%d E #%d P #%d B #%d\n",ssd,currblock->elem_num,
															currblock->plane_num,
															currblock->block_num);
					fprintf(snapshotfile,"%s %ld %ld",deadoralive[currblock->state==SSD_BLOCK_DEAD],
  																					currblock->logical_stresses,
                                            currblock->physical_stresses);
          fprintf(snapshotfile," %lf %d %d %d %lf %.6lf\n",
                                          lrp->delvth,
                                          lrp->logical_stresses,
                                          lrp->physical_stresses,
                                          lrp->stress_increment,
                                          lrp->eqn_cycle,
																					lrp->retention_period/YEAR);
				}
			}
		}
	}
  fprintf(snapshotfile,"#---------------------------\n");
  fflush(snapshotfile);
  fclose(snapshotfile);
  snapshotfile = 0;
  if(disksim->curr_snapshot_output_file)
    fprintf(stderr,"Saved snapshot in %s\n",disksim->curr_snapshot_output_file);
}


// returning 0 to remove warning
double ssd_get_seektime (int devno,
                ioreq_event *req,
                int checkcache,
                double maxtime)
{
  fprintf(stderr, "device_get_seektime not supported for ssd devno %d\n",  devno);
  assert(0);
  return 0;
}

// KJ: deal with reading from SSD snapshots and get the respected info in an array
void get_retention_time_from_snapshot(double *retention, int array_length, FILE *snapshot)
{
	char line[201];
	char line1[201];
	char line2[201];
	int i = 0;
//	char *currline = NULL;
	
	// for now, only the retention time needed
	char dead_or_alive;
	long int blk_phy, blk_logical;
	int page_phy, page_logical, stress_inc;
	double delvth, equ_cycle;
	
	// skip the first 2 lines;
	if(snapshot != NULL)
	{
		assert(fgets(line, 200, snapshot));
		assert(fgets(line, 200, snapshot));
	}
	
	// get the retention into array
	while((fgets(line1, 200, snapshot) != NULL) && (fgets(line2, 200, snapshot) != NULL))
	{
		// line2 is needed; last is retention;
		assert(i < array_length);
		assert(sscanf(line2, "%c %ld %ld %lf %d %d %d %lf %lf\n", &dead_or_alive, &blk_phy, &blk_logical, &delvth, &page_phy, &page_logical, &stress_inc, &equ_cycle, &(retention[i++])) == 9);
		// retention time is in the magnitude of YEAR		
	}
	// TODO: not sure needed
	rewind(snapshot);
	
}

// KJ: get min and max of an array goes here
double get_min(double *array, int length)
{
	int i;
	assert(array != NULL);
	double min = array[0];
	for(i=1; i<length; i++)
	{
		if(array[i] < min)
			min = array[i];
	}
        fprintf(stderr, "The min is %lf\n", min);
	return min;
}

double get_max(double *array, int length)
{
	int i;
	assert(array != NULL);
	double max = array[0];
	for(i=1; i<length; i++)
	{
		if(array[i] > max)
			max = array[i];
	}
        fprintf(stderr, "The max is %lf\n", max);
	return max;
}

// KJ: this function filters the blocks with retention time under max retention limit;
int get_retention_under_limit(double *retention, int retention_array_length, double *retention_under_limit, double max_retention)
{
	int i;
        int length_of_filtered_array = 0;
        for(i=0; i<retention_array_length; i++)
	{
		if(retention[i] <= max_retention)
		{
			retention_under_limit[length_of_filtered_array] = retention[i];
                        length_of_filtered_array++;
		}
	}
	return length_of_filtered_array;
}

// KJ: This is the function to calculate the normalized histogram of given array; 
void get_normed_histogram(double *retention_under_limit, int retention_under_limit_length, int nbins, double *hist, double *bin_edge)
{
	double min = get_min(retention_under_limit, retention_under_limit_length);
	double max = get_max(retention_under_limit, retention_under_limit_length);
	double bin_size = (max - min) / (double)nbins;	

	// get bin edge
        int i;
	for(i=0; i<(nbins+1); i++)
	{
		bin_edge[i] = min + i * bin_size;
	}
	
	// set hist array to 0
	for(i=0; i<nbins; i++)
	{
		hist[i] = 0;
	}      	

	// get hist array
	int index;
	for(i=0; i<retention_under_limit_length; i++)
	{
		index = (int) floor((retention_under_limit[i] - min) / bin_size);
		hist[index]++;
	}
	// normalize it
	for(i=0; i<nbins; i++)
	{
		hist[i] = hist[i] / (double)retention_under_limit_length;
	}   
	
  	
}

// KJ: compute the earth mover distance of 2 given SSD snapshots here
double ssd_calculate_distance(FILE *prev_snapshot, FILE *curr_snapshot)
{

	// devno is 0
	int ssd = 0;
	ssd_t *currdisk = getssd(ssd);
	
	fprintf(stderr, "get rentention time array\n");
	// get the retention time array from 2 snapshots
	int retention_array_length = currdisk->params.nelements * currdisk->params.blocks_per_element;
	double *prev_retention = (double *)malloc(retention_array_length * sizeof(double));
	double *curr_retention = (double *)malloc(retention_array_length * sizeof(double));

        double *prev_retention_under_limit = (double *)malloc(retention_array_length * sizeof(double));
	double *curr_retention_under_limit = (double *)malloc(retention_array_length * sizeof(double));

	get_retention_time_from_snapshot(prev_retention, retention_array_length, prev_snapshot);
	get_retention_time_from_snapshot(curr_retention, retention_array_length, curr_snapshot);

        double max_retention = 35;
        int prev_retention_under_limit_length = get_retention_under_limit(prev_retention, retention_array_length, prev_retention_under_limit, max_retention);
        int curr_retention_under_limit_length = get_retention_under_limit(curr_retention, retention_array_length, curr_retention_under_limit, max_retention);
        // cannot calculate histogram using opencv; do it by hand;
        
        int nbins = 3000;
        double *prev_hist = (double *)malloc(nbins * sizeof(double));
        double *curr_hist = (double *)malloc(nbins * sizeof(double));
        double *prev_bin_edge = (double *)malloc((nbins + 1) * sizeof(double));
        double *curr_bin_edge = (double *)malloc((nbins + 1) * sizeof(double));
        
        get_normed_histogram(prev_retention_under_limit, prev_retention_under_limit_length, nbins, prev_hist, prev_bin_edge);
        get_normed_histogram(curr_retention_under_limit, curr_retention_under_limit_length, nbins, curr_hist, curr_bin_edge);

      
	/*
        // compute histogram
	// creat cvMat for retention array;
	CvMat prev_ma = cvMat(1, retention_array_length, CV_32FC1, prev_retention);
	CvMat curr_ma = cvMat(1, retention_array_length, CV_32FC1, curr_retention);
        CvMat *prev_ma_ptr[] = {&prev_ma};
        CvMat *curr_ma_ptr[] = {&curr_ma};
        // Mat format won't work for C version; turned to image format
    
        fprintf(stderr, "Generate retention image\n");
        IplImage *prev_img = cvCreateImage(cvSize(retention_array_length, 1), IPL_DEPTH_32F, 1); 	
        IplImage *curr_img = cvCreateImage(cvSize(retention_array_length, 1), IPL_DEPTH_32F, 1);
        int j;
        CvScalar *s1 = (CvScalar *)malloc(retention_array_length * sizeof(CvScalar));
        CvScalar *s2 = (CvScalar *)malloc(retention_array_length * sizeof(CvScalar));
        for(j=0; j<retention_array_length; j++)
        {
            (s1[j]).val[0] = prev_retention[j];
            cvSet2D(prev_img, 0, j, s1[j]);
            (s2[j]).val[0] = curr_retention[j];
            cvSet2D(curr_img, 0, j, s2[j]);
        }
        
        IplImage *prev_img_ptr[] = {prev_img};
        IplImage *curr_img_ptr[] = {curr_img};
      */  
       /*
	// create histogram
	int dims = 1;
	// This is the bin number;
	int sizes[] = {3000};
	// TODO: get_min and get_max
        fprintf(stderr, "Calc min and max\n");
	float prev_ranges[] = {get_min(prev_retention, retention_array_length), get_max(prev_retention, retention_array_length)};
	float curr_ranges[] = {get_min(curr_retention, retention_array_length), get_max(curr_retention, retention_array_length)};
        float *prev_ranges_ptr[] = {prev_ranges};
        float *curr_ranges_ptr[] = {curr_ranges};
	fprintf(stderr, "Create histogram\n");
        CvHistogram* prev_hist = NULL;
	CvHistogram* curr_hist = NULL;
	prev_hist = cvCreateHist(dims, sizes, CV_HIST_ARRAY, prev_ranges_ptr, 1);
	curr_hist = cvCreateHist(dims, sizes, CV_HIST_ARRAY, curr_ranges_ptr, 1);
	
        fprintf(stderr, "Calc histogram\n");
	// Calculate histogram
     
	cvCalcHist(prev_img_ptr, prev_hist, 0, 0);
	cvCalcHist(curr_img_ptr, curr_hist, 0, 0);
     */
       /*	
        cvCalcHist(prev_ma_ptr, prev_hist, 0, 0);
	cvCalcHist(curr_ma_ptr, curr_hist, 0, 0);
       */
 
	CvMat *prev_sig = cvCreateMat(nbins, 2, CV_32FC1);
	CvMat *curr_sig = cvCreateMat(nbins, 2, CV_32FC1);
	
        fprintf(stderr, "Calc signiture\n");
	// set both signitures; which has 3000 rows and 2 columns; each row of signiture stores weight followed by coordinates;
	// not sure if should pass in the first parameter as a pointer or not;
	int i;
        for(i=0; i<nbins; i++)
        {
        	cvmSet(prev_sig, i, 0, prev_hist[i]);
            //    fprintf(stderr, "Value associated with bin for prev snapshot %d is %lf\n", i, prev_hist[i]);
                cvmSet(prev_sig, i, 1, (prev_bin_edge[i] + prev_bin_edge[i+1]) / 2.0 );
	//	fprintf(stderr, "Bin range from %lf to %lf for prev snapshot\n", prev_bin_edge[i], prev_bin_edge[i+1]);
      		cvmSet(curr_sig, i, 0, curr_hist[i]);
	//	fprintf(stderr, "Value associated with bin for curr snapshot %d is %lf\n", i, curr_hist[i]);
                cvmSet(curr_sig, i, 1, (curr_bin_edge[i] + curr_bin_edge[i+1]) / 2.0 );
          //	fprintf(stderr, "Bin range from %lf to %lf for curr snapshot\n", curr_bin_edge[i], curr_bin_edge[i+1]);
        }



/*
	for(i=0; i<sizes[0]; i++)
	{
		cvmSet(prev_sig, i, 0, cvQueryHistValue_1D(prev_hist, i));
                fprintf(stderr, "Value associated with bin for prev snapshot %d is %lf\n", i, cvQueryHistValue_1D(prev_hist, i));
		cvmSet(prev_sig, i, 1, ((cvGetHistValue_1D(prev_hist, i))[0] + (cvGetHistValue_1D(prev_hist, i))[1]) / 2.0);
                fprintf(stderr, "Bin range from %lf to %lf for prev snapshot\n", (cvGetHistValue_1D(prev_hist, i))[0], (cvGetHistValue_1D(prev_hist, i))[1]);
		cvmSet(curr_sig, i, 0, cvQueryHistValue_1D(curr_hist, i));
                fprintf(stderr, "Value associated with bin for curr snapshot %d is %lf\n", i, cvQueryHistValue_1D(curr_hist, i));
		cvmSet(curr_sig, i, 1, ((cvGetHistValue_1D(curr_hist, i))[0] + (cvGetHistValue_1D(curr_hist, i))[1]) / 2.0);
                fprintf(stderr, "Bin range from %lf to %lf for curr snapshot\n", (cvGetHistValue_1D(curr_hist, i))[0], (cvGetHistValue_1D(curr_hist, i))[1]);
	}
*/

        fprintf(stderr, "Calc earth mover distance\n");
	// Calculate the earth mover distance
	float earth_mover_distance = cvCalcEMD2(prev_sig, curr_sig, CV_DIST_L1, NULL, NULL, NULL, NULL, NULL);
	fprintf(stderr, "Earth mover distance is %lf\n", earth_mover_distance);

        // free signiture
        CvMat *sig[] = {prev_sig, curr_sig};
        cvReleaseMat(sig);       

	// free memory
        free(prev_retention);
        free(curr_retention);
        free(prev_retention_under_limit);
        free(curr_retention_under_limit);
        free(prev_hist);
        free(curr_hist);
        free(prev_bin_edge);
        free(curr_bin_edge);
	
	return earth_mover_distance;

}

//void ssd_duplicate_prev_round(int trace_repeat);


// KJ: a simple fast_forward strategy;
void ssd_duplicate_prev_round(int trace_repeat)
{
	// try to get the duration if one round of tracefile
	//char line[501];
	int pageno = 0;
	int blkno = 0;
	int elemno = 0;
	double access_time = 0;
	double first_time = 0;
	double last_time = 0;
	double duration = 0;
	int cur_repeat = 0;
//	double cur_time = 0;

        static int is_first_time_to_enter_fastforward = 1;
	int sub_status = 1;
        int ssd = 0;
        ssd_t *currdisk = getssd(ssd);	

	// get the first request time;
	assert(disksim->curr_record);
	first_time = disksim->curr_record->head[0].access_time;
//	assert((disksim->prev_trace_file = fopen(disksim->prev_trace_file_name, "r")) != NULL);
//	assert(fgets(line, 500, disksim->prev_trace_file) != NULL);
	
//	assert(sscanf(line, "P #%d B #%d E #%d T #%lf\n", &pageno, &blkno, &elemno, &first_time) == 4);
	fprintf(stderr, "first_time %lf\n", first_time);
	
	// get the last request time;
//	while(fgets(line, 500, disksim->prev_trace_file) != NULL)
//	{
//		assert(sscanf(line, "P #%d B #%d E #%d T #%lf\n", &pageno, &blkno, &elemno, &last_time) == 4);
	//	fprintf(stderr, "Pageno %d, blkno %d, elemno %d, first_time %lf\n", pageno, blkno, elemno, last_time);
//	}
	last_time = disksim->curr_record->head[disksim->curr_record->index - 1].access_time;
	fprintf(stderr, "last_time %lf\n", last_time);
	// calc the duration
	duration = last_time - first_time;
        if(is_first_time_to_enter_fastforward)
        {
		disksim->duration = duration;
                is_first_time_to_enter_fastforward = 0;
	}
	fprintf(stderr, "Duration of fast forward mode is %lf\n", duration);
//	rewind(disksim->prev_trace_file);
	
	// get the current SSD with the device number to be 0;
	
	ssd_element_metadata *metadata = NULL;
        block_metadata *currblock = NULL;
	ssd_page_metadata *currpage = NULL;

	
	int i = 0;	
	while(cur_repeat < trace_repeat)
	{
		// if not the end of a round
	//	if(fgets(line, 500, disksim->prev_trace_file) != NULL)
		if(i < disksim->curr_record->index)
		{
		//	assert(sscanf(line, "P #%d B #%d E #%d T #%lf\n", &pageno, &blkno, &elemno, &access_time) == 4);
		//	assert(disksim->curr_record->head[i]);
			pageno = disksim->curr_record->head[i].pageno;
			blkno = disksim->curr_record->head[i].blkno;
			elemno = disksim->curr_record->head[i].elemno;
			access_time = disksim->curr_record->head[i].access_time;
			i++;
			simtime = access_time + (cur_repeat + 1) * duration;
			// get element
			metadata = &(currdisk->elements[elemno].metadata);
			// get block
			currblock = &(metadata->block_usage[blkno]);
			// get page
			currpage = &(currblock->page[pageno]);
			// update page stress info;
			// This function updates all the info related to stress; not lpn; the mapping remains the same
			ssd_update_stress_info(currpage, currblock, simtime, currdisk);
			
		}else
		{
			// if it is the end of round
			cur_repeat++;
			fprintf(stderr, "Fast forward Simulation Pass %d has completed!\n", cur_repeat);
		//	rewind(disksim->prev_trace_file);
			i = 0;
		
		}
		
		// checkpoint here:
		// if meet the time to dump output file and snapshots file;
		if(time_to_checkpoint(simtime))
		{
			fprintf(stderr, "Fast Forward Simulation Writing snapshots at %lf\n", simtime);
			// deals with both snapshots and output file;
			set_snapshot_file();
			// dump snapshots;
			ssd_notify_trace_done();
			// change next checkpoint starttime to current time;
			checkpoint_starttime = simtime;
		//	disksim->checkpoint_starttime = cur_time;
			// dump stats in output file;
			disksim_printstats2(sub_status);
                        
			
		}
		
	}
	// TODO: remove the prev_trace_file of the completed round; remember to remove after debugging done;
	//remove(disksim->prev_trace_file_name);
}


// KJ: fast forward; contains all the strategies;
int ssd_fast_forward_simulation(int cur_repeat_trace_count)
{
	// get the actual trace_repeat counts in this fast_forward round
	int trace_repeat_rounds = ((cur_repeat_trace_count+disksim->fast_forward_duration) > disksim->repeat_trace)? 
	(disksim->repeat_trace-cur_repeat_trace_count) : disksim->fast_forward_duration;
	fprintf(stderr, "This fast forward round has %d trace repeat times\n", trace_repeat_rounds);
	
	// browse all the fast forward strategy
	switch(disksim->fast_forward_strategy)
	{
		case 0:
		// the simplest strategy: simply duplicate the nearest detailed simulation round trace_repeat_rounds times;
		ssd_duplicate_prev_round(trace_repeat_rounds);
		break;
	}
	
	return trace_repeat_rounds;
}

/* default ssd dev header */
struct device_header ssd_hdr_initializer = {
  DEVICETYPE_SSD,
  sizeof(struct ssd),
  "unnamed ssd",
  (void*)ssd_copy,
  ssd_set_depth,
  ssd_get_depth,
  ssd_get_inbus,
  ssd_get_busno,
  ssd_get_slotno,
  ssd_get_number_of_blocks,
  ssd_get_maxoutstanding,
  ssd_get_numcyls,
  ssd_get_blktranstime,
  ssd_get_avg_sectpercyl,
  ssd_get_mapping,
  ssd_event_arrive,
  ssd_get_distance,
  ssd_get_servtime,
  ssd_get_seektime,
  ssd_get_acctime,
  ssd_bus_delay_complete,
  ssd_bus_ownership_grant,
  ssd_notify_trace_done,
  // KJ: added for ff simulation
  ssd_fast_forward_simulation,
  ssd_calculate_distance
};
