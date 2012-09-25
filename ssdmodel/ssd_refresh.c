#include "ssd_clean.h"
#include "ssd.h"
#include "ssd_utils.h"
#include "ssd_refresh.h"


int refresh_queue_insert(ssd_t *currdisk,void *data_to_insert) {
  listnode *start = currdisk->refresh_queue;
  switch(currdisk->params.refresh_policy) {
    case SSD_FCFS_REFRESH:
      refresh_queue_fcfs_insert(start,data_to_insert);
      break;
    case SSD_EDF_REFRESH:
      refresh_queue_edf_insert(start,data_to_insert);
      break;
    case SSD_REFRESH_NONE:
      break;
    default:
      assert(0);
  }
  return 0;
}

void ssd_perform_refresh(ssd_t *currdisk, double now)
{
  int size = ll_get_size(currdisk->refresh_queue);
  int i=0,blocks_to_refresh=0;


  //double next_refresh_time = currdisk->params.refresh_interval + currdisk->params.refresh_service_time * currdisk->params.nelements * currdisk->params.blocks_per_element/1000; //This value is an upper bound on the refresh time.
  listnode **clean_blocks_issue_list = (listnode**)malloc(currdisk->params.nelements * sizeof(listnode*));

  for(i=0;i<currdisk->params.nelements;i++)
      ll_create(&clean_blocks_issue_list[i]);

  i=0;
  while(i<size) {

    listnode* currnode = ll_get_nth_node(currdisk->refresh_queue,i);
    block_metadata* currblock = (block_metadata*)currnode->data;
    ssd_page_metadata* least_retention_page = currblock->least_retention_page;
    //assert(now - least_retention_page->time_of_last_stress >0);
    //least_retention_page->retention_period -= (now - least_retention_page->time_of_last_stress); //account for time spent idling.

    if(ssd_block_dead(currblock,currdisk)) //Too late, block is dead. Handle it accordingly
      continue;

    //check if currblock needs refresh *now*
    if(ssd_virtual_retention_expired(currblock,currdisk)) {
	  //pick this block to refresh only if the corresponding chip is not busy.
	  //else just skip it. the next time when the refresh operation is triggered,
	  //the blocks will be refreshed.
	  if (currdisk->elements[currblock->elem_num].media_busy != TRUE ) {
		refresh_queue_free_node(currdisk->refresh_queue,currblock);
		blocks_to_refresh++;size--;
		ll_insert_at_tail(clean_blocks_issue_list[currblock->elem_num],currblock);
	  }
    }
    i++;
  }

  if(blocks_to_refresh!=0) {
    fprintf(stderr, "# of blocks to refresh :%d\n",blocks_to_refresh);
  }

  for(i=0;i<currdisk->params.nelements;i++) {
    if(is_queue_empty(clean_blocks_issue_list[i]))
      continue;
    //fprintf(stderr, "About to refresh %d blocks in elem #%d\n",ll_get_size(clean_blocks_issue_list[i]),i);
    ssd_invoke_element_refresh_fcfs(i,clean_blocks_issue_list[i],currdisk);
    ll_release(clean_blocks_issue_list[i]);
  }
}
void ssd_invoke_element_refresh_fcfs(int elem_num,listnode *blocks_to_refresh,ssd_t *currdisk)
{
  int block_len = 0;
  int i;
  double cost = 0;

  block_len = ll_get_size(blocks_to_refresh);
  block_metadata *currblock = 0;
  listnode *currnode = 0;
  //for all blocks in every element.
  for(i=0;i<block_len;i++) {
    currnode = ll_get_nth_node(blocks_to_refresh,i);
    currblock = (block_metadata*)currnode->data;
    ASSERT(currblock->elem_num  == elem_num);
    cost+= ssd_refresh_block(currblock,currdisk); //sum every cost, because we are not applying refresh in batches
  }

  ssd_element *elem = &currdisk->elements[elem_num];
  if (cost > 0) {
      ioreq_event *tmp;
      elem->media_busy = TRUE;

      // we use the 'blkno' field to store the element number
      tmp = (ioreq_event *)getfromextraq();
      tmp->devno = currdisk->devno;
      tmp->time = simtime + cost;
      tmp->blkno = elem_num;
      tmp->ssd_elem_num = elem_num;
      tmp->type = SSD_REFRESH_ELEMENT;
      tmp->flags = SSD_REFRESH_ELEMENT;
      tmp->busno = -1;
      tmp->bcount = -1;
      stat_update (&currdisk->stat.acctimestats, cost);
      addtointq ((event *)tmp);
      // stat
      elem->stat.tot_refresh_time += cost;
  }
}


double ssd_refresh_block(block_metadata *block_metadata, ssd_t *s)
{
    double cost = 0;
    int refresh_invoked = 0;
    int is_refresh = 1;
    int elem_num = block_metadata->elem_num;
    int plane_num = block_metadata->plane_num;
    ssd_element *elem = &s->elements[elem_num];
    ssd_element_metadata *metadata = &(s->elements[elem_num].metadata);
    plane_metadata *pm = &metadata->plane_meta[plane_num];

    // element must be free
    ASSERT(elem->media_busy == FALSE);

	  // calculate cost of refresh operation.
    // Check if the plane is already being cleaned or refreshed. If so, skip this block.
    ASSERT(pm->clean_in_progress == FALSE);
    metadata->active_page = metadata->plane_meta[plane_num].active_page;
    cost = _ssd_clean_block_fully(block_metadata->block_num, plane_num,
                                          elem_num, metadata, s,is_refresh);
    ASSERT(pm->clean_in_progress == FALSE);
    return cost;
}
