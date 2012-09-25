// DiskSim SSD support
// ©2008 Microsoft Corporation. All Rights Reserved

#include "ssd_utils.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////////////
//                 code for bit manipulation routines
//////////////////////////////////////////////////////////////////////////////


//clears a particular bit.
//pos 0 corresponds to msb
void ssd_clear_bit(unsigned char *c, int pos)
{
    int byte = pos / 8;
    int bit_pos = 7 - (pos % 8);

    c[byte] = (c[byte] & (~(0x1 << bit_pos)));

    return;
}


//sets a particular bit.
//pos 0 corresponds to msb and pos 7 corresponds
//to lsb.
void ssd_set_bit(unsigned char *c, int pos)
{
    int byte = pos / 8;
    int bit_pos = 7 - (pos % 8);

    c[byte] = (c[byte] | (0x1 << bit_pos));

    return;
}

//returns true if a bit is set
int ssd_bit_on(unsigned char *c, int pos)
{
    int byte = pos / 8;
    int bit_pos = 7 - (pos % 8);

    return (c[byte] & (0x1 << bit_pos));
}

//finds the position of the first zero-th
//bit in an unsigned char array. the positions are
//numbered starting from 0 from msb to lsb. returns -1
//if all the bits are already set. 'total' specifies the
//number of bits to consider in the array. 'start' gives
//the location from which to start the search.
int ssd_find_zero_bit(unsigned char *c, int total, int start)
{
    int bit = start;
    int temp = total;

    while (temp) {
        int byte = bit / 8;
        ASSERT(byte<1024);
        if ((c[byte] >> (7 - (bit % 8))) & 0x1) {
		//if (ssd_bit_on(c,bit)) {
            bit = (bit + 1) % total;
            temp --;
        } else {
            return bit;
        }
    }

    return -1;
}

//////////////////////////////////////////////////////////////////////////////
//             adding some code for a linked list module
//////////////////////////////////////////////////////////////////////////////

void ll_create(listnode **start)
{
    *start = (listnode *)malloc(sizeof(listnode));
    memset(*start, 0, sizeof(listnode));

    (*start)->data = (header_data *)malloc(sizeof(header_data));
    ((header_data *)(*start)->data)->size = 0;
    ((header_data *)(*start)->data)->cache_index = -1;
    ((header_data *)(*start)->data)->cache_node = 0;
}

//1 - if empty.
//0- otherwise.

int is_queue_empty(listnode *start) {
  assert(start);
  return !((header_data *)(start->data))->size;
}


// Free all the nodes in the list.
void ll_release(listnode *start)
{
    free(start->data);

    while (start) {
        listnode *next = start->next;

        if (next) {
            next->prev = NULL;

            if (next == start) {
                next = NULL;
            }
        }

        if (start->prev) {
            start->prev->next = NULL;
        }

        /* we just release the node. we don't release the
         * data that is contained in the node. that is the
         * responsibility of the function that is using
         * this linked list */
        free(start);

        start = next;
    }

    return;
}

listnode *_ll_insert_at_tail(listnode *start, listnode *toinsert)
{
    if ((!start) || (!toinsert)) {
        fprintf(stderr, "Error: invalid value to _ll_insert_at_tail\n");
        exit(-1);
    } else {
        if (start->prev) {
            listnode *prevnode = start->prev;

            toinsert->next = prevnode->next;
            prevnode->next->prev = toinsert;
            prevnode->next = toinsert;
            toinsert->prev = prevnode;
        } else {
            start->prev = toinsert;
            start->next = toinsert;
            toinsert->prev = start;
            toinsert->next = start;
        }

        return toinsert;
    }
}

listnode *_ll_insert_at_head(listnode *start, listnode *toinsert)
{
    if ((!start) || (!toinsert)) {
        fprintf(stderr, "Error: invalid value to _ll_insert_at_head\n");
        exit(-1);
    } else {
        if (start->next) {
            listnode *nextnode = start->next;

            toinsert->prev = nextnode->prev;
            nextnode->prev->next = toinsert;
            nextnode->prev = toinsert;
            toinsert->next = nextnode;
        } else {
            start->next = toinsert;
            start->prev = toinsert;
            toinsert->next = start;
            toinsert->prev = start;
        }

        return toinsert;
    }
}

/*
 * Insert the data at the tail.
 */
listnode *ll_insert_at_tail(listnode *start, void *data)
{
    /* allocate a new node */
    listnode *newnode = (listnode*)malloc(sizeof(listnode));
    newnode->data = data;

    /* increment the number of entries */
    ((header_data *)(start->data))->size ++;

    return _ll_insert_at_tail(start, newnode);
}

/*
 * Insert the data at the head.
 */
listnode *ll_insert_at_head(listnode *start, void *data)
{
    /* allocate a new node */
    listnode *newnode = (listnode*)malloc(sizeof(listnode));
    newnode->data = data;

    /* increment the number of entries */
    ((header_data *)(start->data))->size ++;

    return _ll_insert_at_head(start, newnode);
}


listnode *ll_insert_at_index(listnode *start,int index,void *data) {//indices are 0-based.
//index specifies where the new data would be after insertion.
  
    if (index == 0)
      return ll_insert_at_head(start,data);
    else if (index == ((header_data *)(start->data))->size-1 )
      return ll_insert_at_tail(start, data);
    else { //do insert at head after getting to the right index.
      int i = 0;
      /* allocate a new node */
      listnode *newnode = (listnode*)malloc(sizeof(listnode));
      newnode->data = data;
      ((header_data *)(start->data))->size ++;

      listnode *head = start;
      while(i<index-1) {
        head = head->next;
        i++;
      }
      return _ll_insert_at_head(head,newnode);
    }
}
/*
 * Release a node from the linked list.
 */
void ll_release_node(listnode *start, listnode *node)
{
    if ((!node) || (!start)) {
        fprintf(stderr, "Error: invalid items on the list\n");
        exit(-1);
    } else {
        assert( (unsigned)((header_data *)(start->data))->size < 131072);
        ((header_data *)(start->data))->size --;

        if (((header_data *)(start->data))->size > 0) {
            listnode *nodesprev = node->prev;

            nodesprev->next = node->next;
            node->next->prev = nodesprev;
        } else {
            /* do some sanity checking */

            /* there is just one element left in the list
             * and we want to release it. so, check if it
             * matches the prev and next of the start node. */
            if ((start->next != node) || (start->prev != node)) {
                fprintf(stderr, "Error: sanity check failed\n");
                exit(-1);
            }

            /* this is the last element - just release it */
            start->next = NULL;
            start->prev = NULL;
        }

        if(node == ((header_data *)(start->data))->cache_node) {
            ((header_data *)(start->data))->cache_node = NULL;
            ((header_data *)(start->data))->cache_index = -1;
        }
        /* we just release the node. we don't release the
         * data that is contained in the node. that is the
         * responsibility of the function that is using
         * this linked list */
        free(node);
    }
}

/*
 * Release the tail node from the linked list.
 */
void ll_release_tail(listnode *start)
{
    listnode *tail = start->prev;

    if ((!tail) || (tail == start)) {
        fprintf(stderr, "Warning: the list is empty. no element to release\n");
        return;
    } else {
        ll_release_node(start, tail);
    }
}

listnode *ll_get_tail(listnode *start)
{
    return (listnode*)(start->prev->data);
}


/*
 * Return the total number of nodes in the list.
 */
int ll_get_size(listnode *start)
{
    return ((header_data *)(start->data))->size;
}

/*
 * For debugging purposes.
 * We can make some improvements on performance,
 * for e.g., by searching from the tail if n is
 * greater than size/2.
 */
listnode *ll_get_nth_node(listnode *start, int n)
{
    int i=0;
    int reverse=0;
    listnode *node;
    int use_cache = 0;
    int cache_index = ((header_data *)(start->data))->cache_index;
    listnode* cache_node = ((header_data *)(start->data))->cache_node;
    int size = ll_get_size(start);   

    if (n >= size) {
        fprintf(stderr, "Error: n (%d) is greater than size %d\n",
            n, ll_get_size(start));
        return NULL;
    }
  
    if(cache_index!=-1 && cache_node) 
    {
      if( abs(n - cache_index) < abs(n - size/2)) {
        use_cache = 1;
        node = cache_node;
        i = cache_index;
        if(n<cache_index)
          reverse = 1;
      }
    }
    use_cache = 0;
    if (!use_cache)
    {
      if (n > size/2) {
          i = size - 1;
          reverse = 1;
          node = start->prev;
      } else {
          i = 0;
          reverse = 0;
          node = start->next;
      }
    }

    /* go over each node and print it */
    while ((node) && (node != start)) {
        if (i == n) {
            //cache_index = i;
            //cache_node = node;
            return node;
        }

        if (reverse) {
            node = node->prev;
            i --;
        } else {
            node = node->next;
            i ++;
        }
    }

    fprintf(stderr, "Error: cannot find the %d node in list\n", n);
    return NULL;
}

/*
 * REFRESH QUEUE POLICIES.
*/

//returns the length of the queue.
int refresh_queue_fcfs_insert(listnode *start, void *data_to_insert) {
	assert(start && data_to_insert);
  
  block_metadata *data = (block_metadata*)data_to_insert;
  assert(data->entry_in_refresh_queue==NULL);
    

	listnode * added_node = ll_insert_at_tail(start,data_to_insert);
	//store this ptr in the data to random access the queue.
	data->entry_in_refresh_queue = (void*)(added_node);
  return ((header_data *)(start->data))->size;
}

int refresh_queue_edf_insert(listnode *start, void *data_to_insert) {
  listnode *added_node = 0;
	assert(start && data_to_insert);

  //if already in the queue, remove this copy and insert again.
  block_metadata *data = (block_metadata*)data_to_insert;
  if (data->entry_in_refresh_queue!=NULL)
    assert(!refresh_queue_free_node(start, data));
    
  //identify the index where the data should be inserted into the queue.
  int index = 0;
  double deadline = ((block_metadata*)data_to_insert)->least_retention_page->retention_period;
  listnode *current = start->next;
  if(!current)
    added_node = ll_insert_at_index(start,0,data_to_insert);
  else {
    double current_deadline = ((block_metadata*)current->data)->least_retention_page->retention_period;
    while(deadline>current_deadline) { //current->next!=start -> is the condition for reaching the end of the list.
      if ( current->next == start ) {
        index = ((header_data *)(start->data))->size;
        break;
      } else {
        current = current->next;
        index++;
        current_deadline = ((block_metadata*)current->data)->least_retention_page->retention_period;
      }
    }
    added_node = ll_insert_at_index(start,index,data_to_insert);
  }
	//store this ptr in the data to random access the queue.
	data->entry_in_refresh_queue = (void*)(added_node);
  return ((header_data *)(start->data))->size;
}

void* refresh_queue_pop(listnode *start) {
  void *data = start->next->data;
	refresh_queue_free_node(start,data);
  return data;
}

//Refresh queue's are not queue's in the strict sense.
//Whenever a block is erased, it should be removed from the queue.
//So,it is possible to remove nodes from the middle of the queue.
int refresh_queue_free_node(listnode *start, void *data) {
	//search the queue until you find data and delete it when you find one.
	block_metadata *bm = (block_metadata*)data;
	assert(bm->entry_in_refresh_queue);
	assert( ((listnode*)(bm->entry_in_refresh_queue))->data == data);
  listnode* node_to_release = (listnode*)(bm->entry_in_refresh_queue);
#ifdef DEBUG
  static int debug = 0;
  if (debug == 1) {
    listnode* current = start->next;
    //check the list to ensure the node does exist in the list.
    int match =0;
    while(current!=start) {
      if (current == node_to_release) {
        match = 1;
        break;
      }
      current= current->next;
    }
    assert(match!=0);
  }
#endif
	ll_release_node(start,node_to_release);
	bm->entry_in_refresh_queue = 0;
	return 0;
}






