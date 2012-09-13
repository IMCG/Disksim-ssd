/*
 * ssd_refresh.h
 *
 *  Created on: Sep 16, 2011
 *      Author: vb
 */

#ifndef SSD_REFRESH_H_
#define SSD_REFRESH_H_
#include "ssd.h"
#include "ssd_utils.h"

void ssd_perform_refresh(ssd_t *currdisk, double now);
void ssd_invoke_element_refresh_fcfs(int elem_num,listnode *blocks_to_refresh,ssd_t *currdisk);
double ssd_refresh_block(block_metadata *block_metadata, ssd_t *s);
int refresh_queue_insert(ssd_t *currdisk,void *data_to_insert);

#endif /* SSD_REFRESH_H_ */
