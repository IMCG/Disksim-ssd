#include "disksim_iosim_param.h"
#include <libparam/bitvector.h>
#include "../disksim_iosim.h"
static int DISKSIM_IOSIM_IO_TRACE_TIME_SCALE_depend(char *bv) {
return -1;
}

static void DISKSIM_IOSIM_IO_TRACE_TIME_SCALE_loader(int result, double d) { 
 ioscale = d;

}

static int DISKSIM_IOSIM_IO_MAPPINGS_depend(char *bv) {
return -1;
}

static void DISKSIM_IOSIM_IO_MAPPINGS_loader(int result, struct lp_list *l) { 
if (! (iosim_load_mappings(l))) { // foo 
 } 

}

void * DISKSIM_IOSIM_loaders[] = {
(void *)DISKSIM_IOSIM_IO_TRACE_TIME_SCALE_loader,
(void *)DISKSIM_IOSIM_IO_MAPPINGS_loader
};

lp_paramdep_t DISKSIM_IOSIM_deps[] = {
DISKSIM_IOSIM_IO_TRACE_TIME_SCALE_depend,
DISKSIM_IOSIM_IO_MAPPINGS_depend
};

