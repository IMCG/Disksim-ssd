

include .paths

CC = gcc
CFLAGS = -g -O $(DISKSIM_CFLAGS) $(DISKMODEL_CFLAGS) $(LIBPARAM_CFLAGS) $(LIBDDBG_CFLAGS) -Wall -Wno-unused -Wimplicit -Werror

LDFLAGS = $(DISKSIM_LDFLAGS)  $(MEMSMODEL_LDFLAGS) $(SSDMODEL_LDFLAGS) -L.. -ldiskmodel  $(LIBPARAM_LDFLAGS) $(LIBDDBG_LDFLAGS) -lm -lz

TARGETS = g4_skews
all: $(TARGETS)

%: %.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

%: %.c


clean:
	rm -f *.o $(TARGETS)

distclean: clean
	rm -f *~
