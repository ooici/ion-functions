CC=gcc
SRCDIR=extensions
override CFLAGS+=-std=c99 -g -ggdb -Wall -I$(SRCDIR) 
override LDFLAGS+=-lm

test_objects=$(SRCDIR)/test.o $(SRCDIR)/spike.o $(SRCDIR)/utils.o $(SRCDIR)/gradient.o $(SRCDIR)/time_utils.o $(SRCDIR)/GeomagnetismLibrary.o $(SRCDIR)/wmm.o $(SRCDIR)/polycals.o

all: $(SRCDIR)/test

$(SRCDIR)/test: $(test_objects)
	$(CC) -o $(SRCDIR)/test $(LDFLAGS) $(test_objects)

$(test_objects): %.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $<

clean:
	rm -rf $(SRCDIR)/*.o
	rm -rf $(SRCDIR)/test

