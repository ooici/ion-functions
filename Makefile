CC=gcc
SRCDIR=extensions
CFLAGS=-std=c99 -g -ggdb -Wall -I$(SRCDIR)

test_objects=$(SRCDIR)/test.o $(SRCDIR)/spike.o

all: $(SRCDIR)/test

$(SRCDIR)/test: $(test_objects)
	$(CC) -o $(SRCDIR)/test $(test_objects)

$(test_objects): %.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $<

clean:
	rm -rf $(SRCIDR)/*.o
	rm -rf $(SRCDIR)/test

