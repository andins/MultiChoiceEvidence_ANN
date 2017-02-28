TARGETS = darco mfield_darco
# deltanu_deltalambda
OBJECTS = $(TARGETS:.cc=.o)
INCLUDE = /usr/include/gsl
HEADERS = brain.h pool.h setup.h 
CC      = g++
#CFLAGS  = -ggdb -O3 -march=pentium4 -ffast-math -W -Wall -pedantic -ansi -Winline -I$(INCLUDE)
#CFLAGS  = -ggdb -ffast-math -W -Wall -pedantic -ansi -Winline -I$(INCLUDE) #FLAGS TO BE USED
CFLAGS = -O2 -pipe -march=native -fomit-frame-pointer -ffast-math -W -Wall -pedantic -ansi -Winline
LIBDIRS = -L/home/andrea/usr/lib
LFLAGS  =  $(LIBDIRS) -lm -lgsl -lgslcblas
#LFLAGS  =  $(LIBDIRS) -lm -lgsl -lgslcblas -static
DESTDIR=$(HOME)/usr/

INSTALL_BIN = $(DESTDIR)/bin # Place to put the binaries


# $(@) stands for the target of the current rule, $(^) stands for all the
# dependencies of the current rule, and $(<) stands for the first dependency of
# the current rule (the first thing in $(^)).  #spiking meanfield
all: $(TARGETS)
#	install $(TARGETS)  $(INSTALL_BIN)
	ctags *.{cc,h}
	cscope -b

darco: brain.o pool.o setup.o darco.o
	$(CC) $^ -o $@ $(LFLAGS)
mfield_darco: brain.o pool.o setup.o mfield_darco.o
	$(CC) $^ -o $@ $(LFLAGS)

#spiking.o: brain.h pool.h spiking.cc
#	$(CC) $(CFLAGS) -c spiking.cc
#deltanu_deltalambda: brain.o pool.o setup.o deltanu_deltalambda.o;
#	$(CC) $(LFLAGS) brain.o pool.o setup.o deltanu_deltalambda.o -o deltanu_deltalambda $(LFLAGS)

darco.o: darco.cc brain.h pool.h setup.h
# Here I can use setup/setup according to my needs.
setup.o: setup.cc brain.h setup.h 
brain.o: brain.cc brain.h pool.h
pool.o: pool.cc pool.h
mfield_darco.o: mfield_darco.cc brain.h pool.h setup.h
%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -f $(TARGETS)
	rm -f *.o
	rm -f *~
	rm -f *.{log,mpx,dvi,aux}
	rm -f dmarti.tar.bz2
pack:
	tar jcvf dmarti.tar.bz2 *.cc *.h Makefile
