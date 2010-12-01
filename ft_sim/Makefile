# MAKEFILE FOR THE MROI FRINGE TRACKING SIMULATOR

##### May need to edit this section #####

# Library paths
INCFLAGS_CFITSIO = -I/usr/include -I/usr/local/include
LDLIBS_CFITSIO = -L/usr/lib -L/usr/local/lib -lcfitsio
INCFLAGS = -I../common `pkg-config --cflags fftw3` $(INCFLAGS_CFITSIO)
LDLIBS = `pkg-config --libs fftw3` $(LDLIBS_CFITSIO)

# Compilers to use
CC = gcc
CPP = g++

# Compiler flags
CFLAGS = -g -O3 -Wall $(INCFLAGS)
CPPFLAGS = -ggdb -g -O3 -Wall --pedantic -Wextra $(INCFLAGS)

##### end of user-editable section #####


# Rules to make .o file from .c file
%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : ../common/%.c
	$(CC) $(CFLAGS) -c $<

# Rules to make .o file from .cpp file
%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $<
%.o : ../common/%.cpp
	$(CPP) $(CPPFLAGS) -c $<


# DEFAULT TARGET
default: simulator

# SIMULATOR
simulator: simulator.o stations.o source.o matrix.o random.o \
  atmosphere.o zernike.o delaylines.o modulator.o \
  beam.o beamcombiner.o detector.o \
  FTBasden.o FTPedretti.o FTphaseABCD1.o FTphaseABCD2.o
	$(CPP) $(CPPFLAGS) $^ -o $@ $(LDLIBS)

# CLEANER
clean:
	rm -f *.o
