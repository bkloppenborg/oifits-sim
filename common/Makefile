

# Create a list of all of the .cpp files:
CPP_SOURCES = $(wildcard *.cpp)
CPP_OBJECTS = $(patsubst %.cpp, %.o, $(CPP_SOURCES))

# Create a list of all of the .c files:
CC_SOURCES = $(wildcard *.c)
CC_OBJECTS = $(patsubst %.c, %.o, $(CC_SOURCES))

# Compilers:
CPP = g++
CPPFLAGS = -ggdb -g -Wall --pedantic -Wextra

CC = gcc
CCFLAGS = -g -Wall

############# Make Targets below #############

default: cc_common cpp_common libOISim.a

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $<
	
%.o : %.c
	$(CC) $(CCFLAGS) -c $<
	
cc_common: $(CC_OBJECTS)

cpp_common: $(CPP_OBJECTS)
	
libOISim.a:	cc_common cpp_common
	rm -f $@
	ar cq $@ $(CPP_OBJECTS) $(CC_OBJECTS)
	
# CLEANER
.PHONY: clean
clean:
	rm -f *.o
	rm -f *.a
