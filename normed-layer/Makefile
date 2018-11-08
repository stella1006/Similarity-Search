-include ./faiss/makefile.inc

CC=g++

CXXFLAGS += -std=c++11
CXXFLAGS += -fopenmp -g -O3 -Wall -I. -msse2
LDLIBS += -lboost_program_options  -lpthread -lm

.PHONY: all clean

PROGS= exknn

all:    $(PROGS)

%: %.cpp ./faiss/libfaiss.so
	     $(CC) -o $@ $(CXXFLAGS) $(LDLIBS) $^ $(LDFLAGS) $(LIBS)
clean:
	rm -f $(PROGS)
