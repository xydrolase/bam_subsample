EXE_NAME=bam_subsample

CPP=g++

CPPFLAGS=-O3 -std=c++11 -I$(HOME)/boost/include -I$(HOME)/sources/htslib
LDFLAGS=-L$(HOME)/boost/lib -L$(HOME)/lib -lhts -lboost_program_options

all: $(EXE_NAME)

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:.cpp=.o)

$(EXE_NAME): $(OBJS)
	$(CPP) -o $(EXE_NAME) $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm bam_subsample src/*.o
