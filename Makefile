EXE_NAME=bam_subsample

CPP=g++

CPPFLAGS=-O3 -std=c++11 -I$(HOME)/boost/include -I$(HOME)/sources/htslib 
LDFLAGS=-L$(HOME)/boost/lib -L$(HOME)/lib -lhts -lboost_program_options

all: pre $(EXE_NAME) done

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:.cpp=.o)

$(EXE_NAME): $(OBJS)
	$(CPP) -o $(EXE_NAME) $^ $(LDFLAGS)

.PHONY: clean pre done

clean:
	rm bam_subsample src/*.o

pre:
	cd src/

done:
	cd ..
