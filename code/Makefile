CXX=clang++
CFLAGS=-O3 -Wall -Wextra -Werror -std=c++11 -stdlib=libc++ -c
LDFLAGS=-std=c++11 -stdlib=libc++ -lsndfile
SRCS=src/filter.cpp

all: filter

filter: filter.o
	      $(CXX) $(LDFLAGS) filter.o -o filter

filter.o: $(SRCS)
	      $(CXX) $(CFLAGS) $(SRCS)

clean:
	      rm *.o filter
