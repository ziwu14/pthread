CC=g++
CFLAGS=-std=gnu++11 -g -pthread

SRCS=$(wildcard *.cpp)
OBJS=$(patsubst %.cpp,%.o,$(SRCS))

PRGS=rainfall_seq rainfall_pt

#------------ EXECUTABLE FILES ---------------
.PHONY: all clean 

all : $(PRGS)

rainfall_seq:rainfall_seq.o
	$(CC) $(CFLAGS) -o $@ $<
rainfall_pt:rainfall_pt.o
	$(CC) $(CFLAGS) -o $@ $<

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

clean:
	rm -f $(PRGS) *.o *~
