CC      ?= gcc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -I../include

.PHONY: all test clean

all: test_spvsh320

test_spvsh320: test_spvsh320.c ../src/spvsh320.c ../include/spvsh320.h
	$(CC) $(CFLAGS) -o $@ test_spvsh320.c ../src/spvsh320.c

test: test_spvsh320
	./test_spvsh320

clean:
	rm -f test_spvsh320
