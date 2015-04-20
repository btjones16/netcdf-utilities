CC = clang
CFLAGS = -Wall -std=c11 -g

CFLAGS += `nc-config --cflags`
LDFLAGS += `nc-config --libs`

all:
	$(CC) $(CFLAGS) -c -o log.c.o log.c
	$(CC) $(CFLAGS) -c -o subset.c.o subset.c
	$(CC) $(LDFLAGS) *.c.o -o subset
