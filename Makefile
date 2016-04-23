CC=cilk++
CFLAGS=-O3

all: Aloop Bloop Cloop Dloop

Aloop: 
	$(CC) $(CFLAGS) -o Aloop_CPU AloopFW.cilk -lcilkutil

Bloop: 
	$(CC) $(CFLAGS) -o Bloop_CPU BloopFW.cilk -lcilkutil

Cloop: 
	$(CC) $(CFLAGS) -o Cloop_CPU CloopFW.cilk -lcilkutil

Dloop: 
	$(CC) $(CFLAGS) -o Dloop_CPU DloopFW.cilk -lcilkutil

clean:
	rm -rf *_CPU  core* Test* a_* *.csv *.plt
