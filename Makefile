MAKE = make
CC = gcc 
CFLAGS = -O3 -Wall -march=native -std=c99 

SOURCE = reader-el.c

EXE = st-exact

all: $(EXE)

st-exact: $(SOURCE)
	$(CC) $(CFLAGS) -o $@ $< -lm
	
clean:  
	rm -f *.o *.a *~ 
	rm -f $(EXE)
