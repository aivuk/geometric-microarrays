CFLAGS=-O3 -Wall -DHAVE_INLINE
LDFLAGS=-lgsl -lm -lgslcblas -static
COMPILER=gcc

ma_geom.exe: ma_geom.o
	${COMPILER} -o ma_geom.exe ma_geom.o $(LDFLAGS) $(CFLAGS)

clean:
	rm ma_geom.exe

