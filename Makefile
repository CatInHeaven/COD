cod:main.o RefSys.o GPSTime.o matrix.o util.o

	gcc -o cod main.o RefSys.o GPSTime.o matrix.o util.o -lm

main.o:main.c main.h

	gcc -c main.c

RefSys.o:RefSys.c RefSys.h
	gcc -c RefSys.c

GPSTime.o:GPSTime.c GPSTime.h
	gcc -c GPSTime.c

matrix.o:matrix.c matrix.h
	gcc -c matrix.c

util.o:util.c util.h
	gcc -c util.c

clean:

	rm cod main.o RefSys.o GPSTime.o matrix.o util.o