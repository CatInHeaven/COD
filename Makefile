CC = mpicc

objects = main.o RefSys.o GPSTime.o matrix.o util.o
TARGET := cod

$(TARGET):$(objects)

	$(CC) -o $(TARGET) -g $(objects) -lm

main.o:main.h

RefSys.o:RefSys.h

GPSTime.o:GPSTime.h

matrix.o:matrix.h

util.o:util.h

.PHONY : clean
clean:

	-rm $(TARGET) $(objects)