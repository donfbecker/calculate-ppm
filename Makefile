CC=gcc

PROGRAMS = count-pulses calculate-ppm

all: $(PROGRAMS)

count-pulses: count-pulses.c
	$(CC) count-pulses.c -o count-pulses -lm -std=gnu99

calculate-ppm: calculate-ppm.c rcfilter.c rcfilter.h
	$(CC) -I /usr/include/ffmpeg -o calculate-ppm calculate-ppm.c rcfilter.c -L /usr/lib/ffmpeg -lavcodec -lavformat -lavutil -lm -lfftw3 -std=gnu99

clean:
	rm -f $(PROGRAMS)

