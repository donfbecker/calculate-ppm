CC=gcc

PROGRAMS = count-pulses calculate-ppm

all: $(PROGRAMS)

count-pulses: count-pulses.c
	$(CC) count-pulses.c -o count-pulses -lm

calculate-ppm: calculate-ppm.c
	$(CC) -I /usr/include/ffmpeg -o calculate-ppm calculate-ppm.c -L /usr/lib/ffmpeg -lavcodec -lavformat -lavutil -lm -lfftw3

clean:
	rm -f $(PROGRAMS)

