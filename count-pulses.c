#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

float calculate_power(float *samples, int length) {
    float power = 0;
    for(int i = 0; i < length; i++) {
        power += pow(samples[i], 2);
    }
    return sqrt(power / length);
}

float correlate_blocks(float *a, float *b, int length) {
	float correlation = 0;

	for(int i = 0; i < length; i++) {
		for(int j = 0; j < length - i; j++) {
			//correlation += a[i] * b[j + i];
			correlation += fabs(a[i] - b[j + i]);
		}
	}

	// After correlating, copy to the old block
	memcpy(a, b, sizeof(float) * length);
	return fabs(correlation / (length * length));
}

int main(int argc, char **argv) {
    int sample_rate        = 44100;
    int length             = 64;
    int min_pulse_duration = 10;
    float threshold        = 0.075;

    int opt;
    while((opt = getopt(argc, argv, "c:f:l:p:r:t:")) != -1) {
        switch(opt) {
            case 'l': length = atoi(optarg); break;
            case 'p': min_pulse_duration = atoi(optarg); break;
            case 'r': sample_rate = atoi(optarg); break;
            case 't': threshold = atof(optarg); break;
        }
    }

    freopen(NULL, "rb", stdin);
    freopen(NULL, "wb", stdout);

    // Calibrate out phase offset from alternating sampling
    int min_pulse_samples = (sample_rate / 1000) * min_pulse_duration;

    size_t bytes;
    float buffer[length], last_block[length];

    int pulse_samples = 0;
    int sample_number = 0;
    int last_pulse = 0;

    int pulse_count = 0;
    float pulse_delay = 0.0;

    while((bytes = fread(&buffer, sizeof(float), length, stdin)) > 0) {
        //float power = calculate_power(buffer, length);
	float power = correlate_blocks(last_block, buffer, length);
        //fprintf(stderr, "power = %f\n", power);

	//for(int i = 0; i < length; i++) fwrite(&power, sizeof(float), 1, stdout);

        if(power >= threshold) {
            pulse_samples += length;
        } else {
            if(pulse_samples >= min_pulse_samples) {
                if(last_pulse > 0) {
                    pulse_count++;
                    pulse_delay += ((float)(sample_number - last_pulse) / sample_rate);
                }
                last_pulse = sample_number;
            }
            pulse_samples = 0;
        }

        sample_number += length;
    }

    float avg_delay = (pulse_delay / pulse_count);
    float ppm = 60.0f / avg_delay;
    printf("{\"delay\": %f, \"ppm\": %f, \"temp\": %0.2f}\n", avg_delay, ppm, (ppm * (9.0/5.0)) + 32);
    return 0;
}
