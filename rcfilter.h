#ifndef RC_FILTER_H
#define RC_FILTER_H

#define RC_FILTER_LOWPASS 0
#define RC_FILTER_HIGHPASS 1

typedef struct {
	char type;
	float coeff[2];
	float out[2];
} RCFilter;

void RCFilter_init(RCFilter *filt, char type, float cutoffFreqHz, float sampleTimeS);
float RCFilter_update(RCFilter *filt, float inp);

#endif
