#include <stdio.h>
#include <math.h>
#include "rcfilter.h"

void RCFilter_init(RCFilter *filt, char type, float cutoffFreqHz, float sampleTimeS) {
	filt->type = type;

	/* Compute equivalent 'RC' constant from cut-off frequency */
	float RC = 1.0f / (M_PI * 2 * cutoffFreqHz);

	switch(type) {
		case RC_FILTER_LOWPASS:
			/* Pre-compute filter coefficients for first-order low-pass filter */
			filt->coeff[0] = sampleTimeS / (RC + sampleTimeS);
			filt->coeff[1] = RC / (RC + sampleTimeS);
			break;

		case RC_FILTER_HIGHPASS:
			/* Pre-compute filter coefficients for first-order low-pass filter */
			filt->coeff[0] = RC / (RC + sampleTimeS);
			/* Initial value must be 1 */
			filt->coeff[1] = 1.0f;
			break;
	}

	/* Clear output buffer */
	filt->out[0] = 0.0f;
	filt->out[1] = 0.0f;
}

float RCFilter_update(RCFilter *filt, float inp) {
	/* Shift output samples */
	filt->out[1] = filt->out[0];

	switch(filt->type) {
		case RC_FILTER_LOWPASS:
			filt->out[0] = (filt->coeff[0] * inp) + (filt->coeff[1] * filt->out[1]);
			break;

		case RC_FILTER_HIGHPASS:
			// Use coeff[1] to save last input
			filt->out[0] = (filt->coeff[0] * filt->out[1]) + (filt->coeff[0] * (inp - filt->coeff[1]));
			filt->coeff[1] = inp;
			break;
	}

	/* Return filtered sample */
	return filt->out[0];
}
