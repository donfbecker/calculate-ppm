#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <math.h>

#include <fftw3.h>

// Make sure this is after fftw3.h
#include <complex.h>

#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libavutil/avutil.h>

#include "rcfilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX_PULSES 512
#define TOLERANCE 5
bool enable_debug = false;

typedef struct _Pulse {
    uint32_t start;
    uint32_t end;
    float max_power;
    float min_power;
} Pulse;

typedef struct _PulseSummary {
    int count;
    float frequency;
    float delay;
    float ppm;
    float temperature;
    float confidence;
} PulseSummary;

void debug_log(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    if(enable_debug) vfprintf(stderr, fmt, args);
    va_end(args);
}

float calculate_power(float *samples, int length) {
    float power = 0;
    for(int i = 0; i < length; i++) {
        power += pow(samples[i], 2);
    }
    return sqrt(power / length);
}

// This version has been optimized
_Complex float goertzel_dft_kernel(const float *samples, size_t length, int sample_rate, float frequency) {
    float omega = -2.0f * M_PI * (float)frequency / (float)sample_rate;
    float coeff = 2.0f * cosf(omega);
    _Complex float exp = cexpf(I * omega);
    _Complex float s1 = 0.0f + 0.0f * I;
    _Complex float s2 = 0.0f + 0.0f * I;

    for (size_t n = 0; n < length; n++) {
        _Complex float s = samples[n] + coeff * s1 - s2;
        s2 = s1;
        s1 = s;
    }

    return (s1 - s2 * exp) / length;
}


FILE *decode_audio_file(char *path, int *sample_rate) {
    FILE *tmp_fh;
    AVFormatContext *fmt_ctx = NULL;
    AVCodecContext *codec_ctx = NULL;
    const AVCodec *codec = NULL;
    AVPacket *packet = NULL;
    AVFrame *frame = NULL;
    int audio_stream_idx = -1;

    tmp_fh = tmpfile();
    if(tmp_fh == NULL) {
        fprintf(stderr, "Could not create temporary file.");
        return NULL;
    }

    if (avformat_open_input(&fmt_ctx, path, NULL, NULL) < 0) {
        fprintf(stderr, "Could not open %s\n", path);
        return NULL;
    }

    if (avformat_find_stream_info(fmt_ctx, NULL) < 0) {
        fprintf(stderr, "Error looking up stream info.\n");
        return NULL;
    }

    // 2. Find the audio stream
    for (unsigned int i = 0; i < fmt_ctx->nb_streams; i++) {
        if (fmt_ctx->streams[i]->codecpar->codec_type == AVMEDIA_TYPE_AUDIO) {
            audio_stream_idx = i;
            break;
        }
    }

    codec = avcodec_find_decoder(fmt_ctx->streams[audio_stream_idx]->codecpar->codec_id);
    codec_ctx = avcodec_alloc_context3(codec);
    avcodec_parameters_to_context(codec_ctx, fmt_ctx->streams[audio_stream_idx]->codecpar);
    *sample_rate = codec_ctx->sample_rate;

    if (avcodec_open2(codec_ctx, codec, NULL) < 0) {
        fprintf(stderr, "avcodec_open2 error.\n");
        return NULL;
    }

    packet = av_packet_alloc();
    frame  = av_frame_alloc();
    while (av_read_frame(fmt_ctx, packet) >= 0) {
        if (packet->stream_index == audio_stream_idx) {
            if (avcodec_send_packet(codec_ctx, packet) >= 0) {
                while (avcodec_receive_frame(codec_ctx, frame) >= 0) {
                    float *samples = (float *)frame->data[0];
                    fwrite(samples, sizeof(float), frame->nb_samples, tmp_fh);
                }
            }
        }
        av_packet_unref(packet);
    }

    av_packet_free(&packet);
    av_frame_free(&frame);
    avcodec_free_context(&codec_ctx);
    avformat_close_input(&fmt_ctx);

    fseek(tmp_fh, 0, SEEK_SET);
    return tmp_fh;
}

float calculate_ppm(Pulse *pulse, int pulse_count, int sample_rate) {
    int pulse_delay = 0, avg_delay, ppm;

    for(int p = 1; p < pulse_count; p++) pulse_delay += pulse[p].start - pulse[p - 1].start;
    avg_delay = ((float)pulse_delay / (pulse_count - 1)) / sample_rate;
    ppm = 60.0f / avg_delay;

    return ppm;
}

float calculate_power_of_blocks(FILE *fh, int window_size, int block_size, int sample_rate, float target_frequency, float *power) {
    float samples[window_size];
    float max_power = 0.0f;
    float p;
    size_t bytes;
    uint32_t b = 0;

    int block_offset = window_size - block_size;

    fseek(fh, 0, SEEK_SET);
    if(block_offset > 0) fread(samples, sizeof(float), block_offset, fh);
    while((bytes = fread(&samples[block_offset], sizeof(float), block_size, fh)) > 0) {
        if(target_frequency == 0) {
            p = calculate_power(samples, window_size);
        } else {
            p = cabs(goertzel_dft_kernel(samples, window_size, sample_rate, target_frequency));
        }
        if(p > max_power) max_power = p;
        power[b++] = p;

        if(block_offset > 0) memmove(samples, &samples[block_size], sizeof(float) * block_offset);
    }

    return max_power;
}

uint32_t detect_pulses(float *power, uint32_t block_count, float threshold, uint32_t min_pulse_blocks, Pulse *pulse) {
    uint32_t pulse_count  = 0;
    uint32_t pulse_blocks = 0;

    for(int i = 0; i < block_count; i++) {
        if(power[i] >= threshold) {
            if(pulse_blocks == 0) {
                pulse[pulse_count].start = i;
                pulse[pulse_count].max_power = power[i];
                pulse[pulse_count].min_power = power[i];
            }
            pulse_blocks++;

            if(power[i] > pulse[pulse_count].max_power) pulse[pulse_count].max_power = power[i];
            if(power[i] < pulse[pulse_count].min_power) pulse[pulse_count].min_power = power[i];
        } else {
            if(pulse_blocks >= min_pulse_blocks) {
                if(pulse_count < MAX_PULSES) pulse[pulse_count++].end = i;
            }
            pulse_blocks = 0;
        }
    }

    return pulse_count;
}

void pulse_report(FILE *fh, Pulse *pulse, int pulse_count, int total_pulse_count, int sample_rate, int block_size, PulseSummary *summary) {
    size_t file_size = 0;
    int pulse_delay = 0;
    float freq_sum = 0.0f;
    int fft_size = 512;
    float buffer[fft_size];
    double in[fft_size];
    float amplitudes[fft_size];
    fftw_complex out[fft_size];
    fftw_plan fftwp = fftw_plan_dft_r2c_1d(fft_size, in, out, FFTW_ESTIMATE);

    fseek(fh, 0, SEEK_END);
    file_size = ftell(fh);

    for(int p = 0; p < pulse_count; p++) {
        int p2 = p - 1;
        if(p > 0) pulse_delay += pulse[p].start - pulse[p2].start;

        fseek(fh, ((pulse[p].start + 1) * block_size) * sizeof(float), SEEK_SET);
        fread(buffer, sizeof(float), fft_size, fh);
        for(int i = 0; i < fft_size; i++) in[i] = buffer[i] * 0.5 * (1.0 - cos(2 * M_PI * i / fft_size));
        fftw_execute(fftwp);

        float hpb = (float)sample_rate / fft_size;
        int low_bin = floor(500.0f / hpb);
        int high_bin = ceil(4000.0f / hpb);
        float max_power = 0.0f;
        int max_index = 0;
        for(int i = 0; i < fft_size / 2; i++) {
            float power = sqrtf((out[i][0] * out[i][0]) + (out[i][1] * out[i][1]));
            if(i >= low_bin && i <= high_bin && power > max_power) {
                max_power = power;
                max_index = i;
            }
            amplitudes[i] = power;
        }

        float y1 = amplitudes[max_index - 1];
        float y2 = amplitudes[max_index];
        float y3 = amplitudes[max_index + 1];
        float y  = 0.5 * (y1 - y3) / (y1 - (2.0 * y2) + y3);
        float freq = ((max_index + y) * (sample_rate / fft_size));
        freq_sum += freq;

        debug_log("%.1f Hz pulse with power %f/%f starting at %i (%0.3f seconds) ending at %i (%i samples, %0.3fs)", freq, pulse[p].min_power, pulse[p].max_power, pulse[p].start * block_size, (((float)pulse[p].start * block_size) / sample_rate), pulse[p].end * block_size, (pulse[p].end - pulse[p].start) * block_size, ((float)(pulse[p].end - pulse[p].start) * block_size) / sample_rate);
        if(p > 0) debug_log(", %i samples (%0.3f seconds) since last pulse", (pulse[p].end - pulse[p2].end) * block_size, ((float)(pulse[p].end - pulse[p2].end) * block_size) / sample_rate);
        debug_log("\n");
    }

    fftw_destroy_plan(fftwp);

    float avg_delay = (((float)pulse_delay / (pulse_count - 1)) * block_size) / sample_rate;
    float ppm = 60.0f / avg_delay;
    float avg_freq = freq_sum / pulse_count;
    int possible_pulses = floor(((float)file_size - (pulse[0].start * block_size)) / sizeof(float) / sample_rate / avg_delay) + 1;

    debug_log("file_size = %i\n", file_size);

    summary->count       = pulse_count;
    summary->frequency   = avg_freq;
    summary->delay       = avg_delay;
    summary->ppm         = ppm;
    summary->temperature = (ppm * (9.0/5.0)) + 32;
    summary->confidence  = ((float)pulse_count / possible_pulses) * 100;

}

int find_pulse_chain(Pulse *pulse, int pulse_count, int min_chain_length, int max_chain_length, int min_distance, int max_distance, int tolerance) {
    int chain[max_chain_length], longest_chain[max_chain_length];
    int chain_length = 2, longest_chain_length = 0;
    int l1 = 0, l2 = 1, l3 = 2;
    int dt = ((pulse[l2].end - pulse[l1].end) + (pulse[l2].start - pulse[l1].start)) / 2;
    int d = 0;

    chain[0] = l1;
    chain[1] = l2;

    // Can't find a chain if pulse_count is too low
    if(pulse_count < 3 || pulse_count < min_chain_length) return 0;

    while(l1 < pulse_count - (min_chain_length - 1)) {
        bool end_of_chain = false;
        // Using the average distance between starts and ends ensures uniform pulse length also
        d = ((pulse[l3].end - pulse[l2].end) + (pulse[l3].start - pulse[l2].start)) / 2;

        if(d >= min_distance && d <= max_distance && abs(dt - d) <= tolerance) {
            // Valid link found
            chain[chain_length] = l3;
            chain_length++;
            l2 = l3;
            l3++;
        } else if(d > dt) {
            // Further chain is imposible from these nodes
            end_of_chain = true;
        } else {
            // Check next possible link
            l3++;
        }

        if(end_of_chain || l3 == pulse_count) {
            // If we have enough links, see if it's the longest chain
            if(chain_length >= min_chain_length && chain_length > longest_chain_length) {
                debug_log("Found new pulse chain with length %i\n", chain_length);
                memcpy(longest_chain, chain, sizeof(int) * chain_length);
                longest_chain_length = chain_length;
            }

            // Look for a longer chain
            l1++;
            l2 = l1 + 1;
            l3 = l2 + 1;
            dt = ((pulse[l2].end - pulse[l1].end) + (pulse[l2].start - pulse[l1].start)) / 2;
            chain_length = 2;
            chain[0] = l1;
            chain[1] = l2;
        }
    }

    // Did we reach the end of the loop?
    if(chain_length > longest_chain_length) {
        memcpy(longest_chain, chain, sizeof(int) * chain_length);
        longest_chain_length = chain_length;
    }

    // Iterate over longest chain and rearrange pulses
    for(int i = 0; i < longest_chain_length; i++) {
        pulse[i] = pulse[longest_chain[i]];
    }

    // No chain found
    return longest_chain_length;
}

int main(int argc, char **argv) {
    char *filename;
    int block_size         = 64;
    int window_size        = 256;
    int min_pulse_duration = 15;
    float threshold        = 0.05;
    int required_pulses    = 5;
    int sample_rate;
    bool open_raw = false;
    int opt;

    while((opt = getopt(argc, argv, "b:c:d:f:r:t:vw:")) != -1) {
        switch(opt) {
            case 'b': block_size = atoi(optarg); break;
            case 'c': required_pulses = atoi(optarg); break;
            case 'd': min_pulse_duration = atoi(optarg); break;
            case 'f': filename = optarg; break;
            case 'r': open_raw = true; sample_rate = atoi(optarg); break;
            case 't': threshold = atof(optarg); break;
            case 'v': enable_debug = true; break;
            case 'w': window_size = atoi(optarg); break;
        }
    }

    freopen(NULL, "wb", stdout);

    // Decode the m4a file to a temporary file full of floats
    FILE *fh;
    if(open_raw) {
        fh = fopen(filename, "rb");
    } else {
        fh = decode_audio_file(filename, &sample_rate);
    }
    if(fh == NULL) {
        fprintf(stderr, "There was an error opening %s\n", filename);
        return -1;
    }

    // Get block count and power
    uint32_t block_count = 0;

    fseek(fh, 0, SEEK_END);
    block_count = ceil(ftell(fh) / sizeof(float) / block_size);
    fseek(fh, 0, SEEK_SET);

    int32_t min_block_count = floor(((float)min_pulse_duration / 1000.0f) * sample_rate / block_size);
    uint32_t min_distance = (sample_rate / block_size);
    uint32_t max_distance = (sample_rate / block_size) * 15;
    debug_log("block_cnt = %i\n", block_count);
    debug_log("min_block_cnt = %i\n", min_block_count);
    debug_log("min_distance = %i\n", min_distance);
    debug_log("max_distance = %i\n", max_distance);

    // Find max power without knowing frequency
    float power[block_count];
    float max_power   = 0;
    max_power = calculate_power_of_blocks(fh, window_size, block_size, sample_rate, 0, power);
    debug_log("max_power = %f\n", max_power);

    Pulse pulse[MAX_PULSES];
    uint32_t pulse_count = 0, max_pulse_count = 0, total_pulse_count = 0;
    float best_threshold = 0.0f;
    float t_value[] = {20.0f, 19.0f, 18.0f, 17.0f, 16.0f, 15.0f, 14.0f, 13.0f, 12.0f, 11.0f, 10.0f, 9.0f, 8.0f, 7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.9f, 0.8f, 0.7f, 0.6f, 0.5f, 0.4f, 0.3f, 0.2f, 0.1f, 0.09f, 0.08f, 0.07f, 0.06f, 0.05f, 0.04f, 0.03f, 0.02f, 0.01f, 0.009f, 0.008f, 0.007f, 0.006f, 0.005f, 0.004f, 0.003f, 0.002f, 0.001f };
    int t_count = sizeof(t_value) / sizeof(float);

    for(int i = 0; i < t_count; i++) {
        float t = t_value[i];
        if(t > max_power) continue;
        debug_log("t = %f\n", t);

        // Detect pulses
        total_pulse_count = detect_pulses(power, block_count, t, min_block_count, pulse);
        if(total_pulse_count >= required_pulses) {
            pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);
            if(pulse_count >= required_pulses) {
                debug_log("pulse_count = %i\n", pulse_count);
                if(pulse_count > max_pulse_count) {
                    max_pulse_count = pulse_count;
                    best_threshold = t;
                } else if(pulse_count < max_pulse_count) {
                    break;
                }
            }
        }
    }

    total_pulse_count = detect_pulses(power, block_count, best_threshold, min_block_count, pulse);
    pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);

    // Output pulse report
    PulseSummary summary;
    if(pulse_count < required_pulses) pulse_count = 0;
    pulse_report(fh, pulse, pulse_count, total_pulse_count, sample_rate, block_size, &summary);

    // Now check for pulses with known frequency
    max_power = calculate_power_of_blocks(fh, window_size, block_size, sample_rate, summary.frequency, power);
    debug_log("max_power = %f\n", max_power);

    pulse_count = 0;
    max_pulse_count = 0;
    total_pulse_count = 0;
    best_threshold = 0.0f;
    for(int i = 0; i < t_count; i++) {
        // Set threshold
        float t = t_value[i];
        if(t > max_power) continue;
        debug_log("t = %f\n", t);

        // Detect pulses
        total_pulse_count = detect_pulses(power, block_count, t, min_block_count, pulse);
        if(total_pulse_count >= required_pulses) {
            pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);
            if(pulse_count >= required_pulses) {
                debug_log("pulse_count = %i\n", pulse_count);
                if(pulse_count > max_pulse_count) {
                    max_pulse_count = pulse_count;
                    best_threshold = t;
                } else if(pulse_count < max_pulse_count) {
                    break;
                }
            }
        }
    }

    total_pulse_count = detect_pulses(power, block_count, best_threshold, min_block_count, pulse);
    pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);

    PulseSummary summary_g;
    if(pulse_count < required_pulses) pulse_count = 0;
    pulse_report(fh, pulse, pulse_count, total_pulse_count, sample_rate, block_size, &summary_g);

    // Now sweep the frequencies
    pulse_count = 0;
    max_pulse_count = 0;
    total_pulse_count = 0;
    best_threshold = 0.0f;
    bool signal_found = false;
    for(int f = 500; f <= 4000; f += 100) {
        max_power = calculate_power_of_blocks(fh, window_size, block_size, sample_rate, (float)f, power);
        debug_log("max_power = %f at %i hz\n", max_power, f);

        for(int i = 0; i < t_count; i++) {
            // Set threshold
            float t = t_value[i];
            if(t > max_power) continue;
            debug_log("t = %f\n", t);

            // Detect pulses
            total_pulse_count = detect_pulses(power, block_count, t, min_block_count, pulse);
            if(total_pulse_count >= required_pulses) {
                pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);
                if(pulse_count >= required_pulses) {
                    debug_log("pulse_count = %i\n", pulse_count);
                    if(pulse_count > max_pulse_count) {
                        max_pulse_count = pulse_count;
                        best_threshold = t;
                    } else if(pulse_count < max_pulse_count) {
                        break;
                    }
                }
            }
        }

        total_pulse_count = detect_pulses(power, block_count, best_threshold, min_block_count, pulse);
	pulse_count = find_pulse_chain(pulse, total_pulse_count, required_pulses, MAX_PULSES, min_distance, max_distance, TOLERANCE);
        if(pulse_count >= required_pulses) break;
    }

    PulseSummary summary_s;
    if(pulse_count < required_pulses) pulse_count = 0;
    pulse_report(fh, pulse, pulse_count, total_pulse_count, sample_rate, block_size, &summary_s);

    printf("[{\"count\": %i, \"frequency\": %.1f, \"delay\": %f, \"ppm\": %f, \"temp\": %0.4f, \"confidence\": %0.2f},\n", summary.count, summary.frequency, summary.delay, summary.ppm, summary.temperature, summary.confidence);
    printf("{\"count\": %i, \"frequency\": %.1f, \"delay\": %f, \"ppm\": %f, \"temp\": %0.4f, \"confidence\": %0.2f},\n", summary_g.count, summary_g.frequency, summary_g.delay, summary_g.ppm, summary_g.temperature, summary_g.confidence);
    printf("{\"count\": %i, \"frequency\": %.1f, \"delay\": %f, \"ppm\": %f, \"temp\": %0.4f, \"confidence\": %0.2f}]\n", summary_s.count, summary_s.frequency, summary_s.delay, summary_s.ppm, summary_s.temperature, summary_s.confidence);

    //float c = fabs(summary.temperature - summary_g.temperature);
    //printf("c = %f\n", c);

    return 0;
}
