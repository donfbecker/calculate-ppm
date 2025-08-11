#!/bin/sh
ffmpeg -loglevel quiet -i $1 -f f32le -acodec pcm_f32le -ar 44100 -ac 1 -y pipe:1 | ./count-pulses -t 0.3 -l 64
