/*
Copyright (C) 2017 Hector Martin "marcan" <marcan@marcan.st>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include <rtl-sdr.h>
#include <complex.h>
#include <volk/volk.h>
#include <jack/jack.h>
#include <jack/ringbuffer.h>

struct onepole {
    float acc;
    float alpha;
};

struct channel {
    int f;
    int df;

    int ntaps;
    float complex **taps;

    float complex phase;
    float complex phi;
    float complex *sample;
    float complex last;

    struct onepole squelch_lpf;
    float squelch_ctr;
    struct onepole dc_lpf;
    struct onepole env_lpf;
    struct onepole emph_lpf;

    float acc;
    int cnt;
    float *audio_buf1;
    float *p1;

    jack_ringbuffer_t *rb;
    int running;
    jack_port_t *port;
    float history[3];
    float mu;
};

struct state {
    int verbose;

    // Sampling
    int tuner_gain;         // Tuner gain in .1dB increments
    rtlsdr_dev_t *dev;
    int fs;                 // Sampling frequency
    int fc;                 // Center tuning frequency

    // Buffering
    int blockcnt;           // Number of outstanding USB transactions
    int blocksz;            // Buffer size in samples (raw bytes / 2)
    int bufpad;             // Number of history samples to keep
    int bufsz;              // Total buffer size
    float complex *buffer;
    int p;                  // Current buffer pointer

    // Filtering
    struct onepole dc_lpf_i;
    struct onepole dc_lpf_q;
    int width;              // Filter half-width (deviation)
    int twidth;             // Transition band width
    int decimation;         // Decimation after filtering
    float fi;               // Intermediate frequency (fs / decimation)
    int ntaps;              // Number of filter taps
    int align;              // Required platform alignment for buffers in bytes
    int nalign;             // Alignment in complex samples

    // Audio processing
    float fa;               // Audio sampling frequency
    int adecimation;        // Decimation before audio processing
    float squelch_thr;      // Squelch threshold
    float squelch_tau;      // Squelch LPF tau
    float dc_tau;           // DC removal LPF tau
    float exp_tau;          // Expander tau
    float exp_ratio;        // Expander ratio
    float emph_tau;         // De-emphasis tau
    float audio_gain;       // Audio gain factor

    // Audio buffering/resampling
    int fjack;              // JACK sample rate
    int rb_size;            // JACK ringbuffer size
    int iblock;             // Input (SDR) block size at AF frequency;
    int pblock;             // Output (JACK) block size at AF frequency;
    int oblock;             // Output (JACK) block size
    int low_thresh;         // Ring buffer low threshold
    int high_thresh;        // Ring buffer high threshold
    float dmu;              // Downsampling factor
    struct onepole buf_lpf; // Resampling loop lpf

    // Channels
    int nch;
    struct channel *ch;

    jack_client_t *client;
    int jack_alive;
};

struct state st;

static int compute_ntaps(int transition_width, int sample_rate)
{
    float delta_f = transition_width / (float)sample_rate;
    return (((int)(3.3f / delta_f + 0.5f)) & ~1) + 1;
}

static void generate_lpf(float complex *taps, int ntaps, int fs, int width)
{
    int M = ntaps / 2;
    double fwT0 = 2 * M_PI * width / (double)fs;

    // Sinc function with Hamming window
    // Taps are real here but stored as complex for rotation later
    taps[M] = fwT0;
    float gain = taps[M];
    for (int i = 1; i <= M; i++) {
        float c = sin(i * fwT0) / i * (0.54 + 0.46 * cos(M_PI * i / (double)M));
        gain += 2 * c;
        taps[M-i] = taps[M+i] = c;
    }

    // Normalize
    gain = 1.f / gain;
    for(int i = 0; i < ntaps; i++)
        taps[i] *= gain;
}

static void init_onepole(struct onepole *filt, float fs, float tau)
{
    if (tau == 0)
        filt->alpha = 1.0;
    else
        filt->alpha = 1.0 - expf(-(1.0 / fs) / tau);
    filt->acc = 0;
}

static inline float filter_onepole(struct onepole *filt, float v)
{
    filt->acc = filt->acc * (1.f - filt->alpha) + v * filt->alpha;
    return filt->acc;
}

static inline float cubic(float *y, float mu)
{
   float a0,a1,a2,mu2;

   mu2 = mu * mu;
   a0 = y[3] - y[2] - y[0] + y[1];
   a1 = y[0] - y[1] - a0;
   a2 = y[2] - y[0];

   return a0 * mu * mu2 + a1 * mu2 + a2 * mu + y[1];
}

static void init_channel(int c)
{
    struct channel *ch = &st.ch[c];

    ch->df = st.fc - ch->f;
    ch->ntaps = st.ntaps;

    ch->taps = malloc(st.nalign * sizeof(float complex *));
    for (int i = 0; i < st.nalign; i++) {
        ch->taps[i] = volk_malloc(
            sizeof(float complex) * (ch->ntaps + st.nalign), st.align);
        memset(ch->taps[i], 0, sizeof(float complex) * ch->ntaps);
    }

    generate_lpf(ch->taps[0], st.ntaps, st.fs, st.width);

    // Shift LPF into BPF
    double fwT0 = 2 * M_PI * ch->df / (double)st.fs;
    for (int i = 0; i < st.ntaps; i++) {
        ch->taps[0][i] = ch->taps[0][i] * cexpf(lv_cmake(0, i * fwT0));
    }
    // Frequency shifter
    ch->phase = lv_cmake(1.f, 0.f);
    ch->phi = cexp(lv_cmake(0, fwT0 * st.decimation));

    // Build offset tap vectors, padded with zeros, for SIMD alignment
    for (int i = 1; i < st.nalign; i++)
        memcpy(&ch->taps[i][i], ch->taps[0], st.ntaps * sizeof(float complex));

    // Make sure output sample buffer is aligned too
    ch->sample = volk_malloc(sizeof(float complex), st.align);

    init_onepole(&ch->squelch_lpf, st.fi, st.squelch_tau);
    init_onepole(&ch->dc_lpf, st.fa, st.dc_tau);
    init_onepole(&ch->env_lpf, st.fa, st.exp_tau);
    init_onepole(&ch->emph_lpf, st.fa, st.emph_tau);

    ch->audio_buf1 = malloc(st.iblock * sizeof(float));
    ch->p1 = ch->audio_buf1;

    ch->rb = jack_ringbuffer_create(st.rb_size * sizeof(float));

    char name[16];
    sprintf(name, "channel_%d", c + 1);
    ch->port = jack_port_register(st.client, name, JACK_DEFAULT_AUDIO_TYPE,
                                  JackPortIsOutput, 0);
}

static void process_channel_sample(int c, float complex *buf)
{
    struct channel *ch = &st.ch[c];

    uintptr_t pad = (((uintptr_t)buf) & (st.align - 1)) / sizeof(float complex);

    // Bandpass filter
    volk_32fc_x2_dot_prod_32fc_a(ch->sample, buf - pad, ch->taps[pad],
                                 st.ntaps + pad);
    // Shift frequency
    float complex ss = *ch->sample * ch->phase;
    ch->phase *= ch->phi;
    // Demodulate FM
    float s = cargf(lv_conj(ch->last) * ss);
    // Squelch
    if (filter_onepole(&ch->squelch_lpf, cabsf(ss)) < st.squelch_thr) {
        s = 0;
        ch->squelch_ctr = 0;
    } else if (ch->squelch_ctr < 1.0) {
        s *= ch->squelch_ctr;
        ch->squelch_ctr += 0.00005;
    }
    // Simple boxcar filter for audio decimation
    ch->acc += s;
    if (++ch->cnt >= st.adecimation) {
        *(ch->p1++) = ch->acc / st.adecimation;
        ch->cnt = 0;
        ch->acc = 0;
    }

    ch->last = ss;
}

static void process_channel_audio(int c)
{
    struct channel *ch = &st.ch[c];

    jack_ringbuffer_data_t wdata[2];

    jack_ringbuffer_get_write_vector(ch->rb, wdata);

    float *dp = (void*)wdata[0].buf;
    int wcnt = wdata[0].len / sizeof(float);
    int put = 0;

    float *p = ch->audio_buf1;
    while (p != ch->p1) {
        if (!wcnt) {
            if (wdata[1].len) {
                dp = (void*)wdata[1].buf;
                wcnt = wdata[1].len / sizeof(float);
                wdata[1].len = 0;
            } else {
                if (c == 0)
                    fprintf(stderr, "Ring buffer overrun! Left=%d           \n",
                            (int)(ch->p1 - p));
                break;
            }
        }

        float s = *p++;
        // Remove DC
        s -= filter_onepole(&ch->dc_lpf, s);
        // Expander
        s *= powf(filter_onepole(&ch->env_lpf, fabsf(s)), st.exp_ratio - 1.f);
        // De-emphasis
        s = filter_onepole(&ch->emph_lpf, s);
        // Gain
        s *= st.audio_gain;
        *dp++ = s;
        wcnt--;
        put++;
    }
    jack_ringbuffer_write_advance(ch->rb, put * sizeof(float));
    ch->p1 = ch->audio_buf1;

    // Normalize phase to make sure it doesn't go wacky
    ch->phase /= cabsf(ch->phase);
}


void got_samples(unsigned char *buf, uint32_t len, void *ctx)
{
    if (!st.jack_alive)
        return;

    if (len != st.blocksz * 2) {
        printf("Got %d bytes, expected %d!\n", len, st.blocksz);
        exit(1);
    }

    memcpy(st.buffer, &st.buffer[st.blocksz],
           st.bufpad * sizeof(float complex));

    float complex *dst = &st.buffer[st.ntaps];
    len /= 2;

    while (len--) {
        float i = (buf[0] - 127) / 127.f;
        float q = (buf[1] - 127) / 127.f;
        i -= filter_onepole(&st.dc_lpf_i, i);
        q -= filter_onepole(&st.dc_lpf_q, q);
        *dst++ = lv_cmake(i, q);
        buf += 2;
    }

    int p = st.p;
    int max = st.blocksz + st.decimation;
    float complex *pbuf = &st.buffer[p];
    for (; p < max; p += st.decimation)
    {
        for (int i = 0; i < st.nch; i++)
            process_channel_sample(i, pbuf);
        pbuf += st.decimation;
    }
    st.p = p - st.blocksz;

    for (int i = 0; i < st.nch; i++) {
        process_channel_audio(i);
    }

    float fullness = jack_ringbuffer_read_space(st.ch[0].rb) / 
            sizeof(float) / (float)st.rb_size;

    float f2 = filter_onepole(&st.buf_lpf, fullness);

    st.dmu = st.fa / st.fjack * (1.0 + 0.2 * (f2 - 0.5));

    if (st.verbose >= 2)
        fprintf(stderr, "RB fullness: %.02f (%.02f) AF=%.01f     \r",
                100.0f * fullness, 100.0f * f2, st.dmu * st.fjack);

}

void jack_shutdown (void *arg)
{
    exit (1);
}

int jack_process (jack_nframes_t nframes, void *arg)
{
    st.jack_alive = 1;
    for (int c = 0; c < st.nch; c++) {
        struct channel *ch = &st.ch[c];
        float *o = (float *)jack_port_get_buffer(ch->port, nframes);
        jack_ringbuffer_data_t rvec[2];
        jack_ringbuffer_get_read_vector(ch->rb, rvec);
        int avail = (rvec[0].len + rvec[1].len) / sizeof(float);
        int left = nframes;
        int read = 0;

        if (!ch->running && avail < (st.rb_size / 2)) {
            memset(o, 0, nframes * sizeof(float));
            continue;
        }

        avail = rvec[0].len / sizeof(float);
        float *r = (void *)rvec[0].buf;
        ch->running = 1;

        float y[4], mu=ch->mu, dmu=st.dmu;
        memcpy(y, ch->history, sizeof(float) * 4);

        while (left--) {
            *o++ = cubic(y, ch->mu);
            mu += dmu;
            while (mu >= 1.f) {
                y[0] = y[1];
                y[1] = y[2];
                y[2] = y[3];
                if (!avail) {
                    if (rvec[1].len) {
                        avail = rvec[1].len / sizeof(float);
                        rvec[1].len = 0;
                        r = (void *)rvec[1].buf;
                    } else {
                        if (c == 0)
                            fprintf(stderr,
                                    "Ring buffer underrun! left=%d          \n",
                                    left);
                        memset(o, 0, sizeof(float) * left);
                        y[3] = 0;
                        while (mu >= 1.f)
                            mu -= 1.f;
                        ch->running = 0;
                        goto next_channel;
                    }
                }
                y[3] = *r++;
                mu -= 1.f;
                read++;
                avail--;
            }
        }
next_channel:
        jack_ringbuffer_read_advance(ch->rb, read * sizeof(float));
        ch->mu = mu;
        memcpy(ch->history, y, sizeof(float) * 4);
    }

    return 0;
}

static struct option long_options[] =
{
    {"verbose",             no_argument,       0, 'v'},
    {"client-name",         required_argument, 0, 'c'},
    {"fc",                  required_argument, 0, 'f'},
    {"rate",                required_argument, 0, 'r'},
    {"tuner-gain",          required_argument, 0, 'g'},
    {"blocksize",           required_argument, 0, 'b'},
    {"blocks",              required_argument, 0, 'n'},
    {"width",               required_argument, 0, 'w'},
    {"transition-width",    required_argument, 0, 't'},
    {"squelch",             required_argument, 0, 's'},
    {"squelch-tau",         required_argument, 0, 'S'},
    {"dc-tau",              required_argument, 0, 'D'},
    {"expander-ratio",      required_argument, 0, 'e'},
    {"expander-tau",        required_argument, 0, 'E'},
    {"deemph-tau",          required_argument, 0, 'M'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: rtlmic [OPTION]... [FREQUENCY]...\n");
    fprintf(stderr, "rtlmic - Demodulate FM microphones using an RTL-SDR\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -h, --help             this help\n");
    fprintf(stderr, "  -v, --verbose          be verbose\n");
    fprintf(stderr, "  -c, --client-name=NAME JACK client name (default: rtlmic)\n");
    fprintf(stderr, "  -f, --fc=HZ            center frequency to tune to\n");
    fprintf(stderr, "                           (default: auto)\n");
    fprintf(stderr, "  -r, --rate=HZ          capture sample rate\n");
    fprintf(stderr, "                           (default: auto)\n");
    fprintf(stderr, "  -g, --tuner-gain=DB    tuner gain in dB\n");
    fprintf(stderr, "                           (default: 10 dB)\n");
    fprintf(stderr, "  -b, --blocksize=SIZE   capture block size in samples\n");
    fprintf(stderr, "                           must be a multiple of 256\n");
    fprintf(stderr, "                           (default: 8192)\n");
    fprintf(stderr, "  -n, --blocks=NUMBER    number of outstanding transfers\n");
    fprintf(stderr, "                           (default: 4)\n");
    fprintf(stderr, "  -w, --width=HZ         filter half-width (FM deviation)\n");
    fprintf(stderr, "                           (default: 100000 Hz)\n");
    fprintf(stderr, "  -t, --transition-width=HZ\n");
    fprintf(stderr, "                         filter transition bandwidth\n");
    fprintf(stderr, "                           (default: 50000 Hz)\n");
    fprintf(stderr, "  -s, --squelch=DB       squelch level in dB\n");
    fprintf(stderr, "                           (default: -50 dB)\n");
    fprintf(stderr, "  -S, --squelch-tau=MS   squelch time constant in msec\n");
    fprintf(stderr, "                           (default: 0.1 ms)\n");
    fprintf(stderr, "  -D, --dc-tau=S         DC removal time constant in sec\n");
    fprintf(stderr, "                           (default: 2 s)\n");
    fprintf(stderr, "  -e, --expander-ratio=N expander ratio. 1=1:1, 2=2:1, etc.\n");
    fprintf(stderr, "                           (default: 2 (2:1))\n");
    fprintf(stderr, "  -E, --expander-tau=MS  expander time constant in msec\n");
    fprintf(stderr, "                           (default: 10 ms)\n");
    fprintf(stderr, "  -M, --deemph-tau=US    de-emphasis time constant in usec\n");
    fprintf(stderr, "                           (default: 75 us)\n");
    fprintf(stderr, "  -a, --audio-gain=DB    audio gain in dB (default: 0 dB)\n");
}

// "Nice" sample rates for RTL-SDR
int sample_rates[] = {
    240000, 300000, 960000, 1152000, 1200000, 1440000,
    1600000, 1800000, 1920000, 2400000, 2880000, 3200000,
    0
};

int main(int argc, char **argv)
{
    memset(&st, 0, sizeof(st));

    char *client_name = "rtlmic";

    int fmin = INT_MAX;
    int fmax = 0;

    st.tuner_gain = 100;
    st.blocksz = 8192;
    st.blockcnt = 4;
    st.width = 100000;
    st.twidth = 50000; 
    st.squelch_thr = powf(10.f, -50 / 20.f);
    st.squelch_tau = 0.0001f;
    st.dc_tau = 2.f;
    st.exp_tau = 0.01f;
    st.exp_ratio = 2.f;
    st.emph_tau = 75e-6;
    st.audio_gain = 1.f;

    while (1) {
        int c = getopt_long(argc, argv, "vc:f:r:g:b:n:w:t:s:S:D:e:E:M:a:",
                            long_options, NULL);
        if (c == -1)
            break;

        switch (c)
        {
            case 'v':
                st.verbose++;
                break;
            case 'c':
                client_name = optarg;
                break;
            case 'f':
                st.fc = atoi(optarg);
                break;
            case 'r':
                st.fs = atoi(optarg);
                break;
            case 'g':
                st.tuner_gain = atoi(optarg) * 10;
                break;
            case 'b':
                st.blocksz = atoi(optarg);
                break;
            case 'n':
                st.blockcnt = atoi(optarg);
                break;
            case 'w':
                st.width = atoi(optarg);
                break;
            case 't':
                st.twidth = atoi(optarg);
                break;
            case 's':
                st.squelch_thr = powf(10.f, atoi(optarg) / 20.f);
                break;
            case 'S':
                st.squelch_tau = atof(optarg) / 1000.f;
                break;
            case 'D':
                st.dc_tau = atof(optarg);
                break;
            case 'e':
                st.exp_ratio = atof(optarg);
                break;
            case 'E':
                st.exp_tau = atof(optarg) / 1000.0;
                break;
            case 'M':
                st.emph_tau = atof(optarg) / 1000000.0;
                break;
            case 'a':
                st.audio_gain = powf(10.f, atoi(optarg) / 20.f);
                break;
            case 'h':
                usage();
                return 1;
            default:
                usage();
                return 1;
        }
    }

    if (optind == argc)
    {
        fprintf(stderr, "No frequencies specified.\n");
        usage();
        return 1;
    }

    st.ch = malloc((argc - optind) * sizeof(struct channel));
    memset(st.ch, 0, (argc - optind) * sizeof(struct channel));

    while (optind < argc)
    {
        int f = atoi(argv[optind++]);
        st.ch[st.nch].f = f;
        if (f < fmin)
            fmin = f;
        if (f > fmax)
            fmax = f;
        if (st.verbose)
            fprintf(stderr, "Channel %d: %d Hz\n", st.nch, f);
        st.nch++;
    }
    //st.nch = 1;

    st.align = volk_get_alignment();
    st.nalign = st.align / sizeof(float complex);
    if (st.verbose)
        fprintf(stderr, "Platform alignment: %d (%d samples)\n",
                st.align, st.nalign);

    if (st.fc == 0) {
        st.fc = (fmin + fmax) / 2;
    }
    if (st.verbose)
        fprintf(stderr, "Center frequency: %d Hz\n", st.fc);

    if (st.fs == 0) {
        // Give us one extra filter-width worth of padding
        int need_rate = (fmax - fmin) + 3 * st.width;
        for (int i = 0; st.fs < need_rate; i++) {
            if (!sample_rates[i]) {
                fprintf(stderr, "Channels are too far apart! "
                                "Total bandwidth needed is %d Hz\n", need_rate);
                return 1;
            }
            st.fs = sample_rates[i];
        }
    }
    if (st.verbose) {
        fprintf(stderr, "RTL-SDR sample rate: %d Hz\n", st.fs);
        fprintf(stderr, "RTL-SDR block size: %d\n", st.blocksz);
    }

    jack_status_t jack_status;	

    st.client = jack_client_open(client_name, JackNullOption, &jack_status);
    if (!st.client) {
        fprintf (stderr, "jack server not running?\n");
        return 1;
    }
    st.fjack = jack_get_sample_rate(st.client);
    st.oblock = jack_get_buffer_size(st.client);
    if (st.verbose) {
        fprintf(stderr, "JACK sample rate: %d Hz\n", st.fjack);
    }

    st.decimation = st.fs / (st.width * 2);
    st.fi = (float)st.fs / st.decimation;
    st.adecimation = (int)st.fi / st.fjack;
    st.fa = st.fi / st.adecimation;

    st.ntaps = compute_ntaps(st.twidth, st.fs);

    st.bufpad = st.ntaps + st.decimation;
    st.bufsz = st.blocksz + st.bufpad;
    st.buffer = volk_malloc(sizeof(float complex) * st.bufsz, st.align);
    st.p = st.decimation;
    memset(st.buffer, 0, sizeof(float complex) * st.bufsz);

    st.iblock = (st.blocksz / st.decimation + 1) / st.adecimation + 1;
    st.pblock = (int)((float)st.oblock / st.fjack * st.fa + 5);

    st.low_thresh = st.iblock;
    if (st.pblock < st.iblock)
        st.low_thresh = st.pblock;

    st.low_thresh += st.low_thresh / 3;

    st.rb_size = 128;
    if (st.rb_size < 3 * st.iblock)
        st.rb_size = 3 * st.iblock;
    if (st.rb_size < 3 * st.oblock)
        st.rb_size = 3 * st.oblock;

    st.high_thresh = st.rb_size - st.low_thresh;

    st.dmu = st.fa / st.fjack;

    if (st.verbose) {
        fprintf(stderr, "Decimation: %d\n", st.decimation);
        fprintf(stderr, "IF: %.2f Hz\n", st.fi);
        fprintf(stderr, "Audio decimation: %d\n", st.adecimation);
        fprintf(stderr, "AF: %.2f Hz\n", st.fa);
        fprintf(stderr, "Filter taps: %d\n", st.ntaps);
        fprintf(stderr, "Buffer sizes: %d/%d/%d -> %d\n",
                st.iblock, st.pblock, st.oblock, st.rb_size);
        fprintf(stderr, "Buffer low threshold: %d\n", st.low_thresh);
    }

    init_onepole(&st.dc_lpf_i, st.fs, 1.);
    init_onepole(&st.dc_lpf_q, st.fs, 1.);

    init_onepole(&st.buf_lpf, 1., 200.);
    st.buf_lpf.acc = 0.5;

    for (int i = 0; i < st.nch; i++)
        init_channel(i);

    int r = rtlsdr_open(&st.dev, 0);
    if (r < 0) {
        fprintf(stderr, "Failed to open rtlsdr device: error %d.\n", r);
        return 1;
    }

    rtlsdr_set_offset_tuning(st.dev, 1);
    rtlsdr_set_tuner_gain(st.dev, st.tuner_gain);
    rtlsdr_set_agc_mode(st.dev, 0);
    rtlsdr_reset_buffer(st.dev);
    rtlsdr_set_center_freq(st.dev, st.fc);
    rtlsdr_set_sample_rate(st.dev, st.fs);

    jack_set_process_callback(st.client, jack_process, NULL);
    jack_on_shutdown(st.client, jack_shutdown, NULL);
    if (jack_activate(st.client)) {
        fprintf(stderr, "cannot activate client");
        return 1;
    }

    rtlsdr_read_async(st.dev, got_samples, NULL, st.blockcnt, st.blocksz * 2);
    jack_client_close(st.client);
}
