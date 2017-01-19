# rtlmic: Wireless microphone receiver for RTL-SDR

rtlmic is a multichannel FM microphone receiver/demodulator for RTL-SDR cards.
It outputs realtime audio to JACK.

## Installation

Dependencies:

* libjack
* libvolk
* librtlsdr

To build and install, just use `make && sudo make install`.

## Usage

Basic usage is simply:

```shell
$ rtlmic [channel 1 frequency] [channel 2 frequency]...
```

You can capture as many channels as your CPU can handle, as long as they all
fit within the capture bandwidth of the RTL-SDR.

Use `rtlmic --help` to see all available options.

In order to get back correct audio, you should know certain parameters about
your microphones: the companding ratio (`-e`), the companding tau (`-E`), the
deemphasis tau (`-M`), and the FM deviation (`-w`).

Most microphones use 2:1 companding (the default ratio). The companding
tau is usually determined by an R-C filter in the microphone, connected to the
compressor chip. In my case, the component values are 10kΩ and 1µF, which gives
10kΩ * 1µF = 10ms tau. The deemphasis tau just affects the frequency response:
the higher the tau, the less high-end the microphone will have. If you do not
know the exact values for your microphone, you can just try and see what sounds
best.

By default, the RTL-SDR will be tuned to the frequency in between the highest
and lowest channel specified. However, if you have a channel near that point
(e.g. if you are only capturing one channel, or an odd number of evenly spaced
channels) then you may experience additional noise, as the RTL-SDR does not
perform well near DC. You can work around this by choosing a different center
frequency with `-f`.

You may want to play around with other parameters, e.g. adjusting the tuner gain
(`-g`) and squelch threshold (`-s`), as well as the audio gain (`-a`). If your
channels are spaced closely together, lower the transition width (`-t`) and use
a tight deviation (`-w`). Note that lowering the transition width increases CPU
usage significantly.

If you have persistent buffer over/underrun problems, you should try changing
the buffering settings `-b` and `-B`, as well as the JACK period size in the
JACK server. A few under/overruns on startup are normal, as it takes some time
for rtlmic to lock onto the exact ratio between the true SDR and JACK
frequencies.
