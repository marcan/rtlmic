PREFIX ?= /usr
CFLAGS ?= -O3

all: rtlmic

rtlmic: rtlmic.c
	gcc -Wall $(CFLAGS) -o $@ $< -lrtlsdr -lvolk -lm -ljack

install: rtlmic
	install rtlmic $(PREFIX)/bin/rtlmic

clean:
	rm -f rtlmic
