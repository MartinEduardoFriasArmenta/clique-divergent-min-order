# Requires nauty (Debian/Ubuntu: apt install nauty libnauty-dev; macOS: brew install nauty).
NAUTY_INC ?= -I/usr/include/nauty -I/usr/include/x86_64-linux-gnu/nauty
NAUTY_LIB ?= -lnautyL1
CFLAGS = -O3 -march=native -Wall -Wno-misleading-indentation
all: gen_disks gen_disks_ca kdisk
gen_disks: gen_disks.c
	gcc $(CFLAGS) $(NAUTY_INC) -DWORDSIZE=64 -DMAXN=WORDSIZE -o $@ $< $(NAUTY_LIB)
gen_disks_ca: gen_disks_ca.c
	gcc $(CFLAGS) $(NAUTY_INC) -DWORDSIZE=64 -DMAXN=WORDSIZE -o $@ $< $(NAUTY_LIB)
kdisk: kdisk.c
	gcc $(CFLAGS) -o $@ $<
clean:
	rm -f gen_disks gen_disks_ca kdisk
