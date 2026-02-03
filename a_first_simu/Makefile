CFLAGS += -O1 # O2

# for local use only
ifdef MYSANDBOX
CFLAGS += -I$(MYSANDBOX)/
endif

CFLAGS += -disable-dimensions # skip dimensions check for now
LIBS += -lgsl

all: ml_breaking/namelist.toml ml_breaking.tst plot

plot: ml_breaking/plots

ml_breaking.tst: CC = mpicc -D_MPI=16

ml_breaking/namelist.toml: namelist.toml
	mkdir -p ml_breaking
	ln -sf --target-directory=ml_breaking ../namelist.toml


_ml_breaking.c: CFLAGS += -D_MPI=1
#_ml_breaking.c: ml_breaking.c spectrum.h interpolate.h



plot: ml_breaking/plots


include $(BASILISK)/Makefile.defs

