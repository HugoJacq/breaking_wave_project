
NAME=ml_convection


CFLAGS += -O1 # O2

# for local use only
ifdef MYSANDBOX
CFLAGS += -I$(MYSANDBOX)/
endif

CFLAGS += -disable-dimensions # skip dimensions check for now

all: $(NAME)/namelist.toml $(NAME).tst plot

plot: $(NAME)/plots

$(NAME).tst: CC = mpicc -D_MPI=16

$(NAME)/namelist.toml: namelist.toml
	mkdir -p $(NAME)
	ln -sf --target-directory=$(NAME) ../namelist.toml

_$(NAME).c: CFLAGS += -D_MPI=1

hpc: _$(NAME).c


plot: $(NAME)/plots


include $(BASILISK)/Makefile.defs

