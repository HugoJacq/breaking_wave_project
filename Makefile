##
#
# This make is used to compile and run a Basilisk simulation
#
#
#


EXEC=ml_breaking
OUT_DIR=out


# C Compiler options
CC=$(BASILISK)/qcc
CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp

# MPI compiler
MPICC=mpicc
MPICCFLAGS += -std=c99 -02 -g -Wall
LDFLAGS = -lgfortran -L${BASILISK}/ppr -lppr -lm
CFLAGS_MPI = -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1
CFLAGS_HPC = $(CFLAGS_MPI) -source


# todo: put source file in ./src
# and .o in ./build
# make .h real heard files, and move the code to .c


# save
# CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp
# CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp -DDISPLAY=-1



# CFLAGS="-autolink -disable-dimensions -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp"
# CFLAGS= -autolink -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
# CFLAGS= -autolink -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp


#CFLAGS_MPI = -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1


# HPC
# CFLAGS_HPC = -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -source -D_MPI=1

LIBGL= -L$(BASILISK)/gl -lglutils
INCLUDE= "$(MYSANDBOX)"
OPENGLIBS= -lfb_tiny
MATHLIB= -lm
DEPS= spectrum.h interpolate.h
PARAMETERS = namelist
PARAM= $(PARAMETERS).toml
SRC := $(EXEC).c 

TARGET = $(EXEC)/$(EXEC)
TARGET_MPI = $(EXEC)/$(EXEC)_mpi
TARGET_HPC = $(EXEC)/_$(EXEC).c
F_RESTART = $(EXEC)_restart/


# Commande pour la visualisation
# -> utilise le Makefile par défaut de Basilisk
#  CFLAGS='-I/home/jacqhugo/Debut_these/basilisk_sandbox/ -DDISPLAY=-1 -disable-dimensions' make -f ~/Debut_these/basilisk/src/Makefile.defs ml_breaking.tst

all: $(TARGET) $(EXEC)/$(PARAM)

$(EXEC)/$(PARAM): $(PARAM)
	$(info NAMLIST UPDATED !)
	@cp $(PARAM)  $(EXEC)/$(PARAM)
	

$(TARGET): $(SRC) $(DEPS)
	$(info COMPILING THE FILE $(EXEC).c:)
	@mkdir -p $(EXEC)
	$(CC) -I$(INCLUDE) $(CFLAGS) $(EVENTS) -o $(TARGET) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB)
	@chmod +x $(EXEC)/$(EXEC)

$(TARGET_MPI): $(SRC) $(DEPS)
		$(info COMPILING THE FILE $(EXEC).c FOR MPI:)
	@mkdir -p $(EXEC)
	CC99='$(MPICC) $(MPICCFLAGS)' $(CC) -I$(INCLUDE) $(CFLAGS_MPI) $(EVENTS) -o $(TARGET_MPI) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB)
	@chmod +x $(EXEC)/$(EXEC)

$(TARGET_HPC): $(SRC) $(DEPS)
		$(info COMPILING THE FILE $(EXEC).c FOR HPC:)
	@mkdir -p $(EXEC)
	$(CC) -I$(INCLUDE) $(CFLAGS_HPC) $(EVENTS) -o $(TARGET_HPC) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB)
	mv _$(EXEC).c $(TARGET_HPC)

hpc: $(TARGET_HPC)

mpi: $(TARGET_MPI)

save:
	$(info SAVING THE FOLDER)
	@cp -r $(EXEC) "$(EXEC)_save"

clean:
	rm -fr $(EXEC)
	rm -f *.nc runlog


run: $(TARGET_MPI) $(EXEC)/$(PARAM)
	(cd $(EXEC); \
		mpirun -n 16 $(EXEC)_mpi 2>&1 | /usr/bin/tee runlog; \
		mkdir -p $(OUT_DIR); \
		mv runlog out.nc u_profile.dat out ;\
		)

plot:
	cp plot.gp $(EXEC)/$(OUT_DIR)/
	(cd $(EXEC)/$(OUT_DIR); \
		gnuplot plot.gp; \
		mv *.png ../;\
		cd ..)

# restart uses 
# - ncks (with netcdf)
# - tomli (https://github.com/blinxen/tomli)
restart: $(TARGET) 
	@echo "NO YET CODED"
	@#mkdir -p $(F_RESTART)
	@#cp $(TARGET) $(F_RESTART)
	@#ncks -d time,-1 $(exec)/out.nc $(F_RESTART)/restart.nc
	@#cp $(exec)/namlist.toml $(F_RESTART)/namlist_prev.toml
	@#tomli set -f $(F_RESTART)/namlist_prev.toml restart 1 --type int > $(F_RESTART)/namlist.toml




