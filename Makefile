exec=ml_breaking

CC=$(BASILISK)/qcc
# CFLAGS="-autolink -disable-dimensions -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp"
# CFLAGS= -autolink -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
# CFLAGS= -autolink -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
#CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp
CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp -DDISPLAY=-1

LIBGL= -L$(BASILISK)/gl -lglutils
INCLUDE= "$(MYSANDBOX)"
OPENGLIBS= -lfb_tiny
MATHLIB= -lm
DEPS= spectrum.h interpolate.h
PARAM= namlist.toml
SRC := $(exec).c 

TARGET = $(exec)/$(exec)
F_RESTART = $(exec)_restart/


# Commande pour la visualisation
# -> utilise le Makefile par défaut de Basilisk
#  CFLAGS='-I/home/jacqhugo/Debut_these/basilisk_sandbox/ -DDISPLAY=-1 -disable-dimensions' make -f ~/Debut_these/basilisk/src/Makefile.defs ml_breaking.tst

all: $(TARGET) $(exec)/$(PARAM)

$(exec)/$(PARAM): $(PARAM)
	$(info NAMLIST UPDATED !)
	@cp $(PARAM)  $(exec)/$(PARAM)

$(TARGET): $(SRC) $(DEPS)
	$(info COMPILING THE FILE $(exec).c:)
	@mkdir -p $(exec)
	$(CC) -I$(INCLUDE) $(CFLAGS) $(EVENTS) -o $(exec)/$(exec) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB)
	@chmod +x $(exec)/$(exec)

save:
	$(info SAVING THE FOLDER)
	@cp -r $(exec) "$(exec)_save"

clean:
	rm -fr $(exec)
	rm -f *.nc runlog


run: $(TARGET) $(exec)/$(PARAM)
	(cd $(exec); \
		./$(exec) 2>&1 | /usr/bin/tee runlog; \
		cd ..)


plot: 
	cp plot.gp $(exec)/
	(cd $(exec); \
		gnuplot plot.gp; \
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




