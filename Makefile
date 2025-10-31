exec=ml_breaking

CC=$(BASILISK)/qcc
# CFLAGS="-autolink -disable-dimensions -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp"
# CFLAGS= -autolink -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
# CFLAGS= -autolink -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
CFLAGS= -autolink -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -fopenmp
LIBGL= -L$(BASILISK)/gl -lglutils
INCLUDE= "$(SANDBOX)"
OPENGLIBS= -lfb_tiny
MATHLIB= -lm
DEPS= spectrum.h interpolate.h
PARAM= namlist.toml
SRC := $(exec).c 



$(exec): $(SRC) $(DEPS)
	mkdir -p $(exec)
	cp $(PARAM)  $(exec)/$(PARAM)
	$(CC) -I$(INCLUDE) $(CFLAGS) $(EVENTS) -o $(exec)/$(exec) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB) 

save:
	cp -r $(exec) "$(exec)_save"

clean:
	rm -r $(exec)

