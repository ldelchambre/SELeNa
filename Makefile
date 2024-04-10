CC=gcc
CFLAGS=-D_GNU_SOURCE -D_REENTRANT -g -O0 -W -Wall -Wno-unused-parameter
BASE_DIR=src

INCLUDES=-I./includes
LDFLAGS=
LDLIBS=-lm -lpthread
EXEC= cluster contaminant dataset density ert htmid pconvert pgrid simulation split xmatch

.PHONY: clean dist-clean

all: $(EXEC)

base:
	@$(MAKE) -C $(BASE_DIR)

clean:
	@rm -rf $(OBJ)
	@make -C $(BASE_DIR) dist-clean

dist-clean: clean
	@rm -rf $(EXEC)

%.o : %.c
	@echo "CC $^"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

cluster: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ cluster.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

contaminant: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ contaminant.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

dataset: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ dataset.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

density: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ density.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

ert: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ ert.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

htmid : base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ htmid.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

pconvert: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ pconvert.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

pgrid: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ pgrid.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

split: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ split.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

simulation: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ simulation.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)

xmatch: base
	@echo "MK $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ xmatch.c $(BASE_DIR)/*.o $(LDFLAGS) $(LDLIBS)
