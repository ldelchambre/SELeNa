CC=gcc

# -std=c99             : Source comply with the ISO C99 standard
# -D_REENTRANT         : Use reentrant version of the functions since we will use threads
# -D_XOPEN_SOURCE=600  : Use th POSIX.1-2001 standard
# -D_GNU_SOURCE        : GNU compliant source
# -fPIC \              : Generate position independant code
# -Wall \              : Print a lot of warning messages
# -Wno-unused-parameter: Except the unused function parameters
# -g \                 : Include debugging information within the code
# -O0                  : Do not optimize the code (allow an easier debugging)
CFLAGS=-std=c99 \
       -D_REENTRANT \
       -D_XOPEN_SOURCE=600 \
       -D_GNU_SOURCE \
       -fPIC \
       -W \
       -Wall \
       -Wno-unused-parameter \
       -g \
       -O0

INCLUDES=-I../includes/

LDFLAGS=

LDLIBS=-lm -lpthread

SRC= astrostr.c\
     combination.c \
     circstat.c \
     extratree.c \
     ExtraTrees.c \
     fload.c \
     fparse.c \
     gl.c \
     graddlsqr.c \
     htm.c \
     htmdb.c \
     information.c \
     nelder.c \
     nsieg.c \
     pheap.c \
     raster.c \
     rfprintf.c \
     stat.c \
     threadpool.c \
     tree.c \
     treealgorithm.c \
     treecommittee.c \
     utils.c

OBJ= $(SRC:.c=.o)

TESTCASES= astrostr.test \
           combination.test \
           circstat.test\
           gl.test \
           graddlsqr.test \
           htm.test \
           htmdb.test \
           nelder.test \
           nsieg.test \
           pheap.test \
           raster.test \
           threadpool.test \
           tree.test \
           utils.test

TESTCASES_OBJ= $(TESTCASES:.test=_test.o)

runtest = echo "Running $(1)";./$(1);
runmemcheck = echo "Memcheck $(1)";valgrind --leak-check=full --trace-children=yes --track-origins=yes ./$(1);

.PHONY: clean dist-clean

all: $(OBJ)

%.o : %.c
	@echo "CC $^"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

clean:
	@echo "Cleaning object files"
	@rm -rf $(TESTCASES_OBJ) $(OBJ)

clean-test:
	@echo "Cleaning testcases"
	@rm -rf $(TESTCASES)

dist-clean: clean clean-test

install: all

%.test : $(OBJ) %_test.c
	@echo "MKTEST $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

testcases: $(TESTCASES)

check: testcases
	@$(foreach prog,$(TESTCASES),$(call runtest,$(prog)))

memcheck: testcases
	@$(foreach prog,$(TESTCASES),$(call runmemcheck,$(prog)))
