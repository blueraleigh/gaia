TSKFLAGS = $(C99FLAGS) $(CPICFLAGS)

TSKFLAGS += -I.

TSKFLAGS += -I./subprojects/kastore

TSKOBJ = kastore.o \
	tskit_convert.o \
	tskit_core.o \
	tskit_genotypes.o \
	tskit_haplotype_matching.o \
	tskit_stats.o \
	tskit_tables.o \
	tskit_trees.o 

libtskit.a: $(TSKOBJ)
	$(AR) rcs $@ $(TSKOBJ)

kastore.o:
	$(CC99) -c $(TSKFLAGS) ./subprojects/kastore/kastore.c -o kastore.o

tskit_%.o:
	$(CC99) -c $(TSKFLAGS) ./tskit/$*.c -o $@
