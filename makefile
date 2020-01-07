CC=gcc
OPT=-O2
GPROF=#-pg
GDB=#-g
W= -Wall -Wextra \
   -Wno-sign-compare\
   -Wno-pointer-sign\
   -Wformat=0\
   -Warray-bounds\
   -Wfloat-equal\
   -Wimplicit\
   -Wmaybe-uninitialized\
   -Wmissing-braces\
   -Wparentheses\
   -Wsequence-point\
   -Wtype-limits\
   -Wundef\
   -Wuninitialized\
   -Wunused\
   -Wunused-but-set-variable\
   -Wunused-parameter\
   -Wempty-body\
   -Wmemset-elt-size\
   -Wduplicated-branches\
   -Wswitch-unreachable\
   -Winline\

STDFLAG= -std=gnu11   # -std=gnu99 can be used
CFLAGS= -c $(STDFLAG) -MMD $(OPT) $(GPROF) $(W) $(GDB) $(ASM)
OFLAGS= -lm $(GPROF)
INCL= -I$(SRCDIR)/mol -I$(SRCDIR)/math -I$(SRCDIR)/q -I$(SRCDIR)/qap


OBJDIR=./obj
SRCDIR=./src

molsrc=$(wildcard $(SRCDIR)/mol/*.c)
molobj=$(molsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

mathsrc=$(wildcard $(SRCDIR)/math/*.c)
mathobj=$(mathsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

qsrc=$(wildcard $(SRCDIR)/q/*.c)
qobj=$(qsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

qapsrc=$(wildcard $(SRCDIR)/qap/*.c)
qapobj=$(qapsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)


all : q qap

q : $(molobj) $(mathobj) $(qobj) $(OBJDIR)/q.o
	$(CC) $^ -o $@ $(OFLAGS)
qap : $(molobj) $(mathobj) $(qobj) $(qapobj) $(OBJDIR)/qap.o
	$(CC) $^ -o $@ $(OFLAGS) -lpthread

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $< -o $@ $(INCL)

clean :
	rm -f $(OBJDIR)/*/*.o $(OBJDIR)/*.o q qap
cleand :
	rm -f $(OBJDIR)/*/*.d $(OBJDIR)/*.d


include $(wildcard $(OBJDIR)/*/*.d $(OBJDIR)/*.d)

