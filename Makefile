CC = cc
CXX = CC
AS  = as

LIBS=/home/chulwoo/THETA/CPS_Gparity/install0

#DMAPP_CFLAGS = -Wl,--whole-archive,-ldmapp,--no-whole-archive
#DMAPP_LDFLAGS = -dynamic -ldmapp -craympich-mt
#DMAPP_CFLAGS =-craympich-mt
DMAPP_LDFLAGS = -Wl,--whole-archive,-ldmapp,--no-whole-archive

CXXFLAGS=  -w -fpermissive -fopenmp -g -O3 -Wall -std=c++11  \
		${DMAPP_CFLAGS} 
#		-I${LIBS}/include 

LDFLAGS = ${DMAPP_LDFLAGS}  -fopenmp
#-L${LIBS}/lib -lqmp -mkl

clean:
	rm -f *.o *.x

all: BlockScramble.x

.SUFFIXES:
.SUFFIXES:  .o .C .S .c .x

CSRC :=$(wildcard *.c)
CCSRC :=$(wildcard *.C)
SSRC :=$(wildcard *.S)

COBJ=$(CSRC:.c=.o)
CCOBJ=$(CCSRC:.C=.o)
SOBJ=$(SSRC:.S=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))

$(BIN):  $(OBJS) $(LIBLIST)
	@echo OBJS = $(OBJS)
	$(CXX) $(OBJS) $(LIBLIST) $(LDFLAGS) -o $(BIN)

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -E $(CXXFLAGS) $(DFLAGS) $(INCLIST) $< > $*.i
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<
.S.o:
	$(AS) -o $@ $(ASFLAGS) -c $(INCLIST) $<

.o.x: $(LIBLIST)
	$(CXX) $< $(LIBLIST) $(LDFLAGS) -o $@
