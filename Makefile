ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
ROOFITLIB = -lRooFitCore -lRooFit

#USERINC = -I$(shell echo $CPLUS_INCLUDE_PATH | sed "s/:/\ -I/g")
USERINC = -I$(CPLUS_INCLUDE_PATH)

CC = g++

make_exe = $(CC) $(ROOTINC) $(ROOTLIB) $(USERINC)

all: KaonTrack EffVsP

KaonTrack : KaonTrack.C
	$(make_exe) $(ROOFITLIB) $< -o $@

% : %.C
	$(make_exe) $^ -o $@


.PHONY: clean

clean:
	-rm KaonTrack EffVsP


