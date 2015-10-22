all: KaonTrack EffVsP

% : %.C
	g++ `root-config --cflags --libs` -lRooFitCore -lRooFit -lMathMore $^ -o $@
