all: KaonTrack

% : %.C
	g++ `root-config --cflags --libs` -lRooFitCore -lRooFit -lMathMore $^ -o $@
