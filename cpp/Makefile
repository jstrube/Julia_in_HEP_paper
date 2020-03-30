CPPFLAGS = -Wall \
        -I$(LCIO)/include \
        $(shell root-config --cflags)

LDFLAGS = -lm \
        -L$(LCIO)/lib \
	-L . \
        -llcio \
        -lsio \
        $(shell root-config --libs)


all: libFoxWo fjcoreFoxWoAnalysis
	time g++ -o cppAnalysis cppAnalysis.cxx $(CPPFLAGS) $(LDFLAGS) 
	g++ -o cppFoxWo cppFoxWo.cxx -lFoxWo $(CPPFLAGS) $(LDFLAGS) 
	g++ -o fjcoreAnalysis fjcoreAnalysis.cxx fjcore.cc $(CPPFLAGS) $(LDFLAGS) 
	
foxwo: cppFoxWo.cxx
	g++ -o cppFoxWo cppFoxWo.cxx -lFoxWo $(CPPFLAGS) $(LDFLAGS) 

libFoxWo: libFoxWo.cc libFoxWo.hh
	g++ -o libFoxWo.so -c --shared libFoxWo.cc $(CPPFLAGS) $(LDFLAGS)

fjcoreFoxWoAnalysis: fjcore_FoxWoAnalysis.cxx libFoxWo.hh
	g++ -o fjcoreFoxWoAnalysis fjcore_FoxWoAnalysis.cxx fjcore.cc -lFoxWo $(CPPFLAGS) $(LDFLAGS) 

