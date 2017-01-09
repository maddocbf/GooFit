.PHONY: all clean
GOOFITDIR = ${PWD}/../../

include ../../Makefile.common

info:
	@echo "Not working yet! Use make all to test"

# AngleProgram.cpp
all:  Christoph_K1M2K892 Christoph_K1P2K892 DP4KKPiPi FinalGen0 FitKKPiPi892 FitKKPiPi FitKKPiPiDec14 GenD0892 GenD0B892 GenD0Bar GenD0B GenD0 GenDZeroMOSTUPTODATE GoodThingYouHaveMe K1M2Kstar892 K1M2Kstar K1P2Kstar892 K1P2Kstar KstarKstar892 KstarKstar KstarKstarD KstarKstarP KstarM2Kstar892 KstarM2Kstar KstarP2Kstar892 KstarP2Kstar 

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(CXXFLAGS) -c -o $@ $<

%:	%.o $(THRUSTO) $(ROOTUTILLIB) 
	$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) -o $@
	@echo "$@ done"

clean:
		@rm -f *.o $(PROGRAM) $(PROGRAM).o *.png *.root
