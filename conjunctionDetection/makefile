CFLAGS=
CXX=g++
CXXFLAGS=-c -Wall -O3 -std=c++11
LDFLAGS=
RM=rm -f

SOURCES = main.cpp dsyevc3.cpp dsyevv3.cpp Simulation.cpp ExternalSpaceObject.cpp SpaceObject.cpp VectorOperations.cpp FatalError.h MathWarning.h Conjunction.cpp sgp4ext.cpp sgp4unit.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=conjunctionDetection_PA_XYZ_R2_MD_fR2_CoarseSearch

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) *.o *~

dist-clean: clean
	$(RM) *~ .dependtool
