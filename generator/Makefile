
TRGS = generator

CXXFLAGS = --std=c++11 -O3
CXXFLAGS += -Wall -Werror

all: $(TRGS)

generator: generator.cpp

clean:
	rm -f $(TRGS) *.o *.dot
	rm -f file.txt

%: %.cpp ;
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDES) $(LIBS)
