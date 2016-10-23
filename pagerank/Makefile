
TRGS = pagerank

ARGO_INSTALL_DIRECTORY=../argodsm
INCLUDES = -L${ARGO_INSTALL_DIRECTORY}/build/lib -I${ARGO_INSTALL_DIRECTORY}/src
LIBS += -largo -largobackend-mpi

CXX = mpic++

CXXFLAGS = --std=c++11 -O3
CXXFLAGS += -Wall -Werror

all: $(TRGS)

pagerank: pagerank.cpp

clean:
	rm -f $(TRGS) *.o *.dot
	rm -f file.txt

%: %.cpp ;
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDES) $(LIBS)
