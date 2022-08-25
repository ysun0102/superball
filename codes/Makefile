
ALG=REGION_then_BORDER
# ALG=BORDER_then_REGION
# ALG=WEIGHTED_ONE

# ALG=SIZE_then_STRENGTH


MACRO=-D $(ALG)
CXX=g++
CXXFLAGS=-O3 -std=c++14 -fopenmp -pthread $(MACRO)
INC=-I $(GUROBI_HOME)/include/
LIB=-L $(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi95 -lm


all : generate_$(ALG).run select_$(ALG).run

generate_$(ALG).run : src/generate.cpp
	$(CXX) $(CXXFLAGS) -m64 -g -o generate_$(ALG).run src/generate.cpp src/superball.hpp $(INC) $(LIB)

select_$(ALG).run : src/select.cpp
	$(CXX) $(CXXFLAGS) -m64 -g -o select_$(ALG).run src/select.cpp src/superball.hpp $(INC) $(LIB)


.PHONY: clean
clean : 
	rm -f *.run
