EXE=intergenerational_numeric.exe
CPP=intergenerational.cpp main.cpp
HPP=intergenerational.hpp parameters.hpp
CXX=g++
CXXFLAGS=-Wall -ggdb -std=c++20

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)
