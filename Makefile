# Set C++ COMPILER
CXX = g++
FLAGS = -O3 --std=c++17 
EXEC = pvcreator
INST_DIR = /usr/local/bin

# SET INCLUDE AND LIB DIRECTORY
CFITSIOINC = /usr/local/include
CFITSIOLIB = /usr/local/lib
WCSLIB = /usr/local/include
WCSINC = /usr/local/lib

LIBS = -L$(CFITSIOLIB) -L$(WCSLIB) -lwcs -lcfitsio
INCS = -I./src -I$(CFITSIOINC) -I$(WCSINC) 

pvcreator:
	$(CXX) $(FLAGS) -o $(EXEC) src/pvcreator.cpp $(LIBS) $(INCS)

install:
	mv pvcreator $(INST_DIR)

clean:
	rm -rf pvcreator
