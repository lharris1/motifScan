all: motifScan.cpp
	g++ -g -Wall -o motifScan motifScan.cpp

clean:
	$(RM) motifScan 
