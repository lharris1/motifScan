/*
 motifScan: TODO add a program descripton here!!!

 Created by Lincoln Harris
 Swarthmore College
 June 2016

*/

#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <ctype.h>
#include <stdio.h> 

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//////////////////////////MOTIFSCAN IMPLEMENTATION/////////////////////////////
///////////////////////////////////////////////////////////////////////////////

string fileToString(string fileName) {

	//TODO: IMPLEMENT ME!!
	return "dummy";
}

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize) {
	
	//vector<string>* myF1 = fileToString(file1);
	//vector<string>* myF2 = fileToString(file2); 
	//cout << file1; 
	//cout << file2; 
	//cout <<  mot1; 
	//cout <<  mot2; 
	//cout << " ";
	//cout <<  winSize << endl << endl; 

	

	return 1; //default, if no motif conservation found; else, return seq location where conserved motifs found (on ref seq)
}

int main(int argc, char* argv[]) {

	int location; 

	if(argc != 6) {
		cout << endl; 
		cout << "ERROR" << endl;
		cout << "usage: ./motifScan file1.fasta file2.fasta [motif1] [motif2] [windowSize]" << endl;
		cout << "See README file for more" << endl << endl; 
	}
	else {
		string f1 = string(argv[1]);
		string f2 = string(argv[2]);
		string m1 = string(argv[3]);
		string m2 = string(argv[4]); 
		int wSize = atoi(argv[5]); 
		if(wSize == 0) {
			cout << endl; 
			cout << "ERROR" << endl;
			cout << "usage: ./motifScan file1.fasta file2.fasta [motif1] [motif2] [windowSize]" << endl;
			cout << "See README file for more" << endl << endl; 
		} 
		else {
			cout << endl;
			cout << "searching..." << endl; 
			location = scanEngine(f1, f2, m1, m2, wSize);
			cout << endl; 
			printf("Reference Seq:  %s\n", f1.c_str()); 
			printf("Query Seq:  %s\n", f2.c_str());
			printf("Motif1:  %s\n", m1.c_str());
			printf("Motif2:  %s\n", m2.c_str());
			printf("Window Size:  %d\n", wSize); 
			cout << endl; 
			cout << "---------------------------------" << endl << endl; 
			if(location == 1){
				printf("No conservation of motifs found\n"); 
				cout << endl; 
			}
			else {
				printf("Conservation of both motifs found at %d on ref seq\n", location); 
				cout << endl; 
			}
		} 
	}
}


