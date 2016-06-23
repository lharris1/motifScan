/*
 motifScan: This program takes two files in FASTA format and searches for clusters of conserved binding
            sites, as specified by user. 

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

string fileToString(ifstream &inFile) {
	stringstream buffer;
  	buffer << inFile.rdbuf();
  	string infile = buffer.str();
  	return infile;
}

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize) {
	
	ifstream inFile1, inFile2; 

	inFile1.open(file1.c_str());
	if(!inFile1.is_open()) {
		cout << "Unable to open refSeq file." << endl;
		return 0; 
	}

	string refSeq = fileToString(inFile1);
 	size_t mot1_found = refSeq.find(mot1); 
	if(mot1_found!=string::npos){
		cout << "motif1 found on refSeq at: " << mot1_found << endl << endl; 
		size_t mot2_found = refSeq.find(mot2);
		if(mot2_found!=string::npos){              //BOTH MOTIFS FOUND ON REFSEQ!!!
			cout << "motif2 found on refSeq at: " << mot2_found << endl; 
			inFile2.open(file2.c_str());
			if(!inFile2.is_open()) {
				cout << "Unable to open querySeq file." << endl; 
				return 0; 
			}
			string querySeq = fileToString(inFile2);
			size_t mot1_found1 = querySeq.find(mot1); 
			if(mot1_found1!=string::npos){
				cout << "motif1 found on querySeq at: " << mot1_found1 << endl;  
				size_t mot2_found1 = querySeq.find(mot2);
				if(mot2_found1!=string::npos){
					cout << "motif2 found on querySeq at: " << mot2_found1 << endl; 
				}
			}
		}
	}



 	//printf("%s\n", refSeq.c_str()); 
 	inFile1.close(); 
 	inFile2.close(); 
	return 1; //default, if no motif conservation found; else, return seq location where conserved motifs found (on ref seq)
			  //what to do if return 0?? 
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


