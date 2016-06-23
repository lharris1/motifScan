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

int finderFunc(string seq, string motif, int myStart, int myStop) {
	//cout << myStart << endl; 
	//cout << myStop << endl;
	int strLen = myStop - myStart; 
	string substr = seq.substr(myStart, strLen);  
	size_t found = substr.find(motif);
	if(found!=string::npos) {
		return found+myStart; 
	}
	else {
		return 0; 
	}
}

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize) {
	
	int start, stop; 
	ifstream inFile1, inFile2; 

	inFile1.open(file1.c_str());
	if(!inFile1.is_open()) {
		cout << "Unable to open refSeq file." << endl;
		return 0; 
	}

	string refSeq = fileToString(inFile1);
	start = 0; 
	stop = refSeq.size(); 
 	int mot1_found = finderFunc(refSeq, mot1, start, stop); 

	if(mot1_found == 0){
		cout << "motif1 not found on refSeq " << endl; 
		inFile1.close(); 
		return 1; 
	}
	else {
		cout << "motif1 found on refSeq at: " << mot1_found << endl; 
		start = mot1_found - winSize;
		stop = mot1_found + winSize;  
		int mot2_found = finderFunc(refSeq, mot2, start, stop); 
		if(mot2_found == 0){
			cout << "motif2 not found on refSeq within winSize of " << winSize << endl; 
			inFile1.close(); 
			return 1; 
		}
		else {
			cout << "motif2 found on refSeq at: " << mot2_found << endl; 
			inFile2.open(file2.c_str());
			if(!inFile2.is_open()) {
				cout << "Unable to open querySeq file." << endl; 
				return 0; 
			}
			string querySeq = fileToString(inFile2);
			int mot1_found1 = finderFunc(querySeq, mot1, start, stop);
			if(mot1_found1 == 0){
				cout << "motif1 not found on querySeq" << endl; 
				inFile1.close();
				inFile2.close(); 
				return 1; 
			}
			else {
				cout << "motif1 found on querySeq at: " << mot1_found1 << endl; 
				int mot2_found1 = finderFunc(querySeq, mot2, start, stop); 
				if(mot2_found1 == 0){
					cout << "motif2 not found on querySeq" << endl; 
					inFile1.close();
					inFile2.close(); 
					return 1; 
				}
				else {
					cout << "motif2 found on querySeq at: " << mot2_found1 << endl; 
					inFile1.close(); 
					inFile2.close(); 
					return 2; 
				}
			}
		}
	}
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
			printf("Reference Seq:  %s\n", f1.c_str()); 
			printf("Query Seq:  %s\n", f2.c_str());
			printf("Motif1:  %s\n", m1.c_str());
			printf("Motif2:  %s\n", m2.c_str());
			printf("Window Size:  %d\n", wSize); 
			cout << endl; 
			cout << "---------------------------------" << endl << endl; 
			cout << "searching..." << endl << endl; 
			location = scanEngine(f1, f2, m1, m2, wSize);
			if(location == 1){
				cout << endl; 
				printf("No conservation of motifs found\n"); 
				cout << endl; 
			}
			else {
				cout << endl << endl; 
				cout << "Conserved sites found!!" << endl << endl << endl; 
			}
		} 
	}
}


