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

string fileToString(string fileName) {
	string header, seq; 
	ifstream myFile(fileName);
	getline(myFile, header);
	getline(myFile, seq);
	myFile.close();
	return seq; 
}

int finderFunc(string seq, string motif, int myStart, int myStop) {
	int my_start, my_stop; 
	if(myStart < 0){
		my_start = 0; 
	}
	else{
		my_start = myStart; 
	}
	if(myStop > seq.size()){
		my_stop = seq.size();
	}
	else{
		my_stop = myStop; 
	}
	//cout << my_start << endl; //DELETE ME!!!
	//cout << myStop << endl;   //DELETE ME!!!
	int strLen = my_stop - my_start; 
	string substr = seq.substr(my_start, strLen);  
	size_t found = substr.find(motif);
	if(found!=string::npos) {
		return found + my_start; 
	}
	else {
		return 0; 
	}
}

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize) {


	string refSeq = fileToString(file1);
	string querySeq = fileToString(file2);

	int start = 0; 
	int stop = refSeq.size(); 
 	int mot1_found = finderFunc(refSeq, mot1, start, stop); 

	if(mot1_found == 0){
		cout << "motif1 not found on refSeq " << endl; 
		return 1; 
	}
	else {
		cout << "motif1 found on refSeq at: " << mot1_found << endl; 
		start = mot1_found - winSize;
		stop = mot1_found + winSize;  
		int mot2_found = finderFunc(refSeq, mot2, start, stop); 
		if(mot2_found == 0){
			cout << "motif2 not found on refSeq within winSize of " << winSize << endl; 
			return 1; 
		}
		else {
			cout << "motif2 found on refSeq at: " << mot2_found << endl;  
			int mot1_found1 = finderFunc(querySeq, mot1, start, stop);
			if(mot1_found1 == 0){
				cout << "motif1 not found on querySeq" << endl; 
				return 1; 
			}
			else {
				cout << "motif1 found on querySeq at: " << mot1_found1 << endl; 
				int mot2_found1 = finderFunc(querySeq, mot2, start, stop); 
				if(mot2_found1 == 0){
					cout << "motif2 not found on querySeq" << endl; 
					return 1; 
				}
				else {
					cout << "motif2 found on querySeq at: " << mot2_found1 << endl; 
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


