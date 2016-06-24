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
#include <vector>
#include <stdio.h> 

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//////////////////////////MOTIFSCAN IMPLEMENTATION/////////////////////////////
///////////////////////////////////////////////////////////////////////////////

string fileToString(string fileName) {
	string header, line, fullSeq; 
	ifstream myFile(fileName);
	getline(myFile, header);
	while(getline(myFile, line)){
		fullSeq.append(line);
	}
	myFile.close();
	return fullSeq; 
}

vector<string>* getRevComps(string mot1, string mot2) {
	vector<string>* bigList = new vector<string>; 
	string new_mot1, new_mot2; 
	//cout << mot1 << endl;
	//cout << mot2 << endl; 

	for(int i=0; i<mot1.size(); i++) {
		stringstream ss;
		string myLetter; 
		char letter = mot1.at(i);
		cout << "letter: " << letter << endl;
		ss << letter;
		ss >> myLetter; 
		cout << "myLetter: " << myLetter << endl; 
		//string myLetter = string(letter);
		if(myLetter.compare("A")==0){
			new_mot1.append("T");
		}
		else if(myLetter.compare("T")==0){
			new_mot1.append("A");
		}
		else if(myLetter.compare("G")==0){
			new_mot1.append("C");
		}
		else{
			new_mot1.append("G");
		}
		myLetter = " ";
	}
	cout << new_mot1 << endl; 
	bigList->push_back(new_mot1);

	for(int k=0; k<mot2.size(); k++) {
		stringstream kk; 
		string myLetter1; 
		char letter1 = mot2.at(k);
		kk << letter1;
		kk >> myLetter1;
		if(myLetter1.compare("A")==0){
			new_mot2.append("T");
		}
		else if(myLetter1.compare("T")==0){
			new_mot2.append("A");
		}
		else if(myLetter1.compare("G")==0){
			new_mot2.append("C");
		}
		else{
			new_mot2.append("G");
		}
		myLetter1 = " ";
	}
	cout << new_mot2 << endl;
	bigList->push_back(new_mot2);	
	return bigList; 
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

void outputClusterFound(int m1Ref, int m2Ref, int m1Query, int m2Query) {
	cout << "----------------------------------------" << endl; 
	cout << "       CONSERVED CLUSTER FOUND!!        " << endl; 
	cout << "motif1 found on refSeq at: " << m1Ref<< endl;
	cout << "motif2 found on refSeq at: " << m2Ref << endl;
	cout << "motif1 found on querySeq at: " << m1Query << endl; 
	cout << "motif2 found on querySeq at: " << m2Query << endl; 
	cout << "----------------------------------------" << endl; 
	return; 
}

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize, int rev) {

	string refSeq = fileToString(file1);
	string querySeq = fileToString(file2);

	if(rev==1){
		vector<string>* revComps = getRevComps(mot1, mot2); 
	}
	
	int start = 0; 
	int my_stop = refSeq.size(); 
	int stop = my_stop; 
	int index = mot1.size();
	int index1 = mot2.size(); 

	while(stop <= refSeq.size()) {
 		int mot1_found = finderFunc(refSeq, mot1, start, (stop+index1)); 
		if(mot1_found == 0){
			return 1; 
		}
		else {
			start = mot1_found - winSize;
			stop = mot1_found + winSize;  
			int mot2_found = finderFunc(refSeq, mot2, start, stop); 
			if(mot2_found == 0){
				start = mot1_found+index; 
				stop = my_stop; 
			}
			else {
				int mot1_found1 = finderFunc(querySeq, mot1, start, stop);
				if(mot1_found1 == 0){
					start = mot1_found+index; 
					stop = my_stop; 
				}
				else {
					int mot2_found1 = finderFunc(querySeq, mot2, start, stop); 
					if(mot2_found1 == 0){
						start = mot1_found+index; 
						stop = my_stop; 
					}
					else {
						outputClusterFound(mot1_found, mot2_found, mot1_found1, mot2_found1);
						start = mot1_found+index;
						if(stop> my_stop) {
							return 1; 
						}
						stop = my_stop; 

					}
				}
			}	
		}
	} 
	return 1; 
}

int main(int argc, char* argv[]) {

	int location; 

	if(argc != 7) {
		cout << endl; 
		cout << "ERROR" << endl;
		cout << "usage: ./motifScan file1.fasta file2.fasta [motif1] [motif2] [windowSize] [revComp]" << endl;
		cout << "See README file for more" << endl << endl; 
	}
	else {
		string f1 = string(argv[1]);
		string f2 = string(argv[2]);
		string m1 = string(argv[3]);
		string m2 = string(argv[4]); 
		int revComp = atoi(argv[6]);
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
			if(revComp==0) {
				printf("RevComp:    OFF");
			}
			else {
				printf("RevComp:    ON");
			}
			cout << endl; 
			cout << "----------------------------------------" << endl << endl; 
			cout << "searching..." << endl << endl; 
			location = scanEngine(f1, f2, m1, m2, wSize, revComp);
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


