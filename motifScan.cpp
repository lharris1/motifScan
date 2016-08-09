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
#include <algorithm>

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

string getRevComp(string mot1) {
	string new_mot1; 
	for(int i=0; i<mot1.size(); i++) {
		stringstream ss;
		string myLetter; 
		char letter = mot1.at(i);
		ss << letter;
		ss >> myLetter; 

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
	return new_mot1; 
}

vector<string> getWobbleMotifs(string mot) {
	vector<string> bigList, bigList_temp, newMotList, newMotList_temp, toReturn;

	for(int i=0; i<mot.size(); i++) {
		stringstream kk; 
		string myLetter;
		char letter = mot.at(i);
		kk << letter; 
		kk >> myLetter; 

		if(myLetter.compare("W")==0){
			for(int z=0; z<newMotList.size(); z++){
				string newMotW = newMotList.at(z);
				string newMot1 = newMotList.at(z);
				newMotW.append("T");
				newMot1.append("A");
				bigList_temp.push_back(newMotW);
				bigList_temp.push_back(newMot1); 
				newMotList_temp.push_back(newMotW);
				newMotList_temp.push_back(newMot1);
			}
		}
		else if(myLetter.compare("S")==0){
			for(int z=0; z<newMotList.size(); z++){
				string newMotS = newMotList.at(z); 
				string newMot2 = newMotList.at(z); 
				newMotS.append("C");
				newMot2.append("G");
				bigList_temp.push_back(newMotS);
				bigList_temp.push_back(newMot2); 
				newMotList_temp.push_back(newMotS);
				newMotList_temp.push_back(newMot2);
			}
		}
		else if(myLetter.compare("K")==0) {
			for(int z=0; z<newMotList.size(); z++){
				string newMotK = newMotList.at(z); 
				string newMot3 = newMotList.at(z);
				newMotK.append("G");
				newMot3.append("T");
				bigList_temp.push_back(newMotK);
				bigList_temp.push_back(newMot3); 
				newMotList_temp.push_back(newMotK);
				newMotList_temp.push_back(newMot3); 
			}
		}
		else if(myLetter.compare("M")==0) {
			for(int z=0; z<newMotList.size(); z++){
				string newMotM = newMotList.at(z);
				string newMot4 = newMotList.at(z); 
				newMotM.append("A");
				newMot4.append("C");
				bigList_temp.push_back(newMotM);
				bigList_temp.push_back(newMot4); 
				newMotList_temp.push_back(newMotM);
				newMotList_temp.push_back(newMot4);
			}
		}
		else if(myLetter.compare("Y")==0) {
			for(int z=0; z<newMotList.size(); z++){
				string newMotY = newMotList.at(z); 
				string newMot5 = newMotList.at(z);
				newMotY.append("C");
				newMot5.append("T");
				bigList_temp.push_back(newMotY);
				bigList_temp.push_back(newMot5); 
				newMotList_temp.push_back(newMotY);
				newMotList_temp.push_back(newMot5);
			}
		}
		else if(myLetter.compare("R")==0) {
			for(int z=0; z<newMotList.size(); z++){
				string newMotR = newMotList.at(z); 
				string newMot6 = newMotList.at(z);
				newMotR.append("A");
				newMot6.append("G");
				bigList_temp.push_back(newMotR);
				bigList_temp.push_back(newMot6); 
				newMotList_temp.push_back(newMotR);
				newMotList_temp.push_back(newMot6);
			}
		}
		else {
			string newMot0 = myLetter; 
			newMotList_temp.push_back(newMot0); 
			for(int v=0; v<newMotList.size(); v++){
				string newMot = newMotList.at(v);
				newMot.append(myLetter);
				newMotList_temp.push_back(newMot);
				bigList_temp.push_back(newMot); 
			}
		}
		bigList = bigList_temp; 
		bigList_temp.clear(); 
		newMotList = newMotList_temp;
		newMotList_temp.clear(); 
	}
	for(int y=0; y<bigList.size(); y++){
		string item = bigList.at(y);
		if(item.size()==mot.size()){
			toReturn.push_back(item);
		}
	}

	return toReturn; 
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
		return -1; 
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

int scanEngine(string file1, string file2, string mot1, string mot2, int winSize) {

	string refSeq = fileToString(file1);
	string querySeq = fileToString(file2);
	
	int start = 0, memVar = 0, memVar1 = 0; 
	int my_stop = refSeq.size(); 
	int stop = my_stop; 
	int index = mot1.size();
	int index1 = mot2.size(); 
	int clusterCount = 0, myClusterCount = 0; 

	while(stop <= refSeq.size()) {
 		int mot1_found = finderFunc(refSeq, mot1, start, (stop+index1)); 
		if(mot1_found == -1){
			break; 
		}
		else {                                   //mot1 found on ref
			start = mot1_found - winSize;
			stop = mot1_found + winSize;  
			int mot2_found = finderFunc(refSeq, mot2, start, stop); 
			if(mot2_found == -1){
				start = mot1_found+index; 
				stop = my_stop; 
			}
			else {                               //mot2 found on ref
				int mot1_found1 = finderFunc(querySeq, mot1, start, stop);
				if(mot1_found1 == -1){
					start = mot1_found+index; 
					stop = my_stop; 
				}
				else {                           //mot1 found on query...almost there...
					int mot2_found1 = finderFunc(querySeq, mot2, start, stop); 
					if(mot2_found1 == -1){
						start = mot1_found+index; 
						stop = my_stop; 
					}
					else {                    //mot2 found on query....success!!!
						if((mot1_found1 == memVar) && (mot2_found1 == memVar1)){
							break; 
						}  
						else {                //mot1 and mot2 have not been found on query before
							clusterCount += 1; 
							outputClusterFound(mot1_found, mot2_found, mot1_found1, mot2_found1);
							start = mot1_found+index;
							if(stop> my_stop) {
								break; 
							}
							stop = my_stop; 
							memVar = mot1_found1; 
							memVar1 =  mot2_found1; 
						}
					}
				}
			}	
		}
	} 
	myClusterCount = clusterCount;
	clusterCount = 0; 
	return myClusterCount; 
}

void printStartUp(string myFile1, string myFile2, string myMotif1, string myMotif2, int myWinSize, int myRevComp, int myWobble, int myAlt){
		
		if(myWinSize == 0) {
			cout << endl; 
			cout << "ERROR" << endl;
			cout << "usage: ./motifScan file1.fasta file2.fasta [motif1] [motif2] [windowSize] [revComp Y/N] [wobble Y/N] [getAltMotifs Y/N]" << endl;
			cout << "See README file for more" << endl << endl; 
		} 
		else {
			cout << endl;
			printf("Reference Seq:  %s\n", myFile1.c_str()); 
			printf("Query Seq:  %s\n", myFile2.c_str());
			printf("Motif1:  %s\n", myMotif1.c_str());
			printf("Motif2:  %s\n", myMotif2.c_str());
			printf("Window Size:  %d\n", myWinSize); 
			if(myRevComp==0) {
				printf("RevComp:  OFF\n");
			}
			else {
				printf("RevComp:  ON\n");
			}
			if(myWobble==0) {
				printf("WobbleBase:  OFF\n");
			}
			else{
				printf("WobbleBase:  ON\n");
			}
			if(myAlt==0){
				printf("getAltMotifs:  OFF\n");
			}
			else{
				printf("getAltMotifs:  ON\n");
			}
			cout << endl; 
			cout << "----------------------------------------" << endl << endl; 
			cout << "searching..." << endl << endl; 
		}
}

vector<string> getAltMotifs(string myFile1, string myFile2, string myMotif1, int myWinSize){

	string refSeq = fileToString(myFile1);
	string querySeq = fileToString(myFile2);
	
	int start = 0;
	int my_stop = refSeq.size(); 
	int stop = my_stop; 
	int index = myMotif1.size();
	int refStart, refStop, qStart, qStop; 
	vector<string> bigList; 

	while(stop <= refSeq.size()) {
 		int mot1_found = finderFunc(refSeq, myMotif1, start, (stop+index)); 
		if(mot1_found == -1){                     //mot1 not on ref
			break; 
		}
		else {                                   //mot1 found on ref
			start = mot1_found - myWinSize;
			stop = mot1_found + myWinSize;  
			refStart = start; 
			refStop = stop; 

			if(refStop>my_stop){
				refStop = my_stop; 
			}
			int mot1_found1 = finderFunc(querySeq, myMotif1, start, stop); 
			if(mot1_found1 == -1){               //mot1 NOT found on query
				start = mot1_found+index; 
				stop = my_stop; 
				//break;      //not working, should probably leave out
			}
			else {                               //mot1 found on query
				qStart = mot1_found1 - myWinSize; 
				qStop = mot1_found1 + myWinSize; 
				
				for(int n=refStart; n<refStop; n++) {     //driving loop -- splitting and searching
					string altMot = refSeq.substr(n,4); 
					if(finderFunc(querySeq, altMot, qStart, qStop)!=0) {
						if(find(bigList.begin(), bigList.end(), altMot) == bigList.end()) {
							bigList.push_back(altMot);
						}
					}
				}
				if(stop> my_stop) {
					break; 
				}
				start = mot1_found+index; 
				stop = my_stop; 
			}
		}
	}
	return bigList; 
}

int main(int argc, char* argv[]) {

	if(argc != 9) {
		cout << endl; 
		cout << "ERROR" << endl;
		cout << "usage: ./motifScan file1.fasta file2.fasta [motif1] [motif2] [windowSize] [revComp Y/N] [wobble Y/N] [getAltMotifs Y/N]" << endl;
		cout << "See README file for more" << endl << endl; 
	}
	else {
		string f1 = string(argv[1]);
		string f2 = string(argv[2]);
		string m1 = string(argv[3]);
		string m2 = string(argv[4]); 
		int wSize = atoi(argv[5]);
		int revComp = atoi(argv[6]);
		int wobble = atoi(argv[7]); 
		int altMots = atoi(argv[8]);
		int myCount = 0; 

		printStartUp(f1, f2, m1, m2, wSize, revComp, wobble, altMots);

		if(wobble==1){
			vector<string> wobbleMots1 = getWobbleMotifs(m1); 
			vector<string> wobbleMots2 = getWobbleMotifs(m2);
			for(int i=0; i<wobbleMots1.size(); i++){
				string newMot1 = wobbleMots1.at(i); 
				for(int k=0; k<wobbleMots2.size(); k++){
					string newMot2 = wobbleMots2.at(k);
					int temp = scanEngine(f1,f2,newMot1,newMot2, wSize);
					myCount = myCount+temp; 
				}
			}
			if(revComp==1){
				vector<string> revCompList1, revCompList2;
				for(int j=0; j<wobbleMots1.size(); j++) {
					string revComp1 = getRevComp(wobbleMots1.at(j));
					revCompList1.push_back(revComp1);
				}
				for(int l=0; l<wobbleMots2.size(); l++) {
					string revComp2 = getRevComp(wobbleMots2.at(l));
					revCompList2.push_back(revComp2);
				}
				for(int m=0; m<revCompList1.size(); m++){
					string myMot1 = revCompList1.at(m);
					for(int n=0; n<revCompList2.size(); n++){
						string myMot2 = revCompList2.at(n);
						int temp3 = scanEngine(f1,f2,myMot1,myMot2, wSize);
						myCount = myCount+temp3;
					}
				}
			}
		}
		if(revComp==1){
			string revComp1 = getRevComp(m1); 
			string revComp2 = getRevComp(m2);
			int temp4 = scanEngine(f1, f2, revComp1, revComp2, wSize);
			myCount = myCount + temp4; 
		}
		int temp5 = scanEngine(f1, f2, m1, m2, wSize);
		myCount = myCount + temp5; 
		if(myCount == 0){
			cout << endl << "No conserved binding site clusters found" << endl << endl; 
		}
		else {
			cout << endl << myCount << " conserved binding sites found!!!" << endl << endl; 
		}
		if(altMots==1) {
			vector<string> tempList = getAltMotifs(f1,f2,m1,wSize);
			vector<string> finalList = tempList; 
			vector<string> tempList1;
			vector<string> tempList2;

			if(wobble==1){
				vector<string> my_wobbleMots = getWobbleMotifs(m1); 
				for(int i=0; i<my_wobbleMots.size(); i++){
					string newMot1 = my_wobbleMots.at(i); 
					tempList1 = getAltMotifs(f1,f2,newMot1,wSize);
					for(int k=0; k<tempList1.size(); k++){
						string myMot = tempList1.at(k);
						if(find(finalList.begin(), finalList.end(), myMot) == finalList.end()) {
							finalList.push_back(myMot); 
						}
					}
				}
			}
			if(revComp==1){
				string my_revComp = getRevComp(m1);
				tempList2 = getAltMotifs(f1,f2,my_revComp,wSize);
			}
			
			for(int z=0; z<tempList2.size(); z++){
				string myMot1 = tempList2.at(z);
				if(find(finalList.begin(), finalList.end(), myMot1) == finalList.end()) {
					finalList.push_back(myMot1); 
				}
			}
			cout << "----------------------------------------" << endl; 
			cout << "Potential co-motifs to search for: " << endl << endl; 
			cout << finalList.size() << endl; 
			for(int p=0; p<finalList.size(); p++){
				cout << finalList.at(p) << endl; 
			}
			cout << endl; 

			
		}
	} 
}

















