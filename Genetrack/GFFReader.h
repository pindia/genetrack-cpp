// GFFReader.h: GFF3 format file reader 
// Written by Pindi Albert

#pragma once

#include "Common.h"
#include <iostream>
#include <fstream>
#include "GFFRow.h"

class GFFReader{
    
public:
	GFFReader(ifstream* inputFile){
		InputFile = inputFile;
		done = false;
		NextValidRow();
	}
    
	vector<GFFRow> LoadChromosome(){
		vector<GFFRow> rows;
		// Load a single chromosome worth of data from the input file and returns as Vector of GFFRows
		string currentChromosome = CurrentRow.cname;
		while(CurrentRow.cname == currentChromosome){
			rows.push_back(CurrentRow);
			if(!NextValidRow()){
				done = true;
				break;
			}
		}
		return rows;
	}
    
private:
    
	bool NextValidRow(){
		string str;
		do{
			getline(*InputFile, str);
			if((*InputFile).eof())
				return false;
		} while (!GFFRow::IsValidRow(str));
		CurrentRow = GFFRow::ParseRow(str);
		return true;
	}
    
	bool done;
	GFFRow CurrentRow;
	ifstream* InputFile;
    
};