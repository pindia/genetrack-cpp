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
            if(rows.size() > 0 && rows.back().start == CurrentRow.start && rows.back().strand == CurrentRow.strand){
                rows.back().score += CurrentRow.score; // If the current row is at the same position/strand as the last, merge them
            } else{
                rows.push_back(CurrentRow); // Otherwise, push the read row as a new row
            }
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