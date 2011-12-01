// GFFRow.h: representation of a single line of a GFF3 file
// Written by Pindi Albert

#pragma once

#include "Common.h"
#include <map>


enum Strand{
	FORWARD,
	REVERSE,
	NONE
};

struct GFFRow{
	string cname;  // Chromosome name
	string source; // Script that produced feature
	string type;   // Type of feature 
	int start;     // Chromosomal coordinates of feature
	int end;       // INCLUSIVE on both ends
	float score;   // Confidence level of feature
	Strand strand; // Strand feature is located on
	string phase;  // Phase of feature
	map<string, string> attrs; // Miscellaneous attributes
	
	static bool IsValidRow(string str){
		return strcount(str, '\t') == 8;
	}
    
	static GFFRow ParseRow(string str){
		vector<string> vec = strsplit(str, "\t");
		GFFRow row;
		row.cname = vec[0];
		row.source = vec[1];
		row.type = vec[2];
		row.start = str2int(vec[3]);
		row.end = str2int(vec[4]);
		row.score = str2float(vec[5]);
		if(vec[6] == "+")
			row.strand = FORWARD;
		else if(vec[6] == "-")
			row.strand = REVERSE;
		else
			row.strand = NONE;
		row.phase = vec[7];
		// Column 8 has the format "attr1=val1;attr2=val2"
		vector<string> splitattrs = strsplit(vec[8], ";"); // Start by splitting on ;
		if(splitattrs.size() > 1){ // Only if ; is actually found
			for(int i=0; i<splitattrs.size(); i++){
				vector<string> comps = strsplit(splitattrs[i], "="); // Split each attr=val pair on the =
				row.attrs[comps[0]] = comps[1];
			}
		}
		return row;
	}
    
};