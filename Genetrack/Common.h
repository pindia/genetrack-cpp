// Common.h: Genetrack common functions
// Written by Pindi Albert

#pragma once

#include "Options.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <cmath>
using namespace std;

vector<string> strsplit(const string& s, const string& delim, const bool keep_empty = true) {
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;
    while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            result.push_back(temp);
        }
        if (subend == s.end()) {
            break;
        }
        substart = subend + delim.size();
    }
    return result;
}


int str2int (const string &str) {
    stringstream ss(str);
    int num;
    if((ss >> num).fail())
    { 
        //ERROR 
    }
    return num;
}

float str2float (const string &str) {
    stringstream ss(str);
    float num;
    if((ss >> num).fail())
    { 
        //ERROR 
    }
    return num;
}

int strcount(const string& str, const char c){
	int count = 0;
	for(int i=0; i<str.length(); i++)
		if(str[i] == c)
			count++;
	return count;
}

float NormalDistribution(float x, float sigma){
    return exp( (-x * x) /  (2 * sigma * sigma));
}

