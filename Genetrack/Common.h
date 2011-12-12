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

vector<string> strsplit(const string& s, const string& delim, const bool keep_empty = true);

int str2int (const string &str);
float str2float (const string &str);

int strcount(const string& str, const char c);

float NormalDistribution(float x, float sigma);
