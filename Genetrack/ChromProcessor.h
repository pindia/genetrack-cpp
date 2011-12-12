// ChromProcessor.h: main chromosome processing class
// Written by Pindi Albert

#pragma once

#include <fstream>
#include "GFFRow.h"
#include "ChromDist.h"
#include "GFFSearcher.h"

#define PADDING 100
#define NORMAL_BOUNDS 50

bool CompareByStart (const GFFRow& first, const GFFRow& second);

bool CompareByScore (GFFRow* first, GFFRow* second);


class ChromProcessor{
public:
    
    ChromProcessor(ofstream* outputFile);
    
    void ProcessReads(const vector<GFFRow>& reads, int startIndex, int endIndex, const Options& options);

    
private:
    
    void SortByScore(vector<GFFRow*>& peaks);    
    void SortByStart(vector<GFFRow>& peaks);    
    
    void PerformExclusion(vector<GFFRow>*& peaks, const Options& options);
    vector<GFFRow>* CallPeaks(ChromDist* dist, Strand strand, const Options& options) const;    
    void StripZeroPadding(string& cname);    
    void WritePeaks(vector<GFFRow>* reads);    
    
    ofstream* OutputFile;
};

