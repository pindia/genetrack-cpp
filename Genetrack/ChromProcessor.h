// ChromProcessor.h: main chromosome processing class
// Written by Pindi Albert

#pragma once

#include <fstream>
#include "GFFRow.h"
#include "ChromDist.h"

#define PADDING 100
#define NORMAL_BOUNDS 50

class ChromProcessor{
public:
    ChromProcessor(ofstream* outputFile){
        OutputFile = outputFile;
    }
    
    void ProcessReads(const vector<GFFRow>& reads, int startIndex, int endIndex, int sigma, int exclusion){
        ChromDist* dist = ChromDist::AllocateChromDistFitting(reads, PADDING);
        
        // First, set up the chromosome distribution
        for(const GFFRow& read : reads){
            // For each read, calculate a normal distribution centered at its start.
            // Go forward and back NORMAL_BOUNDS and add each element of the normal distribution to the chromosome distribution
            for(int i=-NORMAL_BOUNDS; i<+NORMAL_BOUNDS; i++){
                float x = NormalDistribution(i, sigma);
                x *= read.score; // The higher the score, the more confident we are in the particular distribution
                dist->AddData(i+read.start, x);
            }
        }
        
        // Now, loop through the chromosome and look for maxima to call peaks
        for(int i=dist->GetStart(); i<=dist->GetEnd(); i++){
            
        }
        
        
    }
    
    
private:
    ofstream* OutputFile;
};

