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
        string cname = reads[0].cname;
        ChromDist* forward = ChromDist::AllocateChromDistFitting(reads, PADDING);
        ChromDist* reverse = ChromDist::AllocateChromDistFitting(reads, PADDING);
        
        // First, set up the chromosome distribution
        for(const GFFRow& read : reads){
            // For each read, calculate a normal distribution centered at its start.
            // Go forward and back NORMAL_BOUNDS and add each element of the normal distribution to the chromosome distribution
            for(int i=-NORMAL_BOUNDS; i<+NORMAL_BOUNDS; i++){
                float x = NormalDistribution(i, sigma);
                x *= read.score; // The higher the score, the more confident we are in the particular distribution
                if(read.strand == FORWARD)
                    forward->AddData(i+read.start, x);
                if(read.strand == REVERSE)
                    reverse->AddData(i+read.start, x);
            }
        }
        
        // Call the peaks from the distribution
        vector<GFFRow>* forwardPeaks = CallPeaks(forward, exclusion, FORWARD);
        vector<GFFRow>* reversePeaks = CallPeaks(reverse, exclusion, REVERSE);
        
        // Finally, write the called peaks to the output file
        WritePeaks(forwardPeaks);
        WritePeaks(reversePeaks);
        

        
        
    }
    
    vector<GFFRow>* CallPeaks(ChromDist* dist, int exclusion, Strand strand) const{
        vector<GFFRow>* peaks = new vector<GFFRow>();
        // Loop through the chromosome and look for maxima to call peaks
        for(int i=dist->GetStart()+1; i<dist->GetEnd(); i++){
            if(dist->GetData(i) > dist->GetData(i-1) && dist->GetData(i) > dist->GetData(i+1)){
                GFFRow r;
                r.cname = dist->GetChrom();
                r.source = "genetrack";
                r.start = i - exclusion;
                r.end = i + exclusion;
                r.strand = strand;
                peaks->push_back(r);
            }
        }
        return peaks;
    }
    
    void WritePeaks(vector<GFFRow>* reads){
        for(const GFFRow& read : *reads){
            *OutputFile << read.ToString() << endl;
        }
    }
    
    
private:
    ofstream* OutputFile;
};

