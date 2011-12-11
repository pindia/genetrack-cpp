// ChromProcessor.h: main chromosome processing class
// Written by Pindi Albert

#pragma once

#include <fstream>
#include "GFFRow.h"
#include "ChromDist.h"
#include "GFFSearcher.h"

#define PADDING 100
#define NORMAL_BOUNDS 50

class ChromProcessor{
public:
    ChromProcessor(ofstream* outputFile){
        OutputFile = outputFile;
    }
    
    
    void ProcessReads(const vector<GFFRow>& reads, int startIndex, int endIndex, const Options& options){
        string cname = reads[0].cname;
        ChromDist* forward = ChromDist::AllocateChromDistFitting(reads, PADDING);
        ChromDist* reverse = ChromDist::AllocateChromDistFitting(reads, PADDING);
        
        // First, set up the chromosome distribution
        for(const GFFRow& read : reads){
            // For each read, calculate a normal distribution centered at its start.
            // Go forward and back NORMAL_BOUNDS and add each element of the normal distribution to the chromosome distribution
            for(int i=-NORMAL_BOUNDS; i<+NORMAL_BOUNDS; i++){
                float x = NormalDistribution(i, options.sigma);
                x *= read.score; // The higher the score, the more confident we are in the particular distribution
                if(read.strand == FORWARD)
                    forward->AddData(i+read.start, x);
                if(read.strand == REVERSE)
                    reverse->AddData(i+read.start, x);
            }
        }
        
        // Call the peaks from the distribution
        vector<GFFRow>* forwardPeaks = CallPeaks(forward, FORWARD, options);
        vector<GFFRow>* reversePeaks = CallPeaks(reverse, REVERSE, options);
        
        // Process exclusion on the peaks
        
        PerformExclusion(forwardPeaks, options);
        PerformExclusion(reversePeaks, options);
        
        // Finally, write the called peaks to the output file
        WritePeaks(forwardPeaks);
        
        WritePeaks(reversePeaks);
        
        
    }
    
    void PerformExclusion(vector<GFFRow>* peaks, const Options& options){
        GFFSearcher s(peaks);
        for(const GFFRow& peak : *peaks){
            vector<GFFRow>* otherPeaks = s.GetWindow(peak.start, peak.end);
        }
        
    }
    
    vector<GFFRow>* CallPeaks(ChromDist* dist, Strand strand, const Options& options) const{
        vector<GFFRow>* peaks = new vector<GFFRow>();
        // Loop through the chromosome and look for maxima to call peaks
        for(int i=dist->GetStart()+1; i<dist->GetEnd(); i++){
            if( i > 1000)
                break;
            if(dist->GetData(i) > dist->GetData(i-1) && dist->GetData(i) > dist->GetData(i+1)){
                GFFRow r;
                r.cname = dist->GetChrom();
                r.source = "genetrack";
                r.start = i - options.exclusion;
                r.end = i + options.exclusion;
                r.strand = strand;
                r.score = dist->GetData(i);
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

