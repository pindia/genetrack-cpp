// ChromDist.h: Chromosome distribution header file
// Written by Pindi Albert

#pragma once

#include "GFFRow.h"
#include <vector>

// A chromosome distribution stores a float values for each base pair in a chromosome. A single ChromDist object
// stores just a section of a chromosome so it doesn't have to all be held in memory at the same time.
class ChromDist{
public:
    
    // Returns a ChromDist allocated to fit all of the reads given
    // Reads must be in sorted order!
    static ChromDist* AllocateChromDistFitting(const vector<GFFRow>& reads, int padding){
        string chrom = reads[0].cname;
        int minIndex = reads[0].start;
        int maxIndex = reads[reads.size()-1].end;
        return new ChromDist(chrom, minIndex-padding, maxIndex-minIndex+padding); 
    }
    
    ChromDist(string chrom, int start, int length){
        Chrom = chrom;
        Start = start;
        Length = length;
        Data = new float[length];
    }
    
    ~ChromDist(){
        delete [] Data;
    }
    
    int GetStart(){
        return Start;
    }
    
    int GetEnd(){
        return Start + Length - 1;
    }
    
    // GetData, SetData, and AddData methods take *absolute chromosomal coordinates*
    
    int GetData(int index){
        return Data[index - Start];
    }
    
    void SetData(int index, float value){
        Data[index - Start] = value;
    }
    
    void AddData(int index, float value){
        SetData(index, GetData(index) + value);
    }
    
    
    
private:
    string Chrom; // Name of the chromosome
    int Start; // Start index of the distribution on the chromosome
    int Length; // Length of the distribution
    float* Data; // Actual distribution data
    
    
};