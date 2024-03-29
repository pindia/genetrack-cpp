// ChromDist.h: Chromosome distribution header file
// Written by Pindi Albert

#pragma once

#include "GFFRow.h"
#include <vector>

// A chromosome distribution stores a float values for each base pair in a chromosome. A single ChromDist object
// stores just a section of a chromosome so it doesn't have to all be held in memory at the same time.
class ChromDist{
public:
    
    ChromDist(string chrom, int start, int length){
        Chrom = chrom;
        Start = start;
        Length = length;
        Data = new float[length];
        for(int i=0; i<length; i++) // Make sure data is zeroed
            Data[i] = 0;
    }
    
    ~ChromDist(){
        delete [] Data;
    }
    
    string GetChrom() const{
        return Chrom;
    }
    
    int GetStart() const{
        return Start;
    }
    
    int GetEnd() const{
        return Start + Length - 1;
    }
    
    // GetData, SetData, and AddData methods take *absolute chromosomal coordinates*
    
    int GetData(int index) const{
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