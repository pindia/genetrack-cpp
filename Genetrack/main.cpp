// main.cpp: main Genetrack entry point
// Written by Pindi Albert

#include <iostream>
#include <fstream>
#include <list>

#include "GFFReader.h"
#include "ChromProcessor.h"

int main(){
        
    ifstream input("input.gff");
    GFFReader reader(&input);
    
    vector<GFFRow> chr1 = reader.LoadChromosome();
    
    ofstream output("output.gff");
    ChromProcessor processor(&output);
    
    Options o;
    o.sigma = 5;
    o.exclusion = 20;
    o.width = 20;
    o.filter = 1;
    o.chunkSize = 1000;
    
    int maxIndex = chr1.back().start;
    
    for(int startIndex=0; startIndex<maxIndex; startIndex += o.chunkSize){
        processor.ProcessReads(chr1, startIndex, startIndex + o.chunkSize, o);
    }
    
    
    return 0;
}