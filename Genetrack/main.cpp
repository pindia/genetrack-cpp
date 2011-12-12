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
    
        ofstream output("output.gff");
    ChromProcessor processor(&output);
    
    Options o;
    o.sigma = 5;
    o.exclusion = 20;
    o.width = 20;
    o.filter = 1;
    o.chunkSize = 1000;
    
    do{ // Main loop to process each chromosome in the file

        vector<GFFRow> chr = reader.LoadChromosome();
        
        int maxIndex = chr.back().start;
        
        for(int startIndex=0; startIndex<maxIndex; startIndex += o.chunkSize){
            // Secondary loop to process each chunk in the chromosome
            processor.ProcessReads(chr, startIndex, startIndex + o.chunkSize, o);
        }
        
    } while(!reader.IsDone());
    
    return 0;
}