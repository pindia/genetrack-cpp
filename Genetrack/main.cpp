// main.cpp: main Genetrack entry point
// Written by Pindi Albert

#include <iostream>
#include <fstream>

#include "GFFReader.h"
#include "ChromProcessor.h"

int main(){
    
    
    ifstream input("input.gff");
    GFFReader reader(&input);
    
    vector<GFFRow> chr1 = reader.LoadChromosome();
    
    ofstream output("output.gff");
    ChromProcessor processor(&output);
    
    Options o;
    o.sigma = 10;
    o.exclusion = 10;
    
    processor.ProcessReads(chr1, 0, 0, o);
    
    return 0;
}