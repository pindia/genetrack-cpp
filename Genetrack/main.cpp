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
    
    cout << chr1.size() << endl;
    
    return 0;
}