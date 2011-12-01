// main.cpp: main Genetrack entry point
// Written by Pindi Albert

#include <iostream>
#include <fstream>

#include "GFFReader.h"

int main(){
    
    cout << "Hello world!";
    
    char* s;
    size_t size;
    s = getcwd(s, size);
    cout << s;
    
    ifstream input("input.gff");
    GFFReader reader(&input);
    
    vector<GFFRow> chr1 = reader.LoadChromosome();
    
    cout << chr1.size();
    
    system("pause");
}