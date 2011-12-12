// main.cpp: main Genetrack entry point
// Written by Pindi Albert

#include <iostream>
#include <fstream>
#include <list>

#include "GFFReader.h"
#include "ChromProcessor.h"

int main(int argc, char* argv[]){

    
    if(argc < 2){
        cout << "Usage: Genetrack.exe <input_path> <options>" << endl;
        cout << "Valid command-line options: " << endl;
        cout << "-s: Sigma to use when smoothing reads to call peaks. Default 5." << endl;
        cout << "-e: Exclusion zone around each peak that prevents others from being called. Default 10." << endl;
        cout << "-w: Width of called peaks. Default 10." << endl;
        cout << "-F: Absolute read filter; outputs only peaks with larger read count. Default 1." << endl;
        cout << "-k: Size, in base pairs, to chunk each chromosome into when processing. Default 100,000." << endl;
        return 1;
    }
    
    ifstream input(argv[1]);
    GFFReader reader(&input);
    
    ofstream output("output.gff");
    ChromProcessor processor(&output);
    
    Options o; // Initialize default options
    o.sigma = 5;
    o.exclusion = 10;
    o.width = 10;
    o.filter = 1;
    o.chunkSize = 100000;
    
    for(int i=2; (i+1)<argc; i += 2){ // Loop over possible options switches
        string option = argv[i];
        int value = atoi(argv[i+1]);
        if(option == "-s")
            o.sigma = value;
        else if(option == "-e")
            o.exclusion = value;
        else if(option == "-w")
            o.width = value;
        else if(option == "-f")
            o.filter = value;
        else if(option == "-k")
            o.chunkSize = value;
        else{
            cout << "Unrecognized option " << option;
            return 1;
        }
    }
    
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