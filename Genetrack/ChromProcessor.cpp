// ChromProcessor.h: main chromosome processing implementation
// Written by Pindi Albert

#include "ChromProcessor.h"


bool CompareByStart (const GFFRow& first, const GFFRow& second)
{
    if (first.start<second.start) return true;
    else return false;
}

bool CompareByScore (GFFRow* first, GFFRow* second)
{
    if (first->score>second->score) return true;
    else return false;
}


ChromProcessor::ChromProcessor(ofstream* outputFile){
        OutputFile = outputFile;
}
    
    
void ChromProcessor::ProcessReads(const vector<GFFRow>& reads, int startIndex, int endIndex, const Options& options){
    
    string cname = reads[0].cname;
    
    cout << "Processing " << cname << " " << startIndex << "-" << endIndex << endl;
    
    
    ChromDist* forward = new ChromDist(cname, startIndex - PADDING, endIndex + PADDING);
    ChromDist* reverse = new ChromDist(cname, startIndex - PADDING, endIndex + PADDING);
    
    // First, set up the chromosome distribution
    for(const GFFRow& read : reads){
        
        if(read.start < startIndex || read.end > endIndex)
            continue; // If read is not in process bounds, ignore it
        
        if(read.score <= options.filter)
            continue; // If read does not pass filter, do not process it
        
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

void ChromProcessor::SortByScore(vector<GFFRow*>& peaks){
    sort(peaks.begin(), peaks.end(), CompareByScore);
    
}

void ChromProcessor::SortByStart(vector<GFFRow>& peaks){
    sort(peaks.begin(), peaks.end(), CompareByStart);
    
}


void ChromProcessor::PerformExclusion(vector<GFFRow>*& peaks, const Options& options){
    
    vector<GFFRow>* newPeaks = new vector<GFFRow>();
    
    GFFSearcher s(peaks); // Construct GFFSearcher with the list of peaks (currently sorted by index)
    
    // Copy the peaks into a new vector to sort by score
    // This vector stores pointers to the original peaks so that they can be mutated
    vector<GFFRow*> peaksByScore;
    for(GFFRow& row : *peaks){
        peaksByScore.push_back(&row);
    }
    SortByScore(peaksByScore);
    
    
    // Loop over peaks by score, so that peaks with higher score have priority over peaks with lower score
    for(GFFRow* peak : peaksByScore){
        if(peak->score == 0)
            continue; // This peak has already been excluded; don't let it exclude anything itself
        newPeaks->push_back(*peak); // Add the non-excluded peaks
        // Exclude other peaks
        vector<GFFRow*>* otherPeaks = s.GetWindow(peak->start - options.exclusion, peak->start + options.exclusion);
        for(GFFRow* otherPeak : *otherPeaks){
            if(peak->start != otherPeak->start){
                otherPeak->score = 0; // Zero out the score for excluded peaks so we know they're excluded
            }
        }
        
    }
    
    SortByStart(*newPeaks); // New peaks should be once again sorted by index
    
    peaks = newPeaks;
    
}

vector<GFFRow>* ChromProcessor::CallPeaks(ChromDist* dist, Strand strand, const Options& options) const{
    vector<GFFRow>* peaks = new vector<GFFRow>();
    // Loop through the chromosome and look for maxima to call peaks
    for(int i=dist->GetStart()+1; i<dist->GetEnd(); i++){
        if(dist->GetData(i) > dist->GetData(i-1) && dist->GetData(i) > dist->GetData(i+1)){
            GFFRow r;
            r.cname = dist->GetChrom();
            r.source = "genetrack";
            r.start = i - options.width; // Use peak width from options
            r.end = i + options.width;
            r.strand = strand;
            r.score = dist->GetData(i);
            peaks->push_back(r);
        }
    }
    return peaks;
}

void ChromProcessor::StripZeroPadding(string& cname){
    size_t pos;
    pos = cname.find('0');
    if(pos != string::npos){
        cname.erase(pos, 1); // Erase the single 0 character if found
    }
}

void ChromProcessor::WritePeaks(vector<GFFRow>* reads){
    for(GFFRow& read : *reads){
        StripZeroPadding(read.cname); // Translate chr01 to chr1
        *OutputFile << read.ToString() << endl;
    }
}
    

