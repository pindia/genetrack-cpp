// GFFSearcher.h: wraps a std::vector<GFFRow>, providing fast binary search within it
//

#include "Common.h"

#include <algorithm>

#pragma once

class GFFSearcher{
public:
    GFFSearcher(vector<GFFRow>* data){
        Data = data;
        Keys = new vector<int>();
        for(const GFFRow& row : *data){
            Keys->push_back(row.start); // Build the keys table for binary search
        }
    }
    
    // Returns a sub-vector with all of the rows with start coordinates within the given boundaries
    vector<GFFRow>* GetWindow(int start, int end){
        int startIndex = BisectLeft(start);
        int endIndex = BisectRight(end);
        vector<GFFRow>* subvector = new vector<GFFRow>();
        for(int i=startIndex; i<endIndex; i++){
            subvector->push_back((*Data)[i]);
        }
        return subvector;
    }
    
    
    vector<GFFRow>* GetData() const{
        return Data;
    }
    
private:
    
    // Returns the index of the first contained GFFRow with start coordinate greater than the given value
    int BisectLeft(int value) const{
        vector<int>::iterator it;
        it = lower_bound(Keys->begin(), Keys->end(), value);
        return (int)(it - Keys->begin());
    }
    
    // Returns the index of the last contained GFFRow with start coordinate less then the given value
    int BisectRight(int value) const{
        vector<int>::iterator it;
        it = upper_bound(Keys->begin(), Keys->end(), value);
        return (int)(it - Keys->begin());        
    }
    
    vector<GFFRow>* Data;
    vector<int>* Keys;
    
};