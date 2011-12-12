// GFFSearcher.h: wraps a std::vector<GFFRow>, providing fast binary search within it
//

#pragma once

#include "Common.h"
#include "GFFRow.h"

#include <algorithm>


class GFFSearcher{
public:
    GFFSearcher(vector<GFFRow>* data);
    
    // Returns a sub-vector with all of the rows with start coordinates within the given boundaries
    vector<GFFRow*>* GetWindow(int start, int end);
    
    
private:
    
    // Returns the index of the first contained GFFRow with start coordinate greater than the given value
    int BisectLeft(int value) const;
    
    // Returns the index of the last contained GFFRow with start coordinate less then the given value
    int BisectRight(int value) const;      
    
    vector<GFFRow*>* Data;
    vector<int>* Keys;
    
};