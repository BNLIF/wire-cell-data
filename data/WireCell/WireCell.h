#ifndef WireCell_h
#define WireCell_h

#include <vector>

namespace WireCell {

class Id {
public:
    Id(int cellid = 0, const std::vector<int>& wireids = std::vector<int>());
    ~Id();

    // A number identifying a cell
    int cell;
    
    // A vector of numbers identifying associated wires
    std::vector<int> wire;
};

class Properties {
public:
    Properties(int id = 0, float area_cm2=0.0, 
	       const std::vector<float>& center_cm = std::vector<float>());
    ~Properties();

    int cell;
    float area_cm2;
    std::vector<float> center_cm;
};
}
#endif
