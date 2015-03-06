#ifndef WireCellProperties_h
#define WireCellProperties_h

#include <vector>

namespace WireCell {

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
