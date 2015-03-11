#ifndef CellProperties_h
#define CellProperties_h

#include <vector>

namespace WireCell {

class CellProperties {
public:
    CellProperties(int id = 0, float area_cm2=0.0, 
	       const std::vector<float>& center_cm = std::vector<float>());
    ~CellProperties();

    int id;
    float area_cm2;
    std::vector<float> center_cm;
};
}
#endif
