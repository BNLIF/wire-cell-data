#ifndef CellProperties_h
#define CellProperties_h

#include "WireCell/Units.h"
#include <vector>

namespace WireCell {

class CellProperties {
public:
    CellProperties(int id = 0, float area=0.0 * units::radian, 
	       const std::vector<float>& center = std::vector<float>());
    ~CellProperties();

    int id;
    float area;
    std::vector<float> center;
};
}
#endif
