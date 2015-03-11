#ifndef WireProperties_h
#define WireProperties_h

#include <vector>

namespace WireCell {

class WireProperties {
public:
    WireProperties(int id = 0, float angle_deg = 0.0,
		   const std::vector<float>& location_cm = std::vector<float>());
    ~WireProperties();

    int id;
    float angle_deg;
    std::vector<float> location_cm;
};
}
#endif
