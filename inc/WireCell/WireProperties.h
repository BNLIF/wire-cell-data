#ifndef WireProperties_h
#define WireProperties_h

#include "WireCell/Units.h"
#include <vector>

namespace WireCell {

class WireProperties {
public:
    WireProperties(int id = 0, float angle = 0.0 * units::degree,
		   const std::vector<float>& location = std::vector<float>());
    ~WireProperties();

    int id;
    float angle;
    std::vector<float> location;
};
}
#endif
