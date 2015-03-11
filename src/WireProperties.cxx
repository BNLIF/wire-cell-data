#include "WireProperties.h"

WireCell::WireProperties::WireProperties(int id, float angle_deg,
					 const std::vector<float>& location_cm)
    : id(id)
    , angle_deg(angle_deg)
    , location_cm(location_cm)
{
}
WireCell::WireProperties::~WireProperties()
{
}

