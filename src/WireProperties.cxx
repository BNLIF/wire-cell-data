#include "WireCell/WireProperties.h"

WireCell::WireProperties::WireProperties(int id, float angle,
					 const std::vector<float>& location)
    : id(id)
    , angle(angle)
    , location(location)
{
}
WireCell::WireProperties::~WireProperties()
{
}

