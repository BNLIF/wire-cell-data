#include "WireCell/Ident.h"

WireCell::Ident::Ident(int cellid, const std::vector<int>& wireids)
    : cell(cellid)
    , wire(wireids)
{
}
WireCell::Ident::~Ident()
{
}

