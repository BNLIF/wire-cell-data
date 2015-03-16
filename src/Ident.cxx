#include "WireCellData/Ident.h"

WireCellData::Ident::Ident(int cellid, const std::vector<int>& wireids)
    : cell(cellid)
    , wire(wireids)
{
}
WireCellData::Ident::~Ident()
{
}

