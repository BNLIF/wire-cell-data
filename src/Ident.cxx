#include "WireCellData/Ident.h"

using namespace WireCell;

Ident::Ident(int cellid, const std::vector<int>& wireids)
    : cell(cellid)
    , wire(wireids)
{
}

Ident::~Ident()
{
}

