#ifndef WireCellId_h
#define WireCellId_h

#include <vector>

namespace WireCell {

class Ident {
public:
    Ident(int cellid = 0, const std::vector<int>& wireids = std::vector<int>());
    ~Ident();

    // A number identifying a cell
    int cell;
    
    // A vector of numbers identifying associated wires
    std::vector<int> wire;
};

}
#endif
