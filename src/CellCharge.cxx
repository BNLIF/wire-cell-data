#include "WCPData/CellCharge.h"

using namespace WCP;

double Cell::charge(const Cell::Group& group)
{
    double tot = 0;
    size_t ncells = group.size();
    for (size_t ind=0; ind<ncells; ++ind) {
	tot += group[ind].second;
    }
    return tot;
}


Cell::GroupCollection Cell::singlets(const Cell::Group& group)
{
    Cell::GroupCollection ret;
    size_t count = group.size();
    for (size_t ind=0; ind<count; ++ind) {
	Cell::Group g;
	g.push_back(group[ind]);
	ret.push_back(g);
    }
    return ret;
}

double Cell::cross_section(const Cell::Group& group)
{
    double tot = 0;
    size_t ncells = group.size();
    for (size_t ind=0; ind<ncells; ++ind) {
	tot += group[ind].first->cross_section();
    }
    return tot;

}

Point Cell::center_of_charge(const Cell::Group& group)
{
    double tot_charge = 0;
    double x=0, y=0, z=0;
    size_t ncells = group.size();

    if (!ncells) {
	return Point();
    }

    for (size_t ind=0; ind<ncells; ++ind) {
	double q = group[ind].second;
	tot_charge += q;
	Point p = group[ind].first->center();
	x += q*p.x;
	y += q*p.y;
	z += q*p.z;
    }

    if (tot_charge > 0.0) {
	return Point(x/tot_charge, y/tot_charge, z/tot_charge);
    }
    return Point();    
}
