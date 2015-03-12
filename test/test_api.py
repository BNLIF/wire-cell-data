#!/usr/bin/env python

def WireCell():
    import ROOT
    return ROOT.WireCell

def test_make_id():
    val = WireCell().Id()
    assert val
    print val
def test_make_wire():
    val = WireCell().Wire()
    assert val
    print val
def test_make_cell():
    val = WireCell().Cell()
    assert val
    print val

def test_make_point():
    Point = WireCell().Point
    PointVector = WireCell().PointVector
    assert Point()
    p = Point(1.0,2.0)
    assert p
    assert p.first == 1.0
    assert p.second == 2.0
    p2 = Point(2,4)
    assert p2.first == 2.0
    assert p2.second == 4.0
    pv = PointVector()
    assert pv
    pv.push_back(p)
    pv.push_back(p2)
    assert pv.size() == 2
    pp = pv[0]
    assert pp
    assert pp.first == 1.0
    assert pp.second == 2.0
    pp2 = pv[1]
    assert pp2
    assert pp2.first == 2
    assert pp2.second == 4
    
    pt = tuple([(x.first,x.second) for x in pv])
    assert pt == ((1.0,2.0),(2.0,4.0))


def test_units():
    import ROOT
    assert ROOT.units.mm == 1.0
    assert ROOT.units.radian == 1.0
    assert ROOT.units.ns == 1.0
    assert ROOT.units.kelvin == 1.0
    assert ROOT.units.mole == 1.0
    assert ROOT.units.MeV == 1.0
