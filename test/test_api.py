#!/usr/bin/env python

def WCP():
    import ROOT
    return ROOT.WCP

def test_make_wire():
    val = WCP().GeomWire()
    assert val
    print val
def test_make_cell():
    val = WCP().GeomCell()
    assert val
    print val

def test_make_point():
    Point = WCP().Point
    PointVector = WCP().PointVector
    assert Point()
    p = Point(1.0,2.0,3.0)
    assert p
    assert p[0] == 1.0
    assert p[1] == 2.0
    assert p[2] == 3.0
    p2 = Point(2,4,6)
    assert p2[0] == 2
    assert p2[1] == 4
    assert p2[2] == 6
    pv = PointVector()
    assert pv
    pv.push_back(p)
    pv.push_back(p2)
    assert pv.size() == 2
    pp = pv[0]
    assert pp
    assert pp[0] == 1.0
    assert pp[1] == 2.0
    assert pp[2] == 3.0
    pp2 = pv[1]
    assert pp2
    assert pp2[0] == 2
    assert pp2[1] == 4
    assert pp2[2] == 6
    
    pt = tuple([(x[0],x[1],x[2]) for x in pv])
    assert pt == ((1.0,2.0,3.0),(2,4,6))


def test_units():
    import ROOT.units
    print ROOT.units
    print "mm are", str(ROOT.units.mm)
    assert ROOT.units.mm == 1.0
    assert ROOT.units.radian == 1.0
    assert ROOT.units.ns == 1.0
    assert ROOT.units.kelvin == 1.0
    assert ROOT.units.mole == 1.0
    assert ROOT.units.MeV == 1.0
