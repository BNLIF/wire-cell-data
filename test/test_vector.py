#!/usr/bin/env python

import ROOT

WC = ROOT.WireCell

def dump(pt):
    string = ' '.join([str(pt[i]) for i in range(3)])
    print string

def test_inside_box():
    b1 = WC.D3FloatVector(0,0,0)
    b2 = WC.D3FloatVector(1,1,1)
    ray = WC.D3FloatVector(1,1,2)
    center = WC.D3FloatVector(.5,.5,.5)
    hit1 = WC.D3FloatVector(-111,-111,-111)
    hit2 = WC.D3FloatVector(-222,-222,-222)
    hitmask = WC.box_intersection(b1, b2, center, ray, hit1, hit2);
    print hitmask
    for n in [center, ray, hit1, hit2]:
        dump(n)
    assert hitmask

def test_outside_hitting_box():
    b1 = WC.D3FloatVector(0,0,0)
    b2 = WC.D3FloatVector(1,1,1)
    ray = WC.D3FloatVector(1,1,2)
    center = WC.D3FloatVector(-.5,-.5,-.5)
    hit1 = WC.D3FloatVector(-111,-111,-111)
    hit2 = WC.D3FloatVector(-222,-222,-222)
    hitmask = WC.box_intersection(b1, b2, center, ray, hit1, hit2);
    print hitmask
    for n in [center, ray, hit1, hit2]:
        dump(n)
    assert hitmask

def test_outside_missing_box():
    b1 = WC.D3FloatVector(0,0,0)
    b2 = WC.D3FloatVector(1,1,1)
    ray = WC.D3FloatVector(1,1,2)
    center = WC.D3FloatVector(-2,-2,-.5)
    hit1 = WC.D3FloatVector(-111,-111,-111)
    hit2 = WC.D3FloatVector(-222,-222,-222)
    hitmask = WC.box_intersection(b1, b2, center, ray, hit1, hit2);
    print hitmask
    for n in [center, ray, hit1, hit2]:
        dump(n)
    assert  hitmask == 0

test_inside_box()
test_outside_hitting_box()
test_outside_missing_box()

