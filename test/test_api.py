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
    val = WireCell().Id()
    assert val
    print val
