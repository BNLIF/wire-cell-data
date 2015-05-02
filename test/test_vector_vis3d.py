#!/usr/bin/env python

import ROOT
from ROOT.WireCell import D3FloatVector as Vector, box_intersection


class Display:
    def __init__(self, canvas = None):
        if not canvas:
            canvas = ROOT.TCanvas()
        self.canvas = canvas
        a = Vector(1,2,3)
        b = Vector(11,22,33)
        self.box = ROOT.TPolyLine3D()
        self.box.SetNextPoint(a.x,a.y,a.z)
        self.box.SetNextPoint(b.x,b.y,b.z)
        self.box.SetLineColor(7)
        self.box.Draw()

        off = Vector(2,3,4)
        ray = Vector(1,1,2).norm()
        hit1 = Vector()
        hit2 = Vector()
        hitmask = box_intersection(a, b, off, ray, hit1, hit2);
        
        print 'hitmask %d' % hitmask
        print 'hit1 %f %f %f' % tuple([hit1[i] for i in range(3)])
        print 'hit2 %f %f %f' % tuple([hit2[i] for i in range(3)])

        self.bucket = list()

        if hitmask | 1:
            l1 = ROOT.TPolyLine3D()
            l1.SetNextPoint(off.x, off.y, off.z)
            l1.SetNextPoint(hit1.x, hit1.y, hit1.z)
            l1.SetLineColor(2)
            l1.Draw()
            self.bucket.append(l1)

        if hitmask | 2:
            l2 = ROOT.TPolyLine3D()
            l2.SetNextPoint(off.x, off.y, off.z)
            l2.SetNextPoint(hit2.x, hit2.y, hit2.z)
            l2.SetLineColor(4)
            l2.Draw()
            self.bucket.append(l2)


