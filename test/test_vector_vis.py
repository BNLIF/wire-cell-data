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

        self.frame = self.canvas.DrawFrame(a.z-1, a.y-1, b.z+1, b.y+1)
        self.frame.SetXTitle('Z direction')
        self.frame.SetYTitle('Y direction')
        self.bucket = list()

        rays = [Vector(0,1,1).norm(), Vector(0,1,-1).norm(), Vector(0,1,0)]
        nwires = 10
        for iz in range(nwires):
            for ray in rays:
                z = float(iz)/nwires * (b.z-a.z) + a.z
                off = Vector(a.x+.1, 0.5*(a.y+b.y), z+.1)
                hit1 = Vector()
                hit2 = Vector()
                print '\n off %4.1f %4.1f %4.1f' % tuple([off[i] for i in range(3)])
                print ' ray %4.1f %4.1f %4.1f' % tuple([ray[i] for i in range(3)])
                hitmask = box_intersection(a, b, off, ray, hit1, hit2);
                print 'hitmask %d' % hitmask


                if hitmask & 1:
                    print 'hit1 %4.1f %4.1f %4.1f' % tuple([hit1[i] for i in range(3)])
                    l1 = ROOT.TLine(off.z,off.y, hit1.z,hit1.y)
                    l1.SetLineColor(2)
                    l1.Draw()
                    self.bucket.append(l1)

                if hitmask & 2:
                    print 'hit2 %4.1f %4.1f %4.1f' % tuple([hit2[i] for i in range(3)])
                    l1 = ROOT.TLine(off.z,off.y, hit2.z,hit2.y)
                    l1.SetLineColor(4)
                    l1.Draw()
                    self.bucket.append(l1)
        self.canvas.Update()

