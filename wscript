#!/usr/bin/env python

APPNAME = 'WireCell'

def options(ctx):
    ctx.load('find_package')

def configure(ctx):
    ctx.load('find_package')


def build(bld):
    # main code library
    bld.shared_library()
    bld.api_headers()
    bld.rootcint_dictionary()



