#!/usr/bin/env python

APPNAME = 'WireCellData'

def options(ctx):
    ctx.load('find_package')

def configure(ctx):
    ctx.load('find_package')


def build(bld):
    # main code library
    bld.shared_library()
    bld.api_headers(name='WireCell')
    bld.root_dictionary(headers = 'inc/WireCell/*.h')

    for testsrc in bld.path.ant_glob('test/test_*.cxx'):
        bld.test_program([testsrc])


