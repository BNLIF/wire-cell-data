#!/usr/bin/env python

def options(ctx):
    ctx.load('compiler_cxx')

def configure(ctx):
    ctx.load('compiler_cxx')
    ctx.check_cfg(path='root-config', args='--cflags --libs',
                  package='', uselib_store='ROOT')
    ctx.find_program('rootcint', var='ROOTCINT')

def build(bld):
    bld.recurse('data dict')# test


