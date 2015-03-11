#!/usr/bin/env python

def options(ctx):
    ctx.load('compiler_cxx')

def configure(ctx):
    ctx.load('compiler_cxx')
    ctx.load('find_root', tooldir='.')
#    ctx.check_cfg(path='root-config', args='--cflags --libs',
#                  package='', uselib_store='ROOT')
#    ctx.find_program('rootcint', var='ROOTCINT')


def build(bld):
    # main code library
    bld.shlib(source = bld.path.ant_glob('lib/*.cxx'),
              target = 'WireCell', 
              includes = ["inc"],
              use = 'ROOTSYS')


    # generate rootcint dictionary and build a shared library
    bld(source= bld.path.ant_glob('inc/*.h') + ["dict/LinkDef.h"],
        target='WireCellDict.cxx',
        rule='${ROOTCINT} -f ${TGT} -c ${SRC}')
    bld.shlib(source = 'WireCellDict.cxx',
              target = 'WireCellDict', 
              includes = ["inc"],
              use = 'ROOTSYS')


