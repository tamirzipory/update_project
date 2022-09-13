from setuptools import Extension, setup

module = Extension("myspkmeans",
                   sources=[
                       'spkmeans.c',
                       'spkmeansmodule.c'
                   ])
setup(name='myspkmeans',
      version='1.0',
      author="Adar",
      description='project',
      ext_modules=[module])
