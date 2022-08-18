from setuptools import Extension, setup

module = Extension("myspkmeans",
                   sources=[
                       'spkmeans.c',
                       'spkmeansmodule.c'
                   ])
setup(name='myspkmeans',
      version='1.0',
      author="Tamir Zipory and Yarin Diga",
      description='פרוייקט מזדיין',
      ext_modules=[module])
