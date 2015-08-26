try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from distutils.core import Command
from distutils.dir_util import mkpath, copy_tree
from os import path, chdir
import subprocess
import shutil

def opesci_dir(subdir=''):
    root_dir = path.dirname(path.realpath(__file__))
    return path.join(root_dir, subdir)

class CMakeBuilder(Command):

    user_options = []

    def build_with_cmake(self):
        mkpath(self._build_dir)
        chdir(self._build_dir)
        cmake_cmd = ["cmake", self._root_dir]
        if subprocess.call(cmake_cmd) != 0:
            raise EnvironmentError("error calling cmake")

        if subprocess.call("make") != 0:
            raise EnvironmentError("error calling make")
        chdir(self._root_dir)

    def initialize_options(self):
        self._root_dir = opesci_dir()
        self._build_dir = opesci_dir('build')

    def finalize_options(self):
        pass

    def build_lib(self):
        self.build_with_cmake()

    def run(self):
        self.build_with_cmake()

        # Copy lib and bin directories back to root
        copy_tree(opesci_dir('build/lib'), opesci_dir('lib'))
        copy_tree(opesci_dir('build/bin'), opesci_dir('bin'))


setup(name='opesci',
      version = "0.1", # Need better version management
      description = """Framework for automatic generation of
      Finite Difference models for use in seismic imaging.""",
      author = "Imperial College London and SENAI Cimatec",
      author_email = "opesci@imperial.ac.uk",
      url = "http://opesci.org",
      packages = ['opesci'],
      package_data = {'opesci' : ['templates/*.txt', 'templates/staggered/*.txt']},
      cmdclass = {'build_clib': CMakeBuilder},
)
