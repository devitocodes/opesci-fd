try:
    from setuptools import setup
    from setuptools.command.build import build
    from setuptools.command.build_clib import build_clib
except ImportError:
    from distutils.core import setup
    from distutils.command.build import build
    from distutils.command.build_clib import build_clib
from os import path, chdir
import subprocess
import shutil
import versioneer


class CMakeBuilder(build_clib):

    def run(self):
        root_dir = path.dirname(path.realpath(__file__))
        build_dir = path.realpath(self.build_temp)
        clib_dir = path.realpath(self.build_clib)

        # Run cmake and make in temp dir
        self.mkpath(build_dir)
        chdir(build_dir)
        cmake_cmd = ["cmake", root_dir]
        if subprocess.call(cmake_cmd) != 0:
            raise EnvironmentError("error calling cmake")

        if subprocess.call("make") != 0:
            raise EnvironmentError("error calling make")
        chdir(root_dir)

        # Copy lib and include directories to opesci package
        self.copy_tree(path.join(root_dir, 'include'),
                       path.join(clib_dir, 'opesci/include'))
        self.copy_tree(path.join(build_dir, 'lib'),
                       path.join(clib_dir, 'opesci/lib'))


class PythonBuilder(build):
    def run(self):
        # Pass build_lib to the CMake builder as build_clib
        build_clib = self.get_finalized_command('build_clib')
        build_clib.build_clib = self.build_lib
        build_clib.run()

        # Invoke the actual Python build
        build.run(self)


setup(name='opesci',
      version = versioneer.get_version(),
      description = """Framework for automatic generation of
      Finite Difference models for use in seismic imaging.""",
      author = "Imperial College London and SENAI Cimatec",
      author_email = "opesci@imperial.ac.uk",
      url = "http://opesci.org",
      packages = ['opesci'],
      package_data = {'opesci' : ['lib/*.so', 'include/*.h',
                                  'templates/*.txt',
                                  'templates/staggered/*.txt',
                                  'templates/staggered/*.cpp']},
      cmdclass = {'build_clib': CMakeBuilder, 'build': PythonBuilder}
)
