from distutils.core import setup

setup(name='atmPy',
      version='0.1',
      description='Python Distribution Utilities',
      author='Hagen Telg and Matt Richardson',
      author_email='matt.richardson@msrconsults.com',
      packages=['atmPy'], requires=['numpy', 'pandas', 'matplotlib']
      )
