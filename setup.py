import sys

required_verion = (3,7)
if sys.version_info < required_verion:
    raise ValueError('atm-py needs at least python {}! You are trying to install it under python {}'.format('.'.join(str(i) for i in required_verion), sys.version))

# import ez_setup
# ez_setup.use_setuptools()

from setuptools import setup, find_packages
setup(
    name="atm-py",
    version="0.1",
    packages=find_packages(),
    author="Hagen Telg",
    author_email="hagen@hagnet.net",
    description="This package contains atmospheric science tools",
    license="MIT",
    keywords="atmospheric science tools",
    url="https://github.com/hagne/atm-py",
    install_requires=['numpy','pandas', 'scipy', 'matplotlib'],#pysolar
    # extras_require={'plotting': ['matplotlib'],},
    test_suite='nose.collector',
    tests_require=['nose'],
)