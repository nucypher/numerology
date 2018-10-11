from distutils.core import setup

INSTALL_REQUIRES = []
TESTS_REQUIRE = []

setup(name='numerology',
      version='0.1',
      description='Optimized ECC arithmetic library for secp256k1 in Solidity',
      extras_require={'testing': TESTS_REQUIRE},
      install_requires=INSTALL_REQUIRES,
      packages=['numerology'])
