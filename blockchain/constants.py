import os
from os.path import abspath, dirname

import blockchain

NUCYPHER_GAS_LIMIT = 5000000

DEVELOPMENT_ETH_AIRDROP_AMOUNT = 10 ** 18  # wei -> ether

DEFAULT_NUMBER_OF_URSULAS_IN_DEVELOPMENT_NETWORK = 10

# Base Filepaths
BASE_DIR = abspath(dirname(dirname(blockchain.__file__)))
CONTRACTS_DIR = os.path.join(BASE_DIR, 'numerology', 'contracts')

# Test Constants
TEST_CONTRACTS_DIR = os.path.join(BASE_DIR, 'tests', 'contracts')
