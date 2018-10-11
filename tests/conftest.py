import pytest

from blockchain.chains import TesterBlockchain
from blockchain.compile import SolidityCompiler
from blockchain.constants import TEST_CONTRACTS_DIR
from blockchain.interfaces import BlockchainDeployerInterface
from blockchain.registry import InMemoryEthereumContractRegistry


@pytest.fixture(scope='session')
def solidity_compiler():
    """Doing this more than once per session will result in slower test run times."""
    compiler = SolidityCompiler(test_contract_dir=TEST_CONTRACTS_DIR)
    yield compiler


@pytest.fixture(scope='module')
def testerchain(solidity_compiler):
    _temp_registry = InMemoryEthereumContractRegistry()
    deployer_interface = BlockchainDeployerInterface(compiler=solidity_compiler,  # freshly recompile if not None
                                                     registry=_temp_registry)
    testerchain = TesterBlockchain(interface=deployer_interface)

    origin, *everyone = testerchain.interface.w3.eth.accounts
    deployer_interface.deployer_address = origin

    yield testerchain
