from logging import getLogger
from typing import Tuple

from eth_tester import EthereumTester
from eth_tester import PyEVMBackend
from web3 import Web3
from web3.contract import Contract
from web3.providers.eth_tester.main import EthereumTesterProvider

from blockchain.compile import SolidityCompiler
from blockchain.constants import NUCYPHER_GAS_LIMIT
from blockchain.registry import InMemoryEthereumContractRegistry


class BlockchainDeployerInterface:
    """
    Interacts with a solidity compiler and a registry in order to instantiate compiled
    ethereum contracts with the given web3 provider backend.
    """
    __default_timeout = 10  # seconds
    # __default_transaction_gas_limit = 500000  # TODO: determine sensible limit and validate transactions

    class UnknownContract(Exception):
        pass

    class InterfaceError(Exception):
        pass

    def __init__(self,
                 timeout: int = None,
                 registry: InMemoryEthereumContractRegistry = None,
                 compiler: SolidityCompiler=None,
                 deployer_address: str = None) -> None:

        self.log = getLogger("blockchain-interface")                       # type: Logger

        pyevm_backend = PyEVMBackend.from_genesis_overrides(parameter_overrides={'gas_limit': NUCYPHER_GAS_LIMIT})
        eth_tester = EthereumTester(backend=pyevm_backend, auto_mine_transactions=True)
        self._provider = EthereumTesterProvider(ethereum_tester=eth_tester)
        providers = list()
        providers.append(self._provider)
        web3_instance = Web3(providers=providers)  # Instantiate Web3 object with provider
        self.w3 = web3_instance

        self.timeout = timeout if timeout is not None else self.__default_timeout
        self.__sol_compiler = compiler

        # Setup the registry and base contract factory cache
        registry = registry if registry is not None else InMemoryEthereumContractRegistry()
        self.registry = registry

        interfaces = self.__sol_compiler.compile()
        self.__raw_contract_cache = interfaces

        self.__deployer_address = deployer_address

    def get_contract_factory(self, contract_name: str) -> Contract:
        """Retrieve compiled interface data from the cache and return web3 contract"""
        try:
            interface = self.__raw_contract_cache[contract_name]
        except KeyError:
            raise self.UnknownContract('{} is not a locally compiled contract.'.format(contract_name))
        else:
            contract = self.w3.eth.contract(abi=interface['abi'],
                                            bytecode=interface['bin'],
                                            ContractFactoryClass=Contract)
            return contract

    def _wrap_contract(self, dispatcher_contract: Contract,
                       target_contract: Contract, factory=Contract) -> Contract:
        """Used for upgradeable contracts."""

        # Wrap the contract
        wrapped_contract = self.w3.eth.contract(abi=target_contract.abi,
                                                address=dispatcher_contract.address,
                                                ContractFactoryClass=factory)
        return wrapped_contract

    def get_contract_by_address(self, address: str):
        """Read a single contract's data from the registrar and return it."""
        try:
            contract_records = self.registry.search(contract_address=address)
        except RuntimeError:
            raise self.InterfaceError('Corrupted Registrar')  # TODO: Integrate with Registry
        else:
            if not contract_records:
                raise self.InterfaceError("No such contract with address {}".format(address))
            return contract_records[0]

    def get_contract_by_name(self, name: str, upgradeable=False, factory=Contract) -> Contract:
        """
        Instantiate a deployed contract from registrar data,
        and assemble it with it's dispatcher if it is upgradeable.
        """
        target_contract_records = self.registry.search(contract_name=name)

        if not target_contract_records:
            raise self.InterfaceError("No such contract records with name {}".format(name))

        if upgradeable:
            # Lookup dispatchers; Search fot a published dispatcher that targets this contract record
            dispatcher_records = self.registry.search(contract_name='Dispatcher')

            matching_pairs = list()
            for dispatcher_name, dispatcher_addr, dispatcher_abi in dispatcher_records:

                dispatcher_contract = self.w3.eth.contract(abi=dispatcher_abi,
                                                           address=dispatcher_addr,
                                                           ContractFactoryClass=factory)

                # Read this dispatchers target address from the blockchain
                live_target_address = dispatcher_contract.functions.target().call()

                for target_name, target_addr, target_abi in target_contract_records:
                    if target_addr == live_target_address:
                        pair = dispatcher_addr, target_abi
                        matching_pairs.append(pair)

            else:  # for/else

                if len(matching_pairs) == 0:
                    raise self.InterfaceError("No dispatcher targets known contract records for {}".format(name))

                elif len(matching_pairs) > 1:
                    raise self.InterfaceError("There is more than one dispatcher targeting {}".format(name))

                selected_contract_address, selected_contract_abi = matching_pairs[0]
        else:
            if len(target_contract_records) != 1:  # TODO: Allow multiple non-upgradeable records (UserEscrow)
                m = "Multiple records returned from the registry for non-upgradeable contract {}"
                raise self.InterfaceError(m.format(name))

            selected_contract_name, selected_contract_address, selected_contract_abi = target_contract_records[0]

        # Create the contract from selected sources
        unified_contract = self.w3.eth.contract(abi=selected_contract_abi,
                                                address=selected_contract_address,
                                                ContractFactoryClass=factory)

        return unified_contract

    @property
    def deployer_address(self):
        return self.__deployer_address

    @deployer_address.setter
    def deployer_address(self, checksum_address: str) -> None:
        if self.deployer_address is not None:
            raise RuntimeError("{} already has a deployer address set.".format(self.__class__.__name__))
        self.__deployer_address = checksum_address

    def deploy_contract(self, contract_name: str, *args, **kwargs) -> Tuple[Contract, str]:
        """
        Retrieve compiled interface data from the cache and
        return an instantiated deployed contract
        """
        if self.__deployer_address is None:
            raise self.InterfaceError('No deployer address is configured.')
        #
        # Build the deployment tx #
        #

        deploy_transaction = {'from': self.deployer_address, 'gasPrice': self.w3.eth.gasPrice}
        self.log.info("Deployer address is {}".format(deploy_transaction['from']))

        contract_factory = self.get_contract_factory(contract_name=contract_name)
        deploy_bytecode = contract_factory.constructor(*args, **kwargs).buildTransaction(deploy_transaction)
        self.log.info("Deploying contract: {}: {} bytes".format(contract_name, len(deploy_bytecode['data'])))

        #
        # Transmit the deployment tx #
        #
        txhash = contract_factory.constructor(*args, **kwargs).transact(transaction=deploy_transaction)
        self.log.info("{} Deployment TX sent : txhash {}".format(contract_name, txhash.hex()))

        # Wait for receipt
        receipt = self.w3.eth.waitForTransactionReceipt(txhash)
        address = receipt['contractAddress']
        self.log.info("Confirmed {} deployment: address {}".format(contract_name, address))

        #
        # Instantiate & enroll contract
        #
        contract = contract_factory(address=address)

        self.registry.enroll(contract_name=contract_name,
                             contract_address=contract.address,
                             contract_abi=contract_factory.abi)

        return contract, txhash
