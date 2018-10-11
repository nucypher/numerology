from logging import getLogger
from typing import List

from web3.contract import Contract

from blockchain.constants import (DEVELOPMENT_ETH_AIRDROP_AMOUNT)
from blockchain.interfaces import BlockchainDeployerInterface


class TesterBlockchain:
    """A view of a blockchain through a provided interface"""

    _instance = None
    __default_interface_class = BlockchainDeployerInterface

    class ConnectionNotEstablished(RuntimeError):
        pass

    def __init__(self,
                 interface: BlockchainDeployerInterface = None,
                 airdrop=True) -> None:

        self.log = getLogger("test-blockchain")  # type: Logger

        # Default interface
        if interface is None:
            interface = self.__default_interface_class()
        self.__interface = interface

        # Singleton
        if self._instance is None:
            TesterBlockchain._instance = self
        else:
            raise RuntimeError("Connection already established - Use .connect()")

        if airdrop is True:  # ETH for everyone!
            self.ether_airdrop(amount=DEVELOPMENT_ETH_AIRDROP_AMOUNT)

    def __repr__(self):
        class_name = self.__class__.__name__
        r = "{}(interface={})"
        return r.format(class_name, self.__interface)

    @property
    def interface(self) -> BlockchainDeployerInterface:
        return self.__interface

    def get_contract(self, name: str) -> Contract:
        """
        Gets an existing contract from the registry, or raises UnknownContract
        if there is no contract data available for the name/identifier.
        """
        return self.__interface.get_contract_by_name(name)

    def wait_for_receipt(self, txhash: str, timeout: int = None) -> dict:
        """Wait for a transaction receipt and return it"""
        timeout = timeout if timeout is not None else self.interface.timeout
        result = self.__interface.w3.eth.waitForTransactionReceipt(txhash, timeout=timeout)
        return result

    def ether_airdrop(self, amount: int) -> List[str]:
        """Airdrops ether from creator address to all other addresses!"""

        coinbase, *addresses = self.interface.w3.eth.accounts

        tx_hashes = list()
        for address in addresses:

            tx = {'to': address, 'from': coinbase, 'value': amount}
            txhash = self.interface.w3.eth.sendTransaction(tx)

            _receipt = self.wait_for_receipt(txhash)
            tx_hashes.append(txhash)
            self.log.info("Airdropped {} ETH {} -> {}".format(amount, tx['from'], tx['to']))

        return tx_hashes

    def time_travel(self, hours: int=None, seconds: int=None):
        """
        Wait the specified number of wait_hours by comparing
        block timestamps and mines a single block.
        """

        more_than_one_arg = sum(map(bool, (hours, seconds))) > 1
        if more_than_one_arg:
            raise ValueError("Specify hours or seconds, not a combination")

        if hours:
            duration = hours * (60*60)
            base = 60 * 60
        elif seconds:
            duration = seconds
            base = 1
        else:
            raise ValueError("Specify either hours, seconds, or lock_periods.")

        now = self.interface.w3.eth.getBlock(block_identifier='latest').timestamp
        end_timestamp = ((now+duration)//base) * base

        self.interface.w3.eth.web3.testing.timeTravel(timestamp=end_timestamp)
        self.interface.w3.eth.web3.testing.mine(1)
        self.log.info("Time traveled to {}".format(end_timestamp))
