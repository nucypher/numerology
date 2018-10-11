

class InMemoryEthereumContractRegistry:

    class RegistryError(Exception):
        pass

    class UnknownContract(RegistryError):
        pass

    class IllegalRegistrar(RegistryError):
        """Raised when invalid data is encountered in the registry"""

    def __init__(self) -> None:
        self._registry_data = []

    def enroll(self, contract_name, contract_address, contract_abi):
        """
        Enrolls a contract to the chain registry by writing the name, address,
        and abi information to the filesystem as JSON.

        Note: Unless you are developing NuCypher, you most likely won't ever
        need to use this.
        """
        contract_data = [contract_name, contract_address, contract_abi]
        self._registry_data.append(contract_data)

    def search(self, contract_name: str=None, contract_address: str=None):
        """
        Searches the registry for a contract with the provided name or address
        and returns the contracts.
        """
        if not (bool(contract_name) ^ bool(contract_address)):
            raise ValueError("Pass contract_name or contract_address, not both.")

        contracts = list()
        for name, addr, abi in self._registry_data:
            if contract_name == name or contract_address == addr:
                contracts.append((name, addr, abi))

        if not contracts:
            raise self.UnknownContract
        if contract_address and len(contracts) > 1:
            m = "Multiple records returned for address {}"
            raise self.IllegalRegistrar(m.format(contract_address))

        return contracts if contract_name else contracts[0]
