const FastSecp256k1 = artifacts.require('FastSecp256k1.sol');
var Verifier = artifacts.require("./Verifier.sol");

module.exports = function(deployer) {
	deployer.deploy(FastSecp256k1);
 // deployer.deploy(FastSecp256k1).then(() => {
 //  		deployer.deploy(Verifier, {gas: 5000000});
 //  	});
 //  	deployer.link(FastSecp256k1, Verifier);
};
