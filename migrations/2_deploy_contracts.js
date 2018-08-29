const Numerology = artifacts.require('Numerology.sol');
var Verifier = artifacts.require("Verifier.sol");

module.exports = function(deployer) {
	deployer.deploy(Numerology);
  	deployer.link(Numerology, Verifier);
  	deployer.deploy(Verifier);
};
