# Numerology 
*Fast ECC arithmetic library for secp256k1 in Solidity*

**WARNING: Don't use it with any private input (e.g., a private key), as it's going to be exposed. Don't ever use this to sign something! Also, this is experimental -- Use at your own risk!**

**Numerology** is a Solidity library for elliptic curve arithmetic on secp256k1. Its original goal is to support [NuCypher](https://github.com/nucypher/nucypher) smart contracts for Re-Encryption Verification of [Umbral](https://github.com/nucypher/pyUmbral) ciphertexts, which is done by means of Non-Interactive Zero Knowledge Proofs of Knowledge. In particular, Numerology provides optimized ECC algorithms useful for verifying Schnorr-like NIZKs. These algorithms are honed to spend as little gas as possible. For this reason, and also to bypass some EVM and Solidity limitations, the code can be ugly as hell sometimes (sorry for that!).

### Features:
- Point addition and doubling in Jacobian coordinates with optimized formulas.
- Simultaneous scalar multiplication with several optimizations (interleaving, WNAF form, curve endomorphism).
- Verifying an equation of the form `aP + bQ = R` takes approximately 450,000 gas.




