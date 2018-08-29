pragma solidity ^0.4.24;

import "./Numerology.sol";

contract Verifier {

  function test_proof_verification() public returns (bool) {
    // Verifier.deployed().then(function(inst) { return inst.test_proof_verification.estimateGas(); })

    int256[4] memory k_l = [int256(-89243190524605339210527649141408088119), int256(-53877858828609620138203152946894934485), int256(-185204247857117235934281322466442848518), int256(-7585701889390054782280085152653861472)];

    uint256[4] memory P_Q = [0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798, 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8, 0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5, 0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a];

    uint256[4] memory wnaf;
    uint256 max_count = 0;
    uint256 count = 0;        

    for(uint j=0; j<4; j++){
        (wnaf[j], count) = Numerology._wnaf(k_l[j]);
        if(count > max_count){
          max_count = count;
        }
    }
    
    uint256[3] memory kP_lQ = Numerology._sim_mul_wnaf(wnaf, max_count, P_Q);

    kP_lQ = Numerology._sim_mul_wnaf(wnaf, max_count, P_Q);

    kP_lQ = Numerology._sim_mul_wnaf(wnaf, max_count, P_Q);

    uint256[3] memory expected = [0x7635e27fba8e1f779dcfdde1b1eacbe0571fbe39ecf6056d29ba4bd3ef5e22f2, 0x197888e5cec769ac2f1eb65dbcbc0e49c00a8cdf01f8030d8286b68c1933fb18, 1];

    return Numerology.eq_jacobian(kP_lQ, expected);

  }

}

