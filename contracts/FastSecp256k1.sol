/**
 * @title FastSecp256k1
 *
 * Fast ECC arithmetics for curve secp256k1
 *
 * @author David Nuñez (david@nucypher.com)
 */
library FastSecp256k1 {

    // Field order
    uint256 constant field_order = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

    // // Base point (generator) G
    // uint constant Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798;
    // uint constant Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8;

    // Order of G
    uint constant nn = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;


    // // Maximum value of s
    // uint constant lowSmax = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0;

    // For later
    //uint256 constant lambda = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72;
    //uint256 constant beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee;


    function eq_jacobian(uint256[3] memory P, uint256[3] memory Q) public returns(bool){
        uint p = field_order;

        if(P[2] == 0){
            return Q[2] == 0;   // P and Q are both zero.
        } else if(Q[2] == 0){
            return false;       // Q is zero but P isn't.
        }

        // Now we're sure none of them is zero

        uint256 Q_z_squared = mulmod(Q[2], Q[2], p);
        uint256 P_z_squared = mulmod(P[2], P[2], p);
        if (mulmod(P[0], Q_z_squared, p) != mulmod(Q[0], P_z_squared, p)){
          return false;
        }

        uint256 Q_z_cubed = mulmod(Q_z_squared, Q[2], p);
        uint256 P_z_cubed = mulmod(P_z_squared, P[2], p);
        return mulmod(P[1], Q_z_cubed, p) == mulmod(Q[1], P_z_cubed, p);
    
    }

    // Point addition, P + Q
    // inData: Px, Py, Pz, Qx, Qy, Qz
    // outData: Rx, Ry, Rz
    function _add2001b(uint[3] memory P, uint[3] memory Q) internal constant returns (uint[3] memory R) {

        /* 
        http://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2001-b.op3
      */

        if(P[2] == 0){
            return Q;
        } else if(Q[2] == 0){
            return P;
        }

        uint256 p = field_order;
        uint256 zz1 = mulmod(P[2], P[2], p);
        uint256 zz2 = mulmod(Q[2], Q[2], p);
        uint256 a   = mulmod(P[0], zz2, p);
        uint256 c   = mulmod(P[1], mulmod(Q[2], zz2, p), p);   
        uint256 t0  = mulmod(Q[0], zz1, p);
        uint256 t1  = mulmod(Q[1], mulmod(P[2], zz1, p), p);

        if ((a == t0) && (c == t1)){
            return _double(P);
        }
        uint256 d   = addmod(t1, p-c, p); // d = t1 - c
        uint256[3] memory b;
        b[0] = addmod(t0, p-a, p); // b = t0 - a
        b[1] = mulmod(b[0], b[0], p); // e = b^2
        b[2] = mulmod(b[1], b[0], p);  // f = b^3
        uint256 g = mulmod(a, b[1], p);
        R[0] = addmod(mulmod(d, d, p), p-addmod(mulmod(2, g, p), b[2], p), p);
        R[1] = addmod(mulmod(d, addmod(g, p-R[0], p), p), p-mulmod(c, b[2], p), p);
        R[2] = mulmod(b[0], mulmod(P[2], Q[2], p), p);
    }

    function _add2001bMutates(uint[3] memory P, uint[3] memory Q) internal constant {

        /* 
        http://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2001-b.op3
      */
        //(uint256 Qx, uint256 Qy, uint256 Qz) = (Q[0], Q[1], Q[2]);
        uint256 Pz = P[2];
        uint256 Qz = Q[2];

        if(Pz == 0){
            P[0] = Q[0];
            P[1] = Q[1];
            P[2] = Qz;
            return;
        } else if(Qz == 0){
            return;
        }

        uint256 p = field_order;

        uint256 zz1 = mulmod(Pz, Pz, p);
        uint256 Qzz = mulmod(Qz, Qz, p);
        uint256 a   = mulmod(P[0], Qzz, p);
        uint256 c   = mulmod(P[1], mulmod(Qz, Qzz, p), p);   
        uint256 t0  = mulmod(Q[0], zz1, p);
        uint256 t1  = mulmod(Q[1], mulmod(Pz, zz1, p), p);

        if ((a == t0) && (c == t1)){
            uint256[3] memory R = [P[0], P[1], Pz];
            _doubleM_jarl(R);
            return;
        }
        t1   = addmod(t1, p-c, p); // d = t1 - c
        uint256[3] memory b;
        b[0] = addmod(t0, p-a, p); // b = t0 - a
        b[1] = mulmod(b[0], b[0], p); // e = b^2
        b[2] = mulmod(b[1], b[0], p);  // f = b^3
        t0 = mulmod(a, b[1], p);    // t0 is actually "g"
        P[0] = addmod(mulmod(t1, t1, p), p-addmod(mulmod(2, t0, p), b[2], p), p);
        uint256 jarl = mulmod(t1, addmod(t0, p-P[0], p), p);
        P[1] = addmod(jarl, p-mulmod(c, b[2], p), p);
        P[2] = mulmod(b[0], mulmod(Pz, Qz, p), p);
    }

    // Same as addMixed but params are different and mutates P.
    function _addMixedM2001b(uint[3] memory P, uint[2] memory Q) internal constant {
        if(P[1] == 0) {
            P[0] = Q[0];
            P[1] = Q[1];
            P[2] = 1;
            return;
        }
        if(Q[1] == 0)
            return;

        uint256 p = field_order;
        uint256 zz1  = mulmod(P[2], P[2], p);
        uint256 zzz1 = mulmod(P[2], zz1, p);
        uint256 t0   = mulmod(Q[0], zz1, p);
        uint256 t1   = mulmod(Q[1], zzz1, p);
        uint256 b    = addmod(t0, p-P[0], p); // h
        uint256 d    = addmod(t1, p-P[1], p); // r

        uint256 e    = mulmod(b, b, p); // h2
        uint256 f    = mulmod(b, e, p); //  h3

        uint256 g    = mulmod(P[0], e, p);
        P[0]   = addmod(mulmod(d, d, p), p-addmod(mulmod(2, g, p), f, p), p);
        P[1]   = addmod(mulmod(d, addmod(g, p-P[0], p), p), p-mulmod(P[1], f, p), p);
        P[2]   = mulmod(b, P[2], p);
    }

    // Point doubling, 2*P
    // Params: Px, Py, Pz
    // Not concerned about the 1 extra mulmod.
    function _double(uint[3] memory P) internal constant returns (uint[3] memory Q) {
        uint256 p = field_order;
        if (P[2] == 0)
            return;
        uint Px = P[0];
        uint Py = P[1];
        uint Py2 = mulmod(Py, Py, p);
        uint s = mulmod(4, mulmod(Px, Py2, p), p);
        uint m = mulmod(3, mulmod(Px, Px, p), p);
        var Qx = addmod(mulmod(m, m, p), p - addmod(s, s, p), p);
        Q[0] = Qx;
        Q[1] = addmod(mulmod(m, addmod(s, p - Qx, p), p), p - mulmod(8, mulmod(Py2, Py2, p), p), p);
        Q[2] = mulmod(2, mulmod(Py, P[2], p), p);
    }

    // Same as double but mutates P and is internal only.
    function _doubleM_jarl(uint[3] memory P) internal constant {
        uint256 z = P[2];
        if (z == 0)
            return;
        uint256 p = field_order;
        uint256 x = P[0];
        uint256 _2y = mulmod(2, P[1], p);
        uint256 _4yy = mulmod(_2y, _2y, p);
        uint256 s = mulmod(_4yy, x, p);
        uint256 m = mulmod(3, mulmod(x, x, p), p);
        uint256 t = addmod(mulmod(m, m, p), mulmod(0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2d, s, p),p);
        P[0] = t;

        P[1] = addmod(mulmod(m, addmod(s, p - t, p), p), mulmod(0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffff7ffffe17, mulmod(_4yy, _4yy, p), p), p);

        P[2] = mulmod(_2y, z, p);
    }
    
    function _lookup_sim_mul(uint256[3][4][4] memory iP, uint256[4] memory P_Q) internal constant {
        uint256 p = field_order;
        uint256 beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee;
        iP[0][0] = [P_Q[0], P_Q[1], 1];                     // P1
        iP[1][0] = [mulmod(beta, P_Q[0], p), P_Q[1], 1];    // P2
        iP[2][0] = [P_Q[2], P_Q[3], 1];                     // Q1
        iP[3][0] = [mulmod(beta, P_Q[2], p), P_Q[3], 1];    // Q2

        // Lookup table
        uint256[3] memory double = _double(iP[0][0]);
        iP[0][1] = _add2001b(double, iP[0][0]);
        iP[0][2] = _add2001b(double, iP[0][1]);
        iP[0][3] = _add2001b(double, iP[0][2]);

        double = _double(iP[1][0]);
        iP[1][1] = _add2001b(double, iP[1][0]);
        iP[1][2] = _add2001b(double, iP[1][1]);
        iP[1][3] = _add2001b(double, iP[1][2]);

        double = _double(iP[2][0]);
        iP[2][1] = _add2001b(double, iP[2][0]);
        iP[2][2] = _add2001b(double, iP[2][1]);
        iP[2][3] = _add2001b(double, iP[2][2]);

        double = _double(iP[3][0]);
        iP[3][1] = _add2001b(double, iP[3][0]);
        iP[3][2] = _add2001b(double, iP[3][1]);
        iP[3][3] = _add2001b(double, iP[3][2]);
    }

    function _wnaf(int256 k) internal constant returns (uint256 ptr, uint256 length){
        uint256 d;
        uint256 neg = 0;
        length = 0;

        if(k < 0){
            neg = 1;
            k = -k;
        } else {
            neg = 0;
        }
        d = uint256(k);

        assembly
        {
                let dm := 0
                let dms := 0
                ptr := mload(0x40) // Get free memory pointer
                mstore(0x40, add(ptr, 512)) // Updates free memory pointer to +512 bytes offset
            loop:
                jumpi(loop_end, iszero(d))
                jumpi(even, iszero(and(d, 1)))
                dm := mod(d, 16)
                dms := dm
                if neg {
                    dms := sub(16, dm)
                }
                mstore8(add(ptr, length), dms) // Don't store as signed - convert when reading.
                d := add(sub(d, dm), mul(gt(dm, 8), 16))
            even:
                d := div(d, 2)
                length := add(length, 1)
                jump(loop)
            loop_end:
                //log3(ptr, 1, 0xfabadaacabada, k, length)    
        }

        return (ptr, length);
    }

    function _sim_mul(int256[4] memory k_l, uint256[4] memory P_Q) internal constant returns (uint[3] memory Q) {
        uint256[4] memory wnaf;
        uint256 max_count = 0;
        uint256 count = 0;        

        for(uint j=0; j<4; j++){
            (wnaf[j], count) = _wnaf(k_l[j]);
            if(count > max_count){
                max_count = count;
            }
        }

        Q = _sim_mul_wnaf(wnaf, max_count, P_Q);
    }

    function _sim_mul_wnaf(uint256[4] memory wnaf_ptr, uint256 length, uint256[4] memory P_Q) internal constant returns (uint[3] memory Q) {
        uint256 p = field_order;

        uint256[3][4][4] memory iP;
        _lookup_sim_mul(iP, P_Q);

        // LOOP 
        uint256 i = length;
        while(i > 0) {

            _doubleM_jarl(Q);

            uint dj;
            uint pIdx;

            i--;
            uint256 ptr;
            for(uint j=0; j<4; j++){
                ptr = wnaf_ptr[j];
                assembly {
                    dj := byte(0, mload(add(ptr, i)))
                }

                if (dj > 8) {
                    pIdx = (15 - dj) / 2; // These are the "negative ones", so invert y.
                    _add2001bMutates(Q, [iP[j][pIdx][0], p - iP[j][pIdx][1], iP[j][pIdx][2]]);
                } else if (dj > 0) {
                    pIdx = (dj - 1) / 2;
                    _add2001bMutates(Q, iP[j][pIdx]);
                } 
            }
        }
    }

}