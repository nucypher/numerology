pragma solidity ^0.4.24;


/// @title FastSecp256k1: A Soidity library for fast ECC arithmetics using curve secp256k1
/// @author David Nuñez (david@nucypher.com)
library FastSecp256k1 {

    uint256 constant field_order = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

    /// @notice Equality test of two points in Jacobian coordinates
    /// @param P An EC point in Jacobian coordinates
    /// @param Q An EC point in Jacobian coordinates
    /// @return true if P and Q represent the same point in affine coordinates; false otherwise
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

    /// @notice Addition of two points in Jacobian coordinates
    /// @dev Based on the addition formulas from http://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2001-b.op3
    /// @param P An EC point in Jacobian coordinates
    /// @param Q An EC point in Jacobian coordinates
    /// @return An EC point in Jacobian coordinates with the sum , represented by an array of 3 uint256
    function addJac(uint[3] memory P, uint[3] memory Q) internal constant returns (uint[3] memory R) {

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
            return doubleJac(P);
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

    /// @notice Addition of two points in Jacobian coordinates, placing the result in the first point
    /// @dev Based on the addition formulas from http://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2001-b.op3
    /// @param P An EC point in Jacobian coordinates. The result is returned here.
    /// @param Q An EC point in Jacobian coordinates
    function addJacMutates(uint[3] memory P, uint[3] memory Q) internal constant {

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

        uint256 zz = mulmod(Pz, Pz, p);
        uint256 t0  = mulmod(Q[0], zz, p);
        uint256 t1  = mulmod(Q[1], mulmod(Pz, zz, p), p);

        zz = mulmod(Qz, Qz, p);
        uint256 a   = mulmod(P[0], zz, p);
        uint256 c   = mulmod(P[1], mulmod(Qz, zz, p), p);   
        

        if ((a == t0) && (c == t1)){
            doubleMutates(P);
            return;
        }
        
        t1   = addmod(t1, p-c, p); // d = t1 - c
        uint256 b = addmod(t0, p-a, p); // b = t0 - a
        uint256 e = mulmod(b, b, p); // e = b^2
        t0 = mulmod(a, e, p);    // t0 is actually "g"
        e = mulmod(e, b, p);  // f = b^3  (we will re-use the variable e )
        uint256 temp = addmod(mulmod(t1, t1, p), p-addmod(mulmod(2, t0, p), e, p), p);
        P[0] = temp;
        temp = mulmod(t1, addmod(t0, p-temp, p), p);
        P[1] = addmod(temp, p-mulmod(c, e, p), p);
        P[2] = mulmod(b, mulmod(Pz, Qz, p), p);
    }

    /// @notice Subtraction of two points in Jacobian coordinates, placing the result in the first point
    /// @dev Based on the addition formulas from http://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2001-b.op3
    /// @param P An EC point in Jacobian coordinates. The result is returned here.
    /// @param Q An EC point in Jacobian coordinates
    function subJacMutates(uint[3] memory P, uint[3] memory Q) internal constant {

        uint256 Pz = P[2];
        uint256 Qz = Q[2];
        uint256 p = field_order;

        if(Pz == 0){
            P[0] = Q[0];
            P[1] = p - Q[1];
            P[2] = Qz;
            return;
        } else if(Qz == 0){
            return;
        }

        uint256 zz = mulmod(Pz, Pz, p);
        uint256 t0  = mulmod(Q[0], zz, p);
        uint256 t1  = mulmod(p - Q[1], mulmod(Pz, zz, p), p);

        zz = mulmod(Qz, Qz, p);
        uint256 a   = mulmod(P[0], zz, p);
        uint256 c   = mulmod(P[1], mulmod(Qz, zz, p), p); 

        if ((a == t0) && (c == t1)){
            P[2] = 0;
            return;
        }
        
        t1   = addmod(t1, p-c, p); // d = t1 - c
        uint256 b = addmod(t0, p-a, p); // b = t0 - a
        uint256 e = mulmod(b, b, p); // e = b^2
        t0 = mulmod(a, e, p);    // t0 is actually "g"
        e = mulmod(e, b, p);  // f = b^3  (we will re-use the variable e )
        uint256 temp = addmod(mulmod(t1, t1, p), p-addmod(mulmod(2, t0, p), e, p), p);
        P[0] = temp;
        temp = mulmod(t1, addmod(t0, p-temp, p), p);
        P[1] = addmod(temp, p-mulmod(c, e, p), p);
        P[2] = mulmod(b, mulmod(Pz, Qz, p), p);
    }

    // Point doubling, 2*P
    // Params: Px, Py, Pz
    // Not concerned about the 1 extra mulmod.
    function doubleJac(uint[3] memory P) internal constant returns (uint[3] memory Q) {
        uint256 z = P[2];
        if (z == 0)
            return;
        uint256 p = field_order;
        uint256 Px = P[0];
        uint256 Py = P[1];
        uint256 Py2 = mulmod(Py, Py, p);
        uint256 s = mulmod(4, mulmod(Px, Py2, p), p);
        uint256 m = mulmod(3, mulmod(Px, Px, p), p);
        uint256 Qx = addmod(mulmod(m, m, p), p - addmod(s, s, p), p);
        Q[0] = Qx;
        Q[1] = addmod(mulmod(m, addmod(s, p - Qx, p), p), p - mulmod(8, mulmod(Py2, Py2, p), p), p);
        Q[2] = mulmod(2, mulmod(Py, z, p), p);
    }

    // Same as double but mutates P and is internal only.
    function doubleMutates(uint[3] memory P) internal constant {
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

        uint256[3][4] memory iPj;
        uint256[3] memory double;

        // P1 Lookup Table
        iPj = iP[0];
        iPj[0] = [P_Q[0], P_Q[1], 1];  // P1
        
        double = doubleJac(iPj[0]);
        iPj[1] = addJac(double, iPj[0]);
        iPj[2] = addJac(double, iPj[1]);
        iPj[3] = addJac(double, iPj[2]);

        // P1 Lookup Table
        iPj = iP[1];
        iPj[0] = [mulmod(beta, P_Q[0], p), P_Q[1], 1];    // P2

        double = doubleJac(iPj[0]);
        iPj[1] = addJac(double, iPj[0]);
        iPj[2] = addJac(double, iPj[1]);
        iPj[3] = addJac(double, iPj[2]);

        // Q1 Lookup Table
        iPj = iP[2];
        iPj[0] = [P_Q[2], P_Q[3], 1];                     // Q1

        double = doubleJac(iPj[0]);
        iPj[1] = addJac(double, iPj[0]);
        iPj[2] = addJac(double, iPj[1]);
        iPj[3] = addJac(double, iPj[2]);

        // Q2 Lookup Table
        iPj = iP[3];
        iPj[0] = [mulmod(beta, P_Q[2], p), P_Q[3], 1];    // Q2

        double = doubleJac(iPj[0]);
        iPj[1] = addJac(double, iPj[0]);
        iPj[2] = addJac(double, iPj[1]);
        iPj[3] = addJac(double, iPj[2]);
    }

    function _wnaf(int256 d) internal constant returns (uint256 ptr, uint256 length){
    
        int sign = d < 0 ? -1 : int(1);
        uint256 k = uint256(sign * d);

        length = 0;
        assembly
        {
            let ki := 0
            ptr := mload(0x40) // Get free memory pointer
            mstore(0x40, add(ptr, 300)) // Updates free memory pointer to +512 bytes offset
            for { } gt(k, 0) { } { // while k > 0
                if and(k, 1) {  // if k is odd:
                    ki := mod(k, 16)
                    k := add(sub(k, ki), mul(gt(ki, 8), 16))
                    // if sign = 1, store ki; if sign = -1, store 16 - ki
                    mstore8(add(ptr, length), add(mul(ki, sign), sub(8, mul(sign, 8))))
                }
                length := add(length, 1)
                k := div(k, 2)
            }
            //log3(ptr, 1, 0xfabadaacabada, d, length)    
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
        uint256[3][4][4] memory iP;
        _lookup_sim_mul(iP, P_Q);

        // LOOP 
        uint256 i = length;
        uint256 ki;
        uint256 ptr;
        while(i > 0) {
            i--;

            doubleMutates(Q);

            ptr = wnaf_ptr[0] + i;
            assembly {
                ki := byte(0, mload(ptr))
            }

            if (ki > 8) {
                subJacMutates(Q, iP[0][(15 - ki) / 2]);
            } else if (ki > 0) {
                addJacMutates(Q, iP[0][(ki - 1) / 2]);
            }

            ptr = wnaf_ptr[1] + i;
            assembly {
                ki := byte(0, mload(ptr))
            }

            if (ki > 8) {
                subJacMutates(Q, iP[1][(15 - ki) / 2]);
            } else if (ki > 0) {
                addJacMutates(Q, iP[1][(ki - 1) / 2]);
            } 

            ptr = wnaf_ptr[2] + i;
            assembly {
                ki := byte(0, mload(ptr))
            }

            if (ki > 8) {
                subJacMutates(Q, iP[2][(15 - ki) / 2]);
            } else if (ki > 0) {
                addJacMutates(Q, iP[2][(ki - 1) / 2]);
            } 

            ptr = wnaf_ptr[3] + i;
            assembly {
                ki := byte(0, mload(ptr))
            }

            if (ki > 8) {
                subJacMutates(Q, iP[3][(15 - ki) / 2]);
            } else if (ki > 0) {
                addJacMutates(Q, iP[3][(ki - 1) / 2]);
            } 
            
        }
    }

}