def DFR2(Q, D_u,D_v, N, K, eta1,eta2,k1,k2):
    """
    Adjusted DFR calculation for 1-dimensional v (single polynomial)
    
    Args:
        k1: number of additional public keys summed
        Q: modulus
        P: compression parameter  
        N: ring dimension
        eta: noise parameter for CBD
    """
    
    # 1. sum(pk_i) * r ≈ (k1+1) * (η²/4) * N
    term1_per_coeff_variance = k1*eta1 / 2* k2*eta2/2 * N*K
    
    # 2. sum(s_i) * u  
    term2_per_coeff_variance = k1*eta1 / 2* k2*eta1/2 * N*K
    
    # 3. e' 
    term3_per_coeff_variance = eta2 / 2
    
    # 4. delta_v 
    compression_per_coeff_variance = (Q / (2.0 * D_v))**2 / 3.0
    
    # compress u
    compress_u=(Q / (2.0 * D_u))**2 / 3.0 * N*K*k1*eta1 / 2

    
    # all
    variance_per_coefficient = (term1_per_coeff_variance + 
                              term2_per_coeff_variance + 
                              term3_per_coeff_variance + 
                              compression_per_coeff_variance+compress_u)
    
    total_variance = variance_per_coefficient
    
    sigma = math.sqrt(total_variance)
    z = (Q / 4.0) / sigma
    tail_prob =  math.erfc(z / math.sqrt(2)) / 2
    print("Prob:=",math.log(2*N * tail_prob,2) )

DFR2(7681, 2**11,2**6, 256, 3, 2,2,3,3)
