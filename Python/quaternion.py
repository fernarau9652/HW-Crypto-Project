#######################################################################
# Developed by Elmir Dzaka and Fernando Araujo
# Verification of Post-Quantum Crypto Schemes Using Quaternion Algebra
#######################################################################

from pyquaternion import Quaternion
import random
import hashlib


def modulo(quaternion, q):
    # initalize array for new quaternion computation
    new_quaternion_arr = [];
    for comp in quaternion:
        new_quaternion_arr.append(comp%q);

    # create new quaternion from the arry and return
    new_quaternion = Quaternion(new_quaternion_arr[0], new_quaternion_arr[1], new_quaternion_arr[2], new_quaternion_arr[3])
    
    # return new_quaternion
    return new_quaternion

def vector_product_verify(Q, G, N, q):
    # Compute modulo of all Quaternion arithmetic
    Q_N = modulo(Q*N, q)
    N_Q = modulo(N*Q, q)

    G_N = modulo(G*N, q)
    N_G = modulo(N*G, q)

    Q_G = modulo(Q*G, q)
    G_Q = modulo(G*Q, q)

    # Verify non-invertibility
    if (Q_N != N_Q) and (G_N != N_G) and (Q_G != G_Q):
        print("(A*N != N*A) => ", Q_N, "!=", N_Q)
        print("(B*N != N*B) => ", G_N, "!=", N_G)
        print("(A*B != B*A) => ", Q_G, "!=", G_Q)
        return True
    else:
        return False

def verify_generator(N,q):
    p = (2*q) + 1
    l = 1
    u = 1
    
    for h in range(p):
        for d in range(p):
            # generate ln and ln_prime using formulas 5,22
            ln = [d,h,(N[1]*(1-(l*h)))/(u*N[3]), (N[0]*(1-(u*d)))/(l*N[2])]
            ln_prime = [d, (1-(u*d))/(l), (N[1]/N[3])*d, (N[0]*(1-(u*d)))/(l*N[2])]

            # generate rn and rn_prime using formulas 7,23
            rn = [d,h, (N[0]*(1-(u*d)))/(l*N[3]), (N[1]*(1 - (l*h)))/(u*N[2])]
            rn_prime = [d, (1-(u*d))/(l), (N[0]*(l-u*d))/(l*N[3]), (N[1]/N[2])*d]

            # create Quaternions for lb and ln_prime and perform modulo operation of q
            ln_q = Quaternion(ln[0],ln[1],ln[2],ln[3])
            ln_prime_q = Quaternion(ln_prime[0],ln_prime[1],ln_prime[2],ln_prime[3])
            ln_q = modulo(ln_q, q)
            ln_prime_q = modulo(ln_prime_q,q)
            
            # create Quaternions for rb and rn_prime and perform modulo operation of q
            rn_q = Quaternion(rn[0],rn[1],rn[2],rn[3])
            rn_prime_q = Quaternion(rn_prime[0],rn_prime[1],rn_prime[2],rn_prime[3])
            rn_q = modulo(rn_q, q)
            rn_prime_q = modulo(rn_prime_q,q)

            # verify if left sided units and right sided units equal their prime, and have same generators
            if(ln_q == ln_prime_q and rn_q == rn_prime_q):
                print("(d,h): (",d,",",h,")")
                print("------")
    return True

def find_non_inv_vec(N, q):
    # calulate E_double
    p = (2*q) + 1
    l = 1
    u = 1
    E_double_prime = [(N[0])/((l*N[1])+(u*N[0])), (N[1])/((l*N[1])+(u*N[0])), (N[2])/((l*N[1])+(u*N[0])), (N[3])/((l*N[1])+(u*N[0]))]

    # create quaternion from E_double_prime calculation
    E_double_prime_q = Quaternion(E_double_prime[0],E_double_prime[1],E_double_prime[2],E_double_prime[3])
    E_double_prime_q = modulo(E_double_prime_q,q)

    # find unquie non-inv vector that makes E_prime and E_double_prime equal
    for d in range(p):
        E_prime = [d, ((l*N[1])-(u * N[0]) + ((u ** 2)*N[0]*d))/((l ** 2) * N[1]), (N[0]*(1 - (u*d)))/(l*N[3]), (N[0]*(1-(u*d)))/(l * N[2])]
        E_prime_q = Quaternion(E_prime[0], E_prime[1], E_prime[2], E_prime[3])
        E_prime_q = modulo(E_prime_q,q)

        # check to see if unique non-invertble vector
        if((E_prime_q[0]*E_prime_q[1]) == (E_prime_q[2]*E_prime_q[3])):
            return d
    return -1


def left_sided_unit(A, q, d):
    # Initialize structural coefficients and left handed unit
    l = 1
    u = 1
    
    # Perform Algorithm
    Ln = []
    sub1 = d
    sub2 = (1/l) - ((u*d) / l) 
    sub3 = (A[1]/A[3]) * d
    sub4 = (A[0] * (1-(u*d)))/(l*A[2])
    Ln.append(sub1)
    Ln.append(sub2)
    Ln.append(sub3)
    Ln.append(sub4)
    
    # Create Quaternion and return
    new_quaternion = Quaternion(Ln[0], Ln[1], Ln[2], Ln[3])
    new_quaternion = modulo(new_quaternion,q)
    return new_quaternion


def right_sided_unit(A, q, d):
    l = 1
    u = 1

    # Perform right hand algorithm
    Ln = []
    sub1 = d
    sub2 = (1/l) - ((u*d) / l) 
    sub3 =  (A[0] * (1 - (u*d)))/(l*A[3])
    sub4 = (A[1]/A[2])*d
    Ln.append(sub1)
    Ln.append(sub2)
    Ln.append(sub3)
    Ln.append(sub4)
    
    # Create Quaternion
    new_quaternion = Quaternion(Ln[0], Ln[1], Ln[2], Ln[3])
    new_quaternion = modulo(new_quaternion,q)

    return new_quaternion

def calculate_e_prime(N,q,d):
    l = 1
    u = 1

    E_prime = [d, ((l*N[1])-(u * N[0]) + ((u ** 2)*N[0]*d))/((l ** 2) * N[1]), (N[0]*(1 - (u*d)))/(l*N[3]), (N[0]*(1-(u*d)))/(l * N[2])]
    E_prime_q = Quaternion(E_prime[0], E_prime[1], E_prime[2], E_prime[3])
    E_prime_q = modulo(E_prime_q,q)

    return E_prime_q

def compute_y(A, N, Ln, q, x):
    Y_arr = []

    # raise N to the power of x
    N_raised = N ** x
    sub1 = A*N_raised
    sub2 = sub1*Ln
    Y = sub2*A.inverse
    Y = modulo(Y,q)
    return Y

def compute_z(B, N, Rn, q):
    sub1 = B*Rn
    sub2 = sub1 * N
    Z = sub2 * B.inverse
    Z = modulo(Z,q)
    return Z

def compute_e_double(N, q):
    # structural coefficents
    l = 1
    u = 1

    # perform algorithm
    E_double = []
    sub1 = N[0]/((l*N[1]) + (u*N[0]))
    sub2 = N[1]/((l*N[1]) + (u*N[0]))
    sub3 = N[2]/((l*N[1]) + (u*N[0]))
    sub4 = N[3]/((l*N[1]) + (u*N[0]))

    # append operations to matric and create Quaternion
    E_double.append(sub1)
    E_double.append(sub2)
    E_double.append(sub3)
    E_double.append(sub4)

    quaternion = Quaternion(E_double[0], E_double[1], E_double[2], E_double[3])
    return quaternion

def compute_e_double_w(N, q):
    p = (2*q)+1
    E_double = N ** (p-1)
    return E_double

def compute_T(A, B, E_double, q):
    sub1 = B * E_double
    T = sub1 * A.inverse
    T = modulo(T,q)
    return T

def compute_j(B, rn, q):
    J = B * rn
    J = modulo(J, q)
    return J

def compute_u(A, ln, q):
    U = ln * A.inverse
    U = modulo(U, q)
    return U

def compute_V(k, j, N, u, q):
    V = (j*(N**k)) * u
    V = modulo(V,q)
    return V

def compute_F_h(M, V):
    # convert V into bytes
    V = str(V)
    V_bytes = bytes(V, "utf-8")

    # hash the message using SHA
    M_hashed = hashlib.sha256(M).digest()

    # concatenate hashed M and V to make collision resistant
    concatenated_data = V_bytes + M_hashed

    # hash again for same reason
    final_hash = hashlib.sha256(concatenated_data).digest() 

    # return an integer
    v = int.from_bytes(final_hash, byteorder = 'big')
    return v

def compute_s(k, x, v, q):
    s = (k+ (x*v)) % q
    return s

def main ():
    #####################################
    # KEY GENERATION ALGORITHM 
    #####################################
    print("################################################")
    print("KEY GENERATION ALGORITHM")
    print("################################################")

    # Creating a random Quaternions
    # Q = A
    # G = B
    A = Quaternion(1,2,3,4)
    B = Quaternion(2,3,4,5)

    # Non-Invertable Quaternion N
    N = Quaternion(2,1,1,2)

    # Print Non-commucative mutiplication
    #print(A*B)
    AB = A*B;
    
    # do modulus q operation
    q = 7;
    
    # A % q
    A_mod = modulo(A,q);
    #print(A_mod) 

    # A*B (mod q)
    AB_mod = modulo(AB,q)
    #print(AB_mod)

    # Step 1
    # Verify vectors
    print("-------------------------------------")
    print("STEP 1: Verifying Vector Correctness:")
    print(vector_product_verify(A,B,N,q))
    print("-------------------------------------")

    # step two: Make sure that left_sided unit and its prime have a generator given a prime order qi
    print("STEP 2: Verify generators d,h for ln and rn:")
    verify_generator(N, q)
    print("-------------------------------------")
    #print(verify_generator(ln,ln_prime,q))

    # step 3: find non_invertible vector geneated by d,h ln and rn loops. Only one vector is non-invertible within all ln and rn vectors. NOTE: # of vectors generated  by ln/rn formulas is = p-1 vectors
    print("STEP 3: FInd unique non-inv vector for ln and rn. Find unique d that generates the unique non-inv vec")
    d = find_non_inv_vec(N,q)
    print(d)
    print("Since unique d == 3, we can look at step 2 to find corresponding h. For this example: d == 3 and h == 12. Using these two generators, we input into formulas 22 and 23 to find unique ln,rn to generate public key vectors")
    print("-------------------------------------")

    # step 4: calculate unique public key vectors
    print("STEP 4: Calculate and store unique ln, rn, and en vectors for public key generation using generators found")
    d = 3
    ln = left_sided_unit(N, q, d)
    rn = right_sided_unit(N, q, d)
    e_prime = calculate_e_prime(N, q, d)

    # print out vectors to double check
    print("ln: ", ln)
    print("rn: ", rn)
    print("e_prime_n: ", e_prime)

    print("-------------------------------------")
    print("STEP 5: calculate public key vectors using ln, rn, en. Need random integer x < q")
    x = random.randint(0,q-1)
    y = compute_y(A, N, ln, q, x)
    z = compute_z(B, N, rn, q)
    t = compute_T(A, B, e_prime, q) 
    print("x: ", x)
    print("Y: ", y)
    print("Z: ", z)
    print("T: ", t)
    print("-------------------------------------")
    print("STEP 6: calculate private keys using ln,rn")
    j = compute_j(B, rn, q)
    u = compute_u(A, ln, q)
    print("J: ", j)
    print("U: ", u)
    print("-------------------------------------")

    #####################################################
    # SIGNITURE GENERATION ALGOIRTHM
    #####################################################
    print("################################################")
    print("SIGNITURE GENERATION ALGORITHM")
    print("################################################")
    print("-------------------------------------")
    print("STEP 1: select at random an integer k < q to compute signatures (v,s) NOTE: V and v are different here")
    k = random.randint(0,q-1)
    V = compute_V(k, j, N, u, q)
    print("k: ", k)
    print("V: ", V)
    print("-------------------------------------")
    print("STEP 2: Compute signiture element v = F_h(M,V) where M is the given message and F_h is a collision resistent hash function. For simplicity, we use SHA as our hash protocol")
    M = "Hello World"
    M = bytes(M, "utf-8")
    print("The message decrypted into bytes is as follows: ", M)
    v = compute_F_h(M,V)
    print("v: ", v)
    print("-------------------------------------")
    print("STEP 3: compute second signiture element using first signiture element, v, random integer x, and random integer k")
    s = compute_s(k, x, v, q)
    print("s: ", s)    
    
    #####################################################
    # SIGNITURE VERIFICATION ALGOIRTHM
    #####################################################
    print("-------------------------------------")
    print("################################################")
    print("SIGNITURE VERIFICATION ALGORITHM")
    print("################################################")
    print("-------------------------------------")

    # Calculate v' and compare with v to see if v'=v
    V_prime = (y**(-v))*(t)*(z**s) 
    v_prime = compute_F_h(M,V_prime)
    print("v_prime: ", v_prime)

    if(V_prime == v):
        print("1: encrypted message == message")
    else:
        print("0: encrypted message =/ message")




if __name__ == "__main__":
    main()

