
from pyquaternion import Quaternion
import random


def modulo(quaternion, q):
    # initalize array for new quaternion computation
    new_quaternion_arr = [];
    for comp in quaternion:
        new_quaternion_arr.append(comp%q);

    # create new quaternion from the arry and return
    new_quaternion = Quaternion(new_quaternion_arr[0], new_quaternion_arr[1], new_quaternion_arr[2], new_quaternion_arr[3])
    
    #return new_quaternion
    return new_quaternion

def vector_product_verify(Q, G, N):
    if (Q*N != N*Q) and (G*N != N*G) and (Q*G != G*Q):
        return True
    else:
        return False

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

    return new_quaternion

def compute_y(A, B, N, Ln, q):
    Y_arr = []
    x = random.randint(0,q-1)
    #x = 1

    # raise N to the power of x
    N_raised = N ** x
    sub1 = A*N_raised
    sub2 = sub1*Ln
    Y = sub2*A.inverse

    return Y

def compute_z(A, B, N, Rn):
    sub1 = B*Rn
    sub2 = sub1 * N
    Z = sub2 * B.inverse
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
    E_double = N ** (q-1)
    return E_double

def compute_T(A,B,E_double):
    sub1 = B * E_double
    T = sub1 * A.inverse
    return T

def main ():
    #####################################
    # KEY GENERATION ALGORITHM 
    #####################################

    # Creating a random Quaternions
    A = Quaternion(1,2,3,4)
    B = Quaternion(2,3,4,5)

    # Non-Invertable Quaternion N
    N = Quaternion(2,1,1,2)

    # Print Non-commucative mutiplication
    print(A*B)
    AB = A*B;
    
    # do moudlus q operation
    q = 7;
    
    # A % q
    A_mod = modulo(A,q);
    print(A_mod) 

    # A*B (mod q)
    AB_mod = modulo(AB,q)
    print(AB_mod)

    # Step 1
    # Verify vectors
    print(vector_product_verify(A,B,N))

    # step two: Left sided unit
    ### d = 0
    d = 0
    Ln_0 = left_sided_unit(A,q,d)
    print(Ln_0)

    # Right Sided Unit
    Rn_0 = right_sided_unit(A,q,d)
    print(Rn_0)

    ####################################
    # d = 1
    d = 1
    Ln_1 = left_sided_unit(A,q,d)
    print(Ln_1)

    # maybe do modulo? 
    #print(modulo(Ln_1,q))

    # right-sided unit 
    Rn_1 = right_sided_unit(A,q,d)
    print(Rn_1)

    ####################################
    # d = 2
    d = 2
    Ln_2 = left_sided_unit(A,q,d)
    print(Ln_2)

    # Right Sided Unit
    Rn_2 = right_sided_unit(A,q,d)
    print(Rn_2)

    ####################################
    # d = 3
    d = 3
    Ln_3 = left_sided_unit(A,q,d)
    print(Ln_3)

    # right-sided unit 
    Rn_3 = right_sided_unit(A,q,d)
    print(Rn_3)

    ####################################
    # d = 4
    d = 4
    Ln_4 = left_sided_unit(A,q,d)
    print(Ln_4)

    # Right Sided Unit
    Rn_4 = right_sided_unit(A,q,d)
    print(Rn_4)

    ####################################
    # d = 5
    d = 5
    Ln_5 = left_sided_unit(A,q,d)
    print(Ln_5)

    # right-sided unit 
    Rn_5 = right_sided_unit(A,q,d)
    print(Rn_5)

    ####################################
    # d = 6
    d = 6
    Ln_6 = left_sided_unit(A,q,d)
    print(Ln_6)

    # Right Sided Unit
    Rn_6 = right_sided_unit(A,q,d)
    print(Rn_6)

    ####################################

    # TODO: Figure out what p is and pass into functions
    # step 3: compute Y
    Y = compute_y(A,B,N,Ln_0,q)
    print(Y)

    # Step 4: compute Z
    Z = compute_z(A, B, N, Rn_0)
    print(Z)

    # Step 5: compute E''
    # compute forumula 9
    E_double = compute_e_double(N,q)
    print(E_double)

    # compute formula 10 to compare with 9
    E_double_w = compute_e_double_w(N,q)
    print(E_double_w)

    # Step 6: Compute T
    T = compute_T(A, B, E_double)
    print(T)

    #####################################################
    # SIGNITURE VERIFICATION ALGOIRTHM
    #####################################################


if __name__ == "__main__":
    main()

