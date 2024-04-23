
from pyquaternion import Quaternion

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

def left_sided_unit(A, q):
    # Initialize structural coefficients and left handed unit
    l = 1
    u = 1
    arr = []
    
    # Perform Algorithm
    for d in range(q):
        Ln = []
        sub1 = d
        sub2 = (1/l) - ((u*d) / l) 
        sub3 = (A[1]/A[3]) * d
        sub4 = (A[0] * (1-(u*d)))/(l*A[2])
        Ln.append(sub1)
        Ln.append(sub2)
        Ln.append(sub3)
        Ln.append(sub4)
        arr.append(Ln)
    
    # Create Quaternion and return
    #new_quaternion = Quaternion(arr[0], arr[1], arr[2], arr[3])
    return arr


def right_sided_unit(A, q):
    l = 1
    u = 1
    arr = []

    # Perform right hand algorithm
    for d in range(q):
        Ln = []
        sub1 = d
        sub2 = (1/l) - ((u*d) / l) 
        sub3 =  (A[0] * (1 - (u*d)))/(l*A[3])
        sub4 = (A[1]/A[2])*d
        Ln.append(sub1)
        Ln.append(sub2)
        Ln.append(sub3)
        Ln.append(sub4)
        arr.append(Ln)
    
    # Create Quaternion
    #new_quaternion = Quaternion(arr[0], arr[1], arr[2], arr[3])

    return arr

def main ():
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
    Ln = left_sided_unit(A,q)
    print(Ln)

    # Right Sided Unit
    Rn = right_sided_unit(A,q)
    print(Rn)


if __name__ == "__main__":
    main()

