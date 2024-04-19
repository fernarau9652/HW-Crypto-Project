# This is to test the product operations of vectors.
# Perform mathematical operations
import pandas as pd
def math_operations():
    operations = []
    # Perform addition, subtraction, multiplication, and division
    a = [1,2,3,4]
    b = [0,3,4,2]
    
    # create instance of BVMT table
    df = e_table()

    # do operations to append proper operations defined by A o B formula
    base_a = df.columns
    base_b = [0,1,2,3]

   # print(a)
   # print(b)
    result = []

    for i in range(len(a)):
        row_result = []
        for j in range(len(b)):
            product = a[i]*b[j]
            # append operation for basis vectors
            if isinstance((type(df.at[base_b[j], base_a[i]])), str):
                product.append('*', df.at[base_b[j], base_a[i]])
            elif isinstance((type(df.at[base_b[j], base_a[i]])), int):
                product = product * df.at[base_b[j], base_a[i]]
            row_result.append(product)
        result.append(row_result)
    
    # append operations for bsais table
     

    return result


# Create some sample data
def e_table():
    data = {
    'e0': ['ue0', 0, 0, 'ue3'],
    'e1': ['0', 'le1', 'le2', '0'],
    'e2': ['ue2', '0', '0', 'ue1'],
    'e3': ['0', 'le3', 'le0', '0']
}
    # Create a DataFrame
    df = pd.DataFrame(data)
    #print(type(df.at[2, 'e0']))
    return (df)


# Write results to a file
def write_to_file(operations):
    with open("/home/u1150934/cryptography/HW-Crypto-Project/Python/math_operations.txt", "w") as file:
        for operation in operations:
            file.write(operation + "\n")

if __name__ == "__main__":
    # Perform mathematical operations
    operations = math_operations()
    df = e_table()
   # print(df)
    print(operations)
    # Write results to a file
   # write_to_file(operations)








