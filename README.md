# Verification of Post-Quantum Signature Algorithms using Quaternion Algebra
ECE 6960-010: Hardware Security/Cryptography - Final Project \
Institution: University of Utah \
Instructor: Dr. Priyank Kalla \
Contributors: Elmir Dzaka, Fernando Araujo

## Documentation
Documents developed regarding this project can be found in the following path:
```bash
cd ../HW-Crypto-Project/Documentation
```
On this path, you will find the final paper in the form of a PDF:
```bash
./ECE6960-Final_Project.pdf
```
Additionally, you will discover a folder containing notes used to describe certain aspects of the project: \
(NOTE: these notes are not well organized, but they've helped us understand the beginning parts of this project)
```bash
./Notes/Project-Notes.pdf
```

## Installing Needed Libraries
To run the code properly, a Quaternion package must be included to allow for proper non-commutative algebras. For this project, we use **pyquaternion**. To install pyquaternion, simply do:
```bash
pip install pyquaternion
```
## Running the Code
To get the desired output, we use the **quaternion.py** file. To get to the proper location, run the following command:
```bash
cd ../HW-Crypto-Project/Python/
```
From here, simply run the Python file to get the desired output onto the terminal. This can be done with the following command:
```bash
python3 quaternion.py
```
Once run, the terminal will display all the steps and results that were implemented for the project.
 
## Files Archive
The project also has a few more files that were used in development but not incorporated into the final product. This archive will go over all these files and their contribution to the project's development.

### Singular
When starting the project, we first attempted to implement the project inside a Singular file named **project.sing**. The issue we ran into is that although possible, doing non-commutative algebra inside Singular is very complex. Given the complexity, implementing non-commutativity inside the singular file would take a while, and it would be outside the scope of the given project. The file can be viewed given the following command:
```bash
emacs ../HW-Crypto-Project/Singular/project.sing
```   
### Python
During development, we also created some Python files where we were experimenting with non-commutative algebra. The main file we tested logic in is named **product_test.py**. The goal of these files was to get a better understanding of non-commutative algebras and how non-commutativity can be applied in a crypto scheme to create post-quantum security. The file can be viewed given the following command:
```bash
emacs ../HW-Crypto-Project/Python/product_test.py
``` 
