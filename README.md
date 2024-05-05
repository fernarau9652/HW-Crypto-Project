# HW-Crypto-Project
ECE 6960-010: Hardware Security/Cryptography - Final Project

## Installing Needed Libraries
To be able to run the code properly, a Quaternion package needs to be included to allow for proper non-commutative algebras. For this project, we use **pyquaternion**. To install pyquaternion, simply do:
```bash
pip install pyquaternion
```
## Running the Code
To get the desired output, we use the **quaternion.py** file. To get to the proper location run the following command:
```bash
cd ../HW-Crypto-Project/Python/
```
From here, simply run the python file to get the desired output onto the terminal. This can be done with the following command:
```bash
python3 quaternion.py
```
Once ran, the terminal will display all the steps and results that were implemented for the project.
 
## Files Archive
The project also has a few more files that were used in development, but not incorperated with the final product. This archive will go over all these files, and their contribution to the development of the project.

### Singular
When starting the project, we first attempted to implment the project inside a Singular file named **project.sing**. The issue we ran into is that although possible, doing non-commutative algebra inside Singular is very complex. Given the complexity, implementing non-commutativity inside the singular file would take a while, and it would be outside the scope of the given project. The file can be viewed given the following command:
```bash
emacs ../HW-Crypto-Project/Singular/project.sing
```   
### Python
During development, we also created some python files where we were experimenting with non-commutative algebra. The main file we tested logic in is named **product_test.py** The goal of these files was to get a better understanding of non-commutative algebras, and how non-commutativity can be applied in a crypto scheme to create post-quantum security. The file can be viewed given the following command:
```bash
emacs ../HW-Crypto-Project/Python/product_test.py
``` 
