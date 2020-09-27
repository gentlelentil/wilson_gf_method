Author - Nathaniel Smith

Program accepts atoms via a Cartesian coordinate system, and calculates the B matrix from said coordinates.

B Matrix is the initial step of the GF method developed by Wilson, Decius, and Cross in their book 'Molecular Vibrations' 

Calculates a coefficient matrix by finding the S Vectors of each system.
The S Vectors are better described in 'Molecular Vibrations' pp56.

example input: [[1,2,3],[4,5,6],[7,8,9]]
Representing a 3 atom system.
S Vectors = 3N-6 = 3 S vectors

Returns 3N-6 X 3N matrix.
 
