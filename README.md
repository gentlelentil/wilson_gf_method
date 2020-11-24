Author - Nathaniel Smith

Program accepts atoms via a Cartesian coordinate system, and calculates the B matrix from said coordinates.

B Matrix is the initial step of the GF method developed by Wilson, Decius, and Cross in their book 'Molecular Vibrations' pp54-60

Calculates a coefficient matrix by finding the S Vectors of each system.
The S Vectors are better described in 'Molecular Vibrations' pp56.

An example xyz file is present, load with numpy.loadtxt()

Representing a 3 atom system.
S Vectors = 3N-6 = 3 S vectors

Returns 3N-6 X 3N matrix.
 
