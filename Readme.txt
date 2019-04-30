
1. Read the given IEEE paper thoroughly to gain understanding of the concept and the code.  Fig 5 is implemented in the attached code.

  N. Joachimowicz, C. Pichot and J. P. Hugonin, "Inverse scattering: an iterative numerical method for electromagnetic imaging," in IEEE Transactions on Antennas and Propagation, vol. 39, no. 12, pp. 1742-1753, Dec. 1991.


2. Load all the files( Escattmeasured.mat, inverse.m, object.mat, objectguess.mat) in a folder.

3. In the command window type:
   load objectguess.mat
   Reconstructed= objectguess;

4. Run the inverse.m code

5. check the results.

Tips:

a) You can change the number of iterations in inverse.m code and check the convergence or reconstruction results

b) Change the value of ''beta'' ( regualrization constant) to see its effect on convergence.


