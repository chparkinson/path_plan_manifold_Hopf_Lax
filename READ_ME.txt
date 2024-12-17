The code for the manuscript A SCALABLE METHOD FOR OPTIMAL PATH PLANNING ON
MANIFOLDS VIA A HOPF-LAX TYPE FORMULA

The code which will produce each individual figure is separated into
its own folder. In each folder there are only 2 code files, a driver
and the function titled "HJBSolve" which runs Algorithm 1 from the 
paper. The correspondence is as follows:

ex0 - Figure 2 (table of values)
ex1 - Figure 3 (paths on M(x,y) = a*sin(pi*x)*cos(pi*y))
ex2 - Figure 4 (paths on M(x,y) = 2*exp(-(x^2+y^2));
	        this folder contains one extra function
	        which is only there to help with the 
	        plotting)	
ex3 - Figure 5 (high dimensional / scaling example) 

Most are split into further subexamples for the individual figures. 
To run the code for a particular example, simply download the 
directory and run the driver file. The RNG seeds have been fixed 
in each file to produce the exact image/results from the paper
(in our experience, the results are not particularly sensitive
to the random initialization, so these can be changed without
meaningfully affecting results). 

