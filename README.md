# IPsolver
Optimization using symbolic variables in MATLAB

## About This Project
* This solver is built based on interior point method for convex problems. It accepts non-linear constraints.<br />
* The basic idea of the solver is to calculate the Newton step offline, which is used every time we call the solver.<br />
* The symbolic variables are used to define the follows:<br />
1- The decision variables that we want to minimize.<br />
2- Parameters that are constants with respect to the optimization problem. We feed the exact values of the parameters when we call the solver.
 An example for parameters utilization, let us consider we want to minimize an area of a triangle and there is the minimum length for each side that may varies based on, for example, the application. So we can consider the minimum length as parameters. In example 2 and 3, parameters are clarified.

 
## Prerequisites
MATLAB (with symbolic variables' toolbox).<br />
This project is built in MATLAB 2019a.

## Installing
1- Download the repository on your desktop.<br />
2- Run the file install.m in MATLAB.


## Using the solver
* The main function is (IPsovler.m); after installation type the command (help IPsolver) for more details.<br />
* The solver provides options that contribute to improving the performance and reduce the running time, after installation type the command (options_check) for more details.<br />
* The internal parameters of the solver are defined in function (solver_settings) located at "repository_root/functions/Auxiliary_functions" It is possible to tune the parameters to customize the performance for specific applications.
* Example are provided to give an idea of how to use the solver.



## Authors

* **Anas Abdelkarim**<br />
   For inquiries and suggestions please contact me on eng.anas_arouri@yahoo.com

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.<br />
Please cite this repository if you find it helpful.


