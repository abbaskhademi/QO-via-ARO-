# QO_Via_ARO
Data and codes of the paper _"Quadratic Optimization Through the Lens of Adjustable Robust Optimization"_ by _A. Khademi and A. Marandi (2024)_.
# Introduction
This repository provides the following:

-The exact problem data that is used in the numerical experiments of the paper

# Dependencies
- MATLAB - We use MATLAB (version  2022a) to run the scripts. 

- YALMIP - We use [YALMIP](https://yalmip.github.io/) (through MATLAB) to call optimization solvers. 

- MOSEK - We use [MOSEK](https://www.mosek.com/) to solve subproblems as a (cone) optimization solver. 

- CPLEX - We use [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) to solve nonconvex quadratic problems as a global solver.

- GUROBI -  We use [GUROBI](https://www.gurobi.com/) to solve nonconvex quadratic problems as a global solver.

  


# Description
The following is a guide for using this repository. Please be aware that all data files are in the ".mat" format, which is commonly used by MATLAB. The scripts are provided in ".m" format, which is the standard file type for MATLAB script files.  



# End note
We are grateful for your interest in our work. Please feel free to reach out to Abbas Khademi (mailto:abbaskhademi@ut.ac.ir ) if you encounter any issues while using the scripts. If you have any other comments or questions, do not hesitate to contact us:

A. Khademi: mailto:abbaskhademi@ut.ac.ir 

A. Marandi: mailto:a.marandi@tue.nl 

If you found this work useful in your research and/or applications, please star this repository and cite:

```bibtex
@article{AROQO2024,
  title={Quadratic Optimization Through the Lens of Adjustable Robust Optimization},
  author={Khademi, Abbas and Marandi, Ahmadreza},
  journal={Optimization Online},
  year={2024}
}

