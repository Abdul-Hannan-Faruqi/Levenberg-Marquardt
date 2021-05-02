# Levenberg-Marquardt

1. **Least Squares Minimization**
  * Define the problem in `Le_sq.c`, including the number of variables `n` and the number of equations `nf` 
  * Define the number of variables `n` in `Levenberg_Marquardt.c`
  * Build using `./build.sh`

2. **General Optimization problem with single objective**
  * Define the objective function and the number of variables `n` in `fx.c`
  * Define the number of variables `n` in `Levenberg_Marquardt.c`
  * Build using `./bfx.sh`

After building the desired binary, run the execution script using `./exec.sh`. This will execute the pre-written script in `exec.sh` to generate results for different parameter values. The script may be modified according to requirement. The executable for the project is titled `LM` and may be executed directly using

    ./LM arg1 arg2 arg3 arg4
    
It takes four runtime arguments defined as
* arg1 :  choice of line search - 
  * 0 for exact
  * 1 for inexact
* arg2 :  value of gamma
* arg3 :  value of mu_1
* arg4 :  value of mu_2

The results of the run get appended to the file `log.dat`

***
### Data format in `log.dat`   
The following properties for each run are stored in the file, delimited by a space   
Line-search-startegy &nbsp;&nbsp; gamma &nbsp;&nbsp; mu_1 &nbsp;&nbsp; mu_2 &nbsp;&nbsp; Function-value-at-convergence &nbsp;&nbsp; Number-of-iterations &nbsp;&nbsp; Solution-time
