
# Bayesian microsaccade detection
#### Implementation of the algorithm described in [Mihali, van Opheusden, and Ma, 2017, Journal of Vision](http://jov.arvojournals.org/article.aspx?articleid=2597868)
## Installation

* To run this, you will have to make sure you have a C compiler installed on your computer. If not, install it.

  * MacOS: 
Install Xcode Command Line Tools by entering in Terminal: `xcode-select --install`
  * Windows: 
Install [Windows SDK](https://developer.microsoft.com/en-us/windows/downloads/windows-10-sdk) or [MinGW](http://www.mingw.org/wiki/howto_install_the_mingw_gcc_compiler_suite)

* Install the [boost](http://www.boost.org/) C++ package.

* Compile the `bmd.cpp` code into the `bmd` executable in the Command Prompt.
  * `g++  -std=c++11 -O3 -I$BOOST_INC  -fexpensive-optimizations -Wall -Wextra -o bmd bmd.cpp`
 
## Script usage
* Preprocess raw eye trace `x1.mat` with `preprocess_data.m` to match the isotropy assumption in our generative model (see paper). This writes the output to `x1.txt`.

* Run BMD inference algorithm.
  * `./bmd x1.txt integral_table.txt params.txt changepoints.txt >output.txt`

* Analyze and plot the output. 

    * Run the Matlab script `BMD_results_analysis_and_visualization.m`, which reads in the change points from the `changepoints.txt` file, converts them to the eye state time series C with the helper function `tau_to_C.m` and averages across these C samples to give the probability of the eye being in a microsaccade state across the whole eye position time series `x1.txt`. This can be visualized in `microsaccades_inferred_BMD.pdf`. 

## Output files

* `changepoints.txt`
  * `N` total number of 1’s 
  * `n` total number of change points
  * `t01` indices in the time series for 0 (drift) to 1 (microsaccade) eye state transitions
  * `t10` indices in the time series for 1 to 0 eye state transitions
 
* `params.txt` has the following parameters (columns). Each row represents the values across 6 iterations of the BMD algorithm. 
 * `sigmaz`
 * `sigmax`
 * `d0`
 * `sigma0`
 * `d1`
 * `sigma1`
 
The user might need `params.txt` when running the BMD algorithm on simulated data and checking parameter recovery. Also, the BMD variant BMD reduced plus threshold takes as input the inferred values for `sigmaz` and `sigmax`.

* `output.txt`: logposterior values across iterations and estimated parameter values.

##  BMD variants and additional scripts 

* BMD with parallel tempering: `BMD_pt.cpp` with `BMD_pt.h`
* BMD reduced plus threshold: function `bmd_reduced_thresh.m`. To visualize the output, set the appropriate flags in the `BMD_results_analysis_and_visualization.m` script
* simulated data generation script according to our generative model with several values of motor and measurement noise: `sim_data_create.m` and `sim_data_write.m`
