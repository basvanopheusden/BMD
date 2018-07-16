
# Bayesian microsaccade detection
#### Implementation of the algorithm described in [Mihali, van Opheusden, and Ma, 2017, Journal of Vision](http://jov.arvojournals.org/article.aspx?articleid=2597868)

## Applicability
The BMD algorithm works best on fixational eye movement time series with inferred measurement noise of 0.01 deg or higher. Chances are high data that has been filtered online by the eye tracker will have substantially lower noise. With low noise, the discrepancy between the assumptions of the generative model and data tends to cause mistakes in the inference such as detecting multiple small microsaccades instead of one (Figure A7).

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
* Preprocess raw eye traces `xi.mat` with `preprocess_data.m` to match the isotropy assumption in our generative model (see paper). This writes the output to `xi.txt`.

* Run BMD inference algorithm, either on the command line or on batch with `./run_BMD_all.sh`
  * `./bmd x1.txt integral_table.txt params1.txt changepoints1.txt >output1.txt`

* Analyze and plot the output. 

    * Run the Matlab script `BMD_vis.m`, which reads in the change points from the `changepointsi.txt` file, converts them to the eye state time series C with the helper function `tau_to_C.m` and averages across these C samples to give the probability of the eye being in a microsaccade state across the whole eye position time series `xi.txt`. This can be visualized in `microsaccades_inferred_BMDi.pdf`. 

## Output files

* `changepoints.txt`
  * `N` total number of 1’s 
  * `n` total number of change points
  * `t01` indices in the time series for 0 (drift) to 1 (microsaccade) eye state transitions
  * `t10` indices in the time series for 1 to 0 eye state transitions
 
* `params.txt` has the following parameters (columns). Each row represents the values across N iterations of the BMD algorithm. 
 * `sigmaz` (deg/ sqrt(ms))
 * `sigmax` (deg)
 * `d0` (unitless, fixed to 1)
 * `sigma0`(deg/ms)
 * `d1` (unitless)
 * `sigma1` (deg/ms)
 
The user might need `params.txt` when running the BMD algorithm on simulated data and checking parameter recovery. Also, the BMD variant BMD reduced plus threshold takes as input the inferred values for `sigmaz` and `sigmax`.

* `output.txt`: logposterior values across iterations and estimated parameter values.

## Additional assumptions

BMD might need additional constraints in order to give reasonable inferences on particular datasets. Instead of allowing `sigma0` to take any value in between 0.0001 and 1 (function `estimate_sigma_down`), it can be constrained to [0.0001, 0.0013], which is plausible since this upper limit corresponds (see our parametrization of the generalized gamma distribution) to a drift velocity distribution with mean 1.5 deg/s (head fixed mean drift velocities are typically even lower). Inference might improve further by also also lowering the upper limit on `sigma1`, for instance to 0.1. The mean of the generalized gamma distribution for microsaccade velocity will be: $ sigma1 \dot 1000 \dot sqrt(2) \dot gamma((d1+2)/2)/gamma(d1+1)/2)$.

##  BMD variants and additional scripts 

* BMD with parallel tempering: `BMD_pt.cpp` with `BMD_pt.h`
* `./bmd_pt x1.txt integral_table.txt params_pt1.txt changepoints_pt1.txt paral_temp_pt1.txt chain_swaps_pt1.txt >output_pt1.txt`
* BMD reduced plus threshold: function `bmd_reduced_thresh.m`. To visualize the output, set the appropriate flags in the `BMD_vis.m` script
* simulated data generation script according to our generative model with several values of motor and measurement noise: `sim_data_create.m` and `sim_data_write.m`
