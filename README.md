# SPPT with earlier versions of WRF

The purpose of this project is to shift stochastic parametrization scheme to earlier versions of WRF (3.2.1 and 3.6.1)
This package can use WRF sppt scheme by parallel and have continuous random pattern in DA-cycled.

The setting for namelist as follows:

&stoch

sppt (max_dom)                 = 0        ; Stochastically perturbed physics tendencies (SPPT), 0: off, 1: on

gridpt_stddev_sppt (max_dom)   = 0.5      ; Standard deviation of random perturbation field at each gridpoint. ; Determines amplitude of random perturbations

lengthscale_sppt (max_dom)     = 150000.0 ; Perturbation lengthscale (in m).

timescale_sppt (max_dom)       = 21600.0  ; Temporal decorrelation of random field (in s).

stddev_cutoff_sppt (max_dom)   = 2.0   ; Cutoff tails of perturbation pattern above this threshold standard deviation.

nens                          =1        ; Seed for random number stream. For ensemble forecasts this parameter needs to be different for each member. The seed is a function of initial start time to ensure different random number streams for forecasts starting from different initial times. Changing this seed changes the random number streams for all activated stochastic parameterization schemes.

iseed_sppt                   =17 ; Seed for random number stream for sppt. Will be combined with seed nens signifying ensemble member number and initial start time to ensure different random number streams for forecasts starting from different initial times and for different ensemble members.

start_from_wrfinput                   =.false.   - "switch to control restart in DA-cycled forecasts"

I have already tested all functions in WRF 3.2.1 to make sure it work. If you have any problem please let me know.
