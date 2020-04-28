# IP_COVID-19
**Interval prediction code for COVID-19**

This is the code accompanying the report: https://hal.inria.fr/hal-02517866 (versions > 5)

The script to run is `run_main1.m`, which contains all parameters and the data for all countries.
The variables I0, D0 and H0 correspond to \mathcal{I}, \mathcal{D} and \mathcal{R} in the report, respectively.
The description of other parameters can be found in the file.

The function `script1.m` contains all identification and validation procedures for SEIR model with 4 compartments: *S*, *E*, *I* and *R* (the compartment *R* includes dead and recovered individuals).
See the comments in the file. 

The function `script2.m` contains all identification and validation procedures for SEIR model with 5 compartments: *S*, *E*, *I*, *D* and *R* (dying are considered separately in *D*).

The function `delay_ind.m` is used to calculate the delayed values of the input in the model and predictor.
