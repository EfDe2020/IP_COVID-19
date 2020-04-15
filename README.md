# IP_COVID-19
Interval prediction code for COVID-19

This is the code accompanying the report: https://hal.inria.fr/hal-02517866

The script to run is *"run_main1.m"*, which contains all parameters and the data for all countries.
The variables I0, D0 and H0 correspond to \mathcal{I}, \mathcal{D} and \mathcal{H} in the report, respectively.
The description of other parameters can be found in the file.

The function *"script1.m"* contains all identification and validation procedures.
See the comments in the file. 

The function *"script2.m"* contains all identification and validation procedures for a modified model, which takes into account the community feedback on the number of infectives during confinement.

The function *"delay_ind.m"* is used to calculate the delayed values of the input in the model and predictor.
