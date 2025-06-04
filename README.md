# Population-V1MT
 
This package simulate the neural distributed response of V1-MT pathway in response on plaid. The package is organised with scripts and folders.

FILTER folder -> contains .mat of gabor filters needed by the architecture to simulate an MT neuron.
FUNCTIONS folder -> contains needed functions sed by the main scripts
SIMULATIONS -> contains .mat of simulated response of MT neurons with different weighting methods and different normalisation levels
FIGs -> some results

# Main scripts

_MainScript_ used to made the simulation. Simulate a network under different conditions. First section Population Initialisation with param.
param is a structure which contains different field used to characterise the network (from the desidered spatial frequencies to the temporal profile of the model). 
A documentation on the model parameters wil be defined in the future.

Stim is an another structure needed for the decide the stimulus parameters (kind of stimulus etc.)

# SIMULATIONS
Simulations are runned from script motionPopV1MT (see in FUNCTIONs). It is a method that take stim and param structures
