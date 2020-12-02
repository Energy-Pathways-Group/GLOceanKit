Linear velocity field parameter estimation
==============

This directory contains tools for estimating the flow field of the advection-diffusion of the linear velocity field model using Lagrangian data. There are two basic techniques 1) fit the second moment ellipse of the data to analyatical solutions and 2) do a least squares fit to a Taylor expansion of the velocity field.

### Table of contents
1. [Overview](#overview)
1. [Least squares fits](#least-squares-fits)
1. [Second moment fits](#second-moment-fits)


------------------------

Overview
-----------

All methods attempt to estimate (at least) strain rate, vorticity, and divergence from a cluster of drifters. The `ModelParameter` class contains a list of these (and other) parameters that can be estimated from the drifters.

The simplest example of how this works is to generate trajectories, then use those trajectories to estimate parameters.

```matlab

% Use these parameters to generate a few synthetic trajectories.
sigma = 8e-6;
theta = 30*pi/180;
zeta = 0;
kappa = 1;

velocityField = LinearVelocityField(sigma,theta,zeta);
integrator = AdvectionDiffusionIntegrator(velocityField,kappa);

% Particles placed in a "plus" configuration.
x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];

% Now compute the trajectories
T = 3*86400;
dt = 4*3600;
[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

% And then try to estimate the parameters 
parametersToEstimate = [ModelParameter.strain];
parameterEstimates = EstimateLinearVelocityFieldParameters( x, y, t, parametersToEstimate );
```

Least squares fits
------------

This methodology was published in Oscroft, Sykulski and Early. 

### EstimateLinearVelocityFieldParameters.m

Given Lagrangian trajectories and a list of parameters to estimate, this functional will return the best estimates possible.

The degrees-of-freedom (dof) can optionally be given in order to set how many degrees-of-freedom each parameter is allowed to have. Anything above 1 will generate a B-spline basis to allow for time variation in the parameter estimates, which will then call `EstimateLinearVelocityFieldParametersBSplines`.

### EstimateLinearVelocityFieldParametersBSplines.m

This is a more primative function than `EstimateLinearVelocityFieldParameters`, where you can manually choose the splines (or any basis function really) that you want to use. 

### DecomposeTrajectories.m

Given a Lagrangian trajectories and a set of parameter estimates, this will decompose the velocity time series of each trajectory into background, mesoscale, and submesoscale parts. 

### EstimateSolutionLikelihoodFromBootstraps.m

Given a struct of bootstrap estimates, this will construct PDFs from the estimates, and use that to score the likelihood of each bootstrap estimate. This can be used to determine the most probably solution.



Second moment fits
------------

The following scripts are all related to using the second moment solutions to estimate strain rate, vorticity, and diffusivity.

### MomentTensorEvolutionInStrainVorticityField.m

Contains the analytical solutions to the evolution of the second moment ellipse. These solutions were in the Latmix paper.

### MomentTensorEvolutionInTimeVaryingStrainVorticityField.m

Same as above, but this time the strain, vorticity and diffusivity parameters are allowed to vary with time. This a silly implemention that simply uses the constant in time solution between time points. Should be fine?

### MomentTensorModelError.m

This script takes an *observed* second moment time series and computes the error compared with an *analytical* second moment solution. You can choose different error metrics, include the one in the Latmix paper. The parameters are scaled exponential to allow for efficient search algorithm.


### FitSecondMomentToEllipseModel.m

Given an observed second moment time series, this uses MomentTensorModelError.m to find the minimum solution, constant in time.

### FitTrajectoriesToEllipseModel.m

Very lightweight script that simply wraps around FitSecondMomentToEllipseModel.m by computing the second moment for you, directly from trajectories.

### FitTrajectoriesToEllipseModelWithJacknife.m

Essentially the same as FitTrajectoriesToEllipseModel.m, but performs the fits using N-1 drifter combinations..

### MomentTensorModelErrorForBSplines.m and FitSecondMomentToTimeVaryingEllipseModel.m

Attempts to the use B-splines to allow for slowy varying fit parameters. While it does work, it doesn't work well. The primary problem seems to be finding the minimum with fminsearch. Overall, I suspect this just isn't a good approach.





```matlab

jet = MeanderingJet();
```




