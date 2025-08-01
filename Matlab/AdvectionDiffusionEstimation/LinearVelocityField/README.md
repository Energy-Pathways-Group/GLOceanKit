Linear velocity field parameter estimation
==============

This directory contains tools for estimating the flow field of the advection-diffusion of the linear velocity field model using Lagrangian data. There are two basic techniques 1) do a least squares fit to a Taylor expansion of the velocity field and, 2)  fit the second moment ellipse of the data to analyatical solutions.

### Table of contents
1. [Overview](#overview)
1. [Least squares fits](#least-squares-fits)
1. [Second moment fits](#second-moment-fits)


------------------------

Overview
-----------

Both estimation methods estimate linear velocity field parameters such as  strain rate, vorticity, and divergence from a cluster of drifters. The `ModelParameter` class contains a list of these (and other) parameters that can be estimated from the drifters.

The simplest example of how this works is to either a) generate trajectories using the [advection-diffusion models](../../AdvectionDiffusionModels), or b) use your own trajectories and then use those trajectories to estimate parameters.

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

The structure `parameterEstimates` now contains values for `sigma_n` and `sigma_s` which, hopefully, give you back something close to what you put in.

Recommended approach for least squares fits
------------

To reproduce the results from the Latmix Site 1 drifter fits, try the following approach using the code in the folder `AdvectionDiffusionEstimation/LinearVelocityField/FluidsPaperFigures/`. The data is found in `AdvectionDiffusionEstimation/Fluid Paper Code`.

1. Run `GenerateBootstrapFits.m`, which will read the `smoothedGriddedRho1Drifters.mat` file and then perform estimation of all the models, with 1000 different drifter permutations, with time variation from 1 to 6 degrees-of-freedom. The data will be stored in 6 different `.mat` files in the BootstrapData folder. Note that the velocity is computed from the positions with a second-order finite difference matrix, i.e., line 80 in `EstimateLinearVelocityFieldParametersBSpline.m` calls `D = FiniteDifferenceMatrix(1,t,1,1,2);`.

2. `CompareBootstrapFits.m` will produce a Latex table that will let you assess the quality of fit, and choose the best model. Refer to [section 5.1 in the manuscript](https://www.mdpi.com/2311-5521/6/1/14).

3. The script `MakeFigureBestSplineFitSite1` will produce a plot showing the time varying parameters of any of the resulting models, including the error bars. We chose the model with 4 dof.  This reproduces [figure 13 in the manuscript](https://www.mdpi.com/2311-5521/6/1/14).

4. To actually use the decomposition, you will want to find the most probable set of parameters and then apply the decomposition,
```matlab
[~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
...
[u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(x, y, t, parameterEstimates);
```
which is done in, e.g., the script `MakeFigureBestSplineFitSite1.m` or `Website/MakeSite1DecompositionMovie`.



Least squares fits
------------

This methodology was published in Oscroft, Sykulski and Early.

### [EstimateLinearVelocityFieldParameters.m](EstimateLinearVelocityFieldParameters.m)

This function estimates linear velocity field parameters strain, vorticity and divergence from a cluster of Lagrangian particles.

The degrees-of-freedom (dof) can optionally be given in order to set how many degrees-of-freedom each parameter is allowed to have. Anything above 1 will generate a B-spline basis to allow for time variation in the parameter estimates, which will then call `EstimateLinearVelocityFieldParametersBSplines`.

The unit test [EstimateLinearVelocityFieldParametersUnitTest](UnitTests/EstimateLinearVelocityFieldParametersUnitTest.m) demonstrates the estimation procedure on a range of different parameters.

### [EstimateLinearVelocityFieldParametersBSplines.m](EstimateLinearVelocityFieldParametersBSplines.m)

This is a more primative function than `EstimateLinearVelocityFieldParameters`, where you can manually choose the splines (or any basis function really) that you want to use for the time variation in the parameter estimates.

### [DecomposeTrajectories.m](DecomposeTrajectories.m)

Given a Lagrangian trajectories and a set of parameter estimates, this will decompose the velocity time series of each trajectory into background, mesoscale, and submesoscale parts, following the approach in the manuscript.

The unit test [EstimateLinearVelocityFieldParametersUnitTest](UnitTests/EstimateLinearVelocityFieldParametersUnitTest.m) demonstrates how estimated parameters can be used to decompose the signal, and then estimate submesoscale diffusivity.

### [EstimateSolutionLikelihoodFromBootstraps.m](EstimateSolutionLikelihoodFromBootstraps.m)

Given a struct of bootstrap estimates, this will construct PDFs from the estimates, and use that to score the likelihood of each bootstrap estimate. This can be used to determine the most probable solution.



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




