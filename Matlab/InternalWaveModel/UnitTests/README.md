InternalWaveModel Unit Tests
============================

This folder contains a collection of unit tests for the InternalWaveModel class. Some of these tests are more useful or more automated than others.

CheckAllPlaneWaveSolutions.m
------------------------------
This script initializes the model with a plane wave in *each* resolved wavenumber and compares to the analytical solution. Needs to be manually switched to run through the Arbitrary and Constant stratification classes.

InternalWaveModelGMSpectrumUnitTest.m
-------------------------------
Runs the internal wave model in time with constant stratification, and compares the frequency content to exactly analytical solutions. Note that blurring and aliasing will contaminate the results. This also computes the horizontally averaged variance quantities as a function of depth.

InternalExternalGMSpectrumUnitTest.m
-------------------------------
This script tests to see if the gridded IW solutions are identical to the external IW solutions with the same amplitude and phases. It does this by initializing a randomized Garrett-Munk spectrum on the gridded field, copying the velocity field, then moving those same waves to the external modes and copying the new velocity field.