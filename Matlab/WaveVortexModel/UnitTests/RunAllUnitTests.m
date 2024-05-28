% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationXY);
diffTest = matlab.unittest.TestSuite.fromClass(?TestOrthogonalSolutionGroups);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestRadialTransformation);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestWVTransformInitialization);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationZ);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestNonlinearFlux);
 % diffTest = matlab.unittest.TestSuite.fromClass(?TestInertialOscillationMethods);
 % diffTest = matlab.unittest.TestSuite.fromClass(?TestGeostrophicMethods);
result = run(diffTest);
rt = table(result)

% suite = matlab.unittest.TestSuite.fromFile('TestSpectralDifferentiation.m');
% {suite.Name}'
% suite.run