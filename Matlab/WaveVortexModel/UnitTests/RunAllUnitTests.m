% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationXY);
 % diffTest = matlab.unittest.TestSuite.fromClass(?TestOrthogonalSolutionGroups);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestRadialTransformation);
diffTest = matlab.unittest.TestSuite.fromClass(?TestWVTransformInitialization);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationZ);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestNonlinearFlux);
result = run(diffTest);
rt = table(result)

% suite = matlab.unittest.TestSuite.fromFile('TestSpectralDifferentiation.m');
% {suite.Name}'
% suite.run