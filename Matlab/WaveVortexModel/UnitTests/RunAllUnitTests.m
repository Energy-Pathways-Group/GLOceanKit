% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationXY);
diffTest = matlab.unittest.TestSuite.fromClass(?TestOrthogonalSolutionGroups);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationZ);
result = run(diffTest);
rt = table(result)

% suite = matlab.unittest.TestSuite.fromFile('TestSpectralDifferentiation.m');
% {suite.Name}'
% suite.run