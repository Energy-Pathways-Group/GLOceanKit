import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageResult

runner = testrunner("textoutput");
% format = CoverageResult();
% plugin = CodeCoveragePlugin.forFile("../FlowComponents/WVGeostrophicMethods.m",Producing=format);
% runner.addPlugin(plugin)

% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationXY);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestSpectralDifferentiationZ);
diffTest = matlab.unittest.TestSuite.fromClass(?TestOrthogonalSolutionGroups);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestRadialTransformation);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestWVTransformInitialization);
% 
% diffTest = matlab.unittest.TestSuite.fromClass(?TestNonlinearFlux);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestInertialOscillationMethods);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestGeostrophicMethods);
% diffTest = matlab.unittest.TestSuite.fromClass(?TestInternalGravityWaveMethods);
result = runner.run(diffTest);

% generateHTMLReport(result);
% generateHTMLReport(format.Result);
rt = table(result)

% suite = matlab.unittest.TestSuite.fromFile('TestSpectralDifferentiation.m');
% {suite.Name}'
% suite.run