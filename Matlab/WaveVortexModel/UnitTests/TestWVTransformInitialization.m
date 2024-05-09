classdef TestWVTransformInitialization < matlab.unittest.TestCase

    properties (TestParameter)
        latitude = {0,5,10,90}
    end

    methods (Test)
        function testInitWithLatitude(testCase,latitude)
            if latitude < 5
            testCase.verifyError(@() WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5],latitude=latitude),'MATLAB:validators:mustBeGreaterThanOrEqual' );
            elseif latitude > 85
                testCase.verifyError(@() WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5],latitude=latitude),'MATLAB:validators:mustBeLessThanOrEqual' );
            else
                testCase.verifyWarningFree(@() WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5],latitude=latitude) );
            end
        end


    end

end