classdef TestNetCDF < matlab.unittest.TestCase
    properties
        ncfile
        path = "test.nc"
        x = 0:9
        f_a = @(x) 4*x + sqrt(-1)*2*x;
        b = logical([0 1 0 1 1 0 1 0 1 1]);
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            if isfile(testCase.path)
                delete(testCase.path);
            end
            testCase.ncfile = NetCDFFile(testCase.path,shouldOverwriteExisting=1);
        end
    end

    methods (TestClassTeardown)
        function classTeardown(testCase)
            testCase.ncfile.close();
        end
    end

    methods (Test)
        function testAddDimension(testCase)
            testCase.verifyWarningFree(@() testCase.ncfile.addDimension('x',testCase.x) );
        end

        function testReadDimension(testCase)
            x_back = testCase.ncfile.readVariables('x');
            testCase.verifyEqual(x_back,testCase.x);
        end

        function testAddComplexVariable(testCase)
            testCase.verifyWarningFree(@() testCase.ncfile.addVariable('a',{'x'},testCase.f_a(testCase.x)) );
        end

        function testReadComplexVariable(testCase)
            a_back = testCase.ncfile.readVariables('a');
            testCase.verifyEqual(a_back,testCase.f_a(testCase.x));
        end

        function testAddLogicalVariable(testCase)
            testCase.verifyWarningFree(@() testCase.ncfile.addVariable('b',{'x'},testCase.b) );
        end

        function testReadLogicalVariable(testCase)
            a_back = testCase.ncfile.readVariables('b');
            testCase.verifyEqual(a_back,testCase.b);
        end

        function testAddColumnVector(testCase)
            testCase.verifyWarningFree(@() testCase.ncfile.addVariable('y',{'x'},reshape(testCase.x,[],1)) );
        end

        function testReadColumnVector(testCase)
            y = reshape(testCase.x,[],1);
            y_back = testCase.ncfile.readVariables('y');
            testCase.verifyEqual(y_back,y);
            testCase.verifyNotEqual(y_back,testCase.x);
        end

        function testUnlimitedDimension(testCase)
            s = [1 2.5 3.14].';
            [~,var] = testCase.ncfile.addDimension('t',length=Inf,type='double');
            for i=1:length(s)
                var.setValueAlongDimensionAtIndex(s(i),{'t'},i);
            end

            % Now test a bunch of different ways of reading the data back
            for i=1:length(s)
                val = var.valueAlongDimensionAtIndex({'t'},i);
                testCase.verifyEqual(val,s(i));
            end
            for i=1:length(s)
                val = testCase.ncfile.readVariablesAtIndexAlongDimension({'t'},i,'t');
                testCase.verifyEqual(val,s(i));
            end
            testCase.verifyEqual(var.value,s);
        end

        function testAddGroup(testCase)
            grp = testCase.ncfile.addGroup("MyGroup");
            grp.addAttribute('MyAttribute',"Hello group!")
            s = 0:9;
            grp.addDimension('s',s);
            func = 4*s + sqrt(-1)*2*s;
            grp.addVariable('AFunction',{'s'},func);

            func_back = testCase.ncfile.readVariables('MyGroup/AFunction');
            testCase.verifyEqual(func_back,func);
        end
    end

end