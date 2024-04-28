function ReportErrors(f,Df_analytical,Df_numerical,testname,parameter)
%ReportErrors Summary of this function goes here
% - Parameter f: function (of parameter omega), to be differentiated
% - Parameter Df_analytical: analytical derivative of function (of parameter omega)
% - Parameter Df_numerical: numerical differential operator that takes a function
% - Parameter testname: name of the test
% - Parameter parameter: array of parameters to feed the functions (omega)
failures = 0;
for iParameter = 1:length(parameter)
    x = f(parameter(iParameter));
    Dx_expected = Df_analytical(parameter(iParameter));
    Dx = Df_numerical(x);
    
    % Check the absolute and relative error
    dDx = Dx(:)-Dx_expected(:);
    if any(abs(dDx)>1e-10 & abs(dDx)./abs(Dx_expected(:)) > 1e-10)
        failures = failures + 1;
    end
end

if failures == 0
    fprintf('\t%s passed.\n',testname)
else
    fprintf('\t%s failed with %d failures.\n',testname,failures)
end
end

