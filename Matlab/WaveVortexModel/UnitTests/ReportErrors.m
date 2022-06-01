function ReportErrors(f,Df_analytical,Df_numerical,testname,parameter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
failures = 0;
for iParameter = 1:length(parameter)
    x = f(parameter(iParameter));
    Dx_expected = Df_analytical(parameter(iParameter));
    Dx = Df_numerical(x);
    
    if (max(abs(Dx(:)-Dx_expected(:))) > 1e-7)
        failures = failures + 1;
    end
end

if failures == 0
    fprintf('\t%s passed.\n',testname)
else
    fprintf('\t%s failed with %d failures.\n',testname,failures)
end
end

