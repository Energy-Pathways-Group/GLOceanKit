%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CheckWithAnalyticalSolutions
%
% This script validates all states of motion of the wave vortex transform
% against the analytical solutions.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% November 22th, 2023   Version 2.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isHydrostatic = 0;
wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [4, 8, 5],isHydrostatic=isHydrostatic);

solutionGroups{3} = WVMeanDensityAnomalySolutionGroup(wvt);
solutionGroups{4} = WVInertialOscillationSolutionGroup(wvt);
solutionGroups{1} = WVGeostrophicSolutionGroup(wvt);
solutionGroups{2} = WVInternalGravityWaveSolutionGroup(wvt);

for iGroup = 1:length(solutionGroups)
    totalErrors = 0;
    totalTests = 0;
    solnGroup = solutionGroups{iGroup};
    fprintf('\n***************************************************\n');
    fprintf('Testing %s solution group:\n',solnGroup.name);
    for iSoln = 1:solnGroup.nUniqueSolutions
        soln = solnGroup.uniqueSolutionAtIndex(iSoln,amplitude='random');
        wvt.initWithUVEta(soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t), soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t),soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t));
        [totalTests,totalErrors] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln);
    end
    summarizeTestResults(totalErrors,totalTests);
    fprintf('***************************************************\n');
end

function summarizeTestResults(totalErrors,totalTests)
if totalErrors > 0
    fprintf('FAILED %d of %d tests.\n',totalErrors, totalTests);
else
    fprintf('PASSED all %d tests.\n', totalTests);
end
end

function [totalTests,totalErrors] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln)
[u_error,v_error,w_error,eta_error,p_error,coeff_error] = errorsFromSolution(wvt,soln);
max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(eta_error))) round((log10(p_error))) round((log10(coeff_error)))],[],'includenan');

totalTests = totalTests + 1;
if isnan(max_error) || max_error > -10
    totalErrors = totalErrors + 1;
    fprintf('\nFound at large error at (k,l,j)=(%d,%d,%d):\n',soln.kMode,soln.lMode,soln.jMode);
    fprintf('The model solution for (u,v,w,eta,p,A) matches the analytical solution to 1 part in (10^%d, 10^%d, 10^%d, 10^%d, 10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))), round((log10(w_error))), round((log10(eta_error))), round((log10(p_error))), round((log10(coeff_error))),wvt.t);
end
end

function [u_error,v_error,w_error,eta_error,p_error,coeff_error] = errorsFromSolution(wvt,soln)
u_unit = soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t);
v_unit = soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t);
w_unit = soln.w(wvt.X,wvt.Y,wvt.Z,wvt.t);
eta_unit = soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t);
p_unit = soln.p(wvt.X,wvt.Y,wvt.Z,wvt.t);
u_error = max(max(abs(wvt.u(:)-u_unit(:))/wvt.uMax),1e-16);
v_error = max(max(abs(wvt.v(:)-v_unit(:))/wvt.uMax),1e-16);
w_error = max(max(abs(wvt.w(:)-w_unit(:))),1e-16);
eta_error = max(max(abs(wvt.eta(:)-eta_unit(:))/max(abs(eta_unit(:)))),1e-16);
p_error = max(max(abs(wvt.p(:)-p_unit(:))/max(abs(p_unit(:)))),1e-16);
coeff_error = coefficientErrorFromSolution(wvt,soln);
end

function coeff_error = coefficientErrorFromSolution(wvt,soln)
    switch(soln.coefficientMatrix)
        case WVCoefficientMatrix.Ap
            primaryError = abs(wvt.Ap(soln.coefficientMatrixIndex)-soln.coefficientMatrixAmplitude)/abs(soln.coefficientMatrixAmplitude);
        case WVCoefficientMatrix.Am
            primaryError = abs(wvt.Am(soln.coefficientMatrixIndex)-soln.coefficientMatrixAmplitude)/abs(soln.coefficientMatrixAmplitude);
        case WVCoefficientMatrix.A0
            primaryError = abs(wvt.A0(soln.coefficientMatrixIndex)-soln.coefficientMatrixAmplitude)/abs(soln.coefficientMatrixAmplitude);
    end
    if ~isempty(soln.conjugateCoefficientMatrix)
        switch(soln.conjugateCoefficientMatrix)
            case WVCoefficientMatrix.Ap
                conjugateError = abs(wvt.Ap(soln.conjugateCoefficientMatrixIndex)-soln.conjugateCoefficientMatrixAmplitude)/abs(soln.conjugateCoefficientMatrixAmplitude);
            case WVCoefficientMatrix.Am
                conjugateError = abs(wvt.Am(soln.conjugateCoefficientMatrixIndex)-soln.conjugateCoefficientMatrixAmplitude)/abs(soln.conjugateCoefficientMatrixAmplitude);
            case WVCoefficientMatrix.A0
                conjugateError = abs(wvt.A0(soln.conjugateCoefficientMatrixIndex)-soln.conjugateCoefficientMatrixAmplitude)/abs(soln.conjugateCoefficientMatrixAmplitude);
        end
    else
        conjugateError = 1e-16;
    end
    coeff_error = max([primaryError conjugateError]);
end