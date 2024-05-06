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

isHydrostatic = 1;
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [4, 8, 5],isHydrostatic=isHydrostatic);
% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [4, 8, 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [4, 8, 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

solutionGroups{1} = WVGeostrophicComponent(wvt);
solutionGroups{2} = WVInternalGravityWaveComponent(wvt);
solutionGroups{3} = WVMeanDensityAnomalyComponent(wvt);
solutionGroups{4} = WVInertialOscillationComponent(wvt);

for iGroup = 1:length(solutionGroups)
    totalErrors = 0;
    totalTests = 0;
    log_max_spatial_error = -Inf;
    log_max_spectral_error = -Inf;

    solnGroup = solutionGroups{iGroup};
    fprintf('\n***************************************************\n');
    fprintf('Testing %s solution group:\n',solnGroup.name);
    for iSoln = 1:solnGroup.nModes
        soln = solnGroup.solutionForModeAtIndex(iSoln,amplitude='random');
        wvt.initWithUVEta(soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t), soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t),soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t));
        [totalTests,totalErrors,log_spatial_error,log_spectral_error] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln);
        if log_spatial_error > log_max_spatial_error
            log_max_spatial_error = log_spatial_error;
        end
        if log_spectral_error > log_max_spectral_error
            log_max_spectral_error = log_spectral_error;
        end
    end
    summarizeTestResults(totalErrors,totalTests,log_max_spatial_error,log_max_spectral_error);
    fprintf('***************************************************\n');
end

function summarizeTestResults(totalErrors,totalTests,max_spatial_error,max_spectral_error)
if totalErrors > 0
    fprintf('FAILED %d of %d tests.\n',totalErrors, totalTests);
else
    fprintf('PASSED all %d tests. Maximum (spatial,spectral) error is 1 part in (10^%.1f, 10^%.1f)\n', totalTests,max_spatial_error,max_spectral_error);
end
end

function [totalTests,totalErrors,log_max_spatial_error,log_max_spectral_error] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln,options)
arguments
    totalTests (1,1) double {mustBeNonnegative}
    totalErrors (1,1) double {mustBeNonnegative}
    wvt WVTransform {mustBeNonempty}
    soln WVOrthogonalSolution {mustBeNonempty}
    options.spatialErrorTolerance (1,1) double = -9
    options.spectralErrorTolerance (1,1) double = -3
end
[u_error,v_error,w_error,eta_error,p_error,coeff_error] = errorsFromSolution(wvt,soln);
log_max_spatial_error = max([((log10(u_error)))  ((log10(v_error))) ((log10(w_error))) ((log10(eta_error))) ((log10(p_error)))],[],'includenan');
log_max_spectral_error = ((log10(coeff_error)));

totalTests = totalTests + 1;
if isnan(log_max_spatial_error) || log_max_spatial_error > options.spatialErrorTolerance || log_max_spectral_error > options.spectralErrorTolerance
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
uMax = max(max(sqrt(u_unit(:).^2 + v_unit(:).^2)),1);
etaMax = max(max(abs(eta_unit(:))),1);
pMax = max(max(abs(p_unit(:))),1);
u_error = max(max(abs(wvt.u(:)-u_unit(:))/uMax),1e-16,"includemissing");
v_error = max(max(abs(wvt.v(:)-v_unit(:))/uMax),1e-16,"includemissing");
w_error = max(max(abs(wvt.w(:)-w_unit(:))),1e-16,"includemissing");
eta_error = max(max(abs(wvt.eta(:)-eta_unit(:))/etaMax),1e-16,"includemissing");
p_error = max(max(abs(wvt.p(:)-p_unit(:))/pMax),1e-16,"includemissing");
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