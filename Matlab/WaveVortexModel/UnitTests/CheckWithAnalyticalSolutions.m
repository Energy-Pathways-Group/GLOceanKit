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

%%
totalErrors = 0;
totalTests = 0;
solnGroup = WVGeostrophicSolutionGroup(wvt);
for iSoln = 1:solnGroup.nUniqueSolutions
    soln = solnGroup.uniqueSolutionAtIndex(iSoln,amplitude='random');
    wvt.initWithUVEta(soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t), soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t),soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t));
    [totalTests,totalErrors] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln);
end
summarizeTestResults(totalErrors,totalTests);

%%

[ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = wvt.generateRandomFlowState();
[IO,SGW,IGW,MDA,SG,IG] = wvt.masksForAllFlowConstituents();
RedundantMask = wvt.maskForRedundantHermitianCoefficients();

uniqueGeostrophicModes = IG .* ~RedundantMask;
linearIndices = find(uniqueGeostrophicModes == 1);


totalErrors = 0;
totalTests = 0;

for iMode = 1:length(linearIndices)
    [kIndex,lIndex,jIndex] = ind2sub([wvt.Nk,wvt.Nl,wvt.Nj],linearIndices(iMode));
    wvt.removeAll;
    wvt.A0(linearIndices(iMode)) = A0G(linearIndices(iMode));
    wvt.A0 = wvt.makeHermitian(wvt.A0);
    if kIndex-1 >= wvt.Nx/2
        kMode = (kIndex-1) - wvt.Nx;
    else
        kMode = kIndex-1;
    end
    [u,v,w,eta,p] = wvt.analyticalSolutionA0(kMode, lIndex-1, jIndex-1);

    [u_error,v_error,w_error,eta_error,p_error] = errorsFromAnalyticalSolutions(wvt,u,v,w,eta,p);
    [totalTests,totalErrors] = recordAndReportTestErrors(totalTests,totalErrors, wvt, kIndex, lIndex, jIndex, u_error,v_error,w_error,eta_error,p_error);

    wvt.removeAll;
    wvt.initWithUVEta(u(wvt.X,wvt.Y,wvt.Z,wvt.t), v(wvt.X,wvt.Y,wvt.Z,wvt.t),eta(wvt.X,wvt.Y,wvt.Z,wvt.t));
    [u_error,v_error,w_error,eta_error,p_error] = errorsFromAnalyticalSolutions(wvt,u,v,w,eta,p);
    [totalTests,totalErrors] = recordAndReportTestErrors(totalTests,totalErrors, wvt, kIndex, lIndex, jIndex, u_error,v_error,w_error,eta_error,p_error);
end
summarizeTestResults(totalErrors,totalTests);

%%
uniqueApIGWModes = IGW .* ~RedundantMask;
linearIndices = find(uniqueApIGWModes == 1);


totalErrors = 0;
totalTests = 0;

for iMode = 1:length(linearIndices)
    [kIndex,lIndex,jIndex] = ind2sub([wvt.Nk,wvt.Nl,wvt.Nj],linearIndices(iMode));
    wvt.removeAll;
    wvt.Ap(linearIndices(iMode)) = ApIGW(linearIndices(iMode));
    wvt.Ap = wvt.makeHermitian(wvt.Ap);
    if kIndex-1 >= wvt.Nx/2
        kMode = (kIndex-1) - wvt.Nx;
    else
        kMode = kIndex-1;
    end
    [u,v,w,eta,p] = wvt.analyticalSolutionApm(kMode, lIndex-1, jIndex-1,1);

    [u_error,v_error,w_error,eta_error,p_error] = errorsFromAnalyticalSolutions(wvt,u,v,w,eta,p);
    [totalTests,totalErrors] = recordAndReportTestErrors(totalTests,totalErrors, wvt, kIndex, lIndex, jIndex, u_error,v_error,w_error,eta_error,p_error);

    wvt.initWithUVEta(u(wvt.X,wvt.Y,wvt.Z,wvt.t), v(wvt.X,wvt.Y,wvt.Z,wvt.t),eta(wvt.X,wvt.Y,wvt.Z,wvt.t));
    [u_error,v_error,w_error,eta_error,p_error] = errorsFromAnalyticalSolutions(wvt,u,v,w,eta,p);
    [totalTests,totalErrors] = recordAndReportTestErrors(totalTests,totalErrors, wvt, kIndex, lIndex, jIndex, u_error,v_error,w_error,eta_error,p_error);
end
summarizeTestResults(totalErrors,totalTests);


function summarizeTestResults(totalErrors,totalTests)
fprintf('\n***************************************************\n');
if totalErrors > 0
    fprintf('FAILED %d of %d tests.\n',totalErrors, totalTests);
else
    fprintf('PASSED all %d tests.\n', totalTests);
end
end

function [totalTests,totalErrors] = recordAndReportTestErrors(totalTests,totalErrors, wvt, kIndex, lIndex, jIndex, u_error,v_error,w_error,eta_error,p_error)
max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(eta_error))) round((log10(p_error)))],[],'includenan');

totalTests = totalTests + 1;
if isnan(max_error) || max_error > -10
    totalErrors = totalErrors + 1;
    fprintf('\nFound at large error at (k,l,j)=(%d,%d,%d):\n',kIndex-1, lIndex-1, jIndex-1);
    fprintf('The model solution for (u,v,w,eta,p) matches the analytical solution to 1 part in (10^%d, 10^%d, 10^%d, 10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))), round((log10(w_error))), round((log10(eta_error))), round((log10(p_error))),wvt.t);
end
end

function [u_error,v_error,w_error,eta_error,p_error] = errorsFromAnalyticalSolutions(wvt,u,v,w,eta,p)
u_unit = u(wvt.X,wvt.Y,wvt.Z,wvt.t);
v_unit = v(wvt.X,wvt.Y,wvt.Z,wvt.t);
w_unit = w(wvt.X,wvt.Y,wvt.Z,wvt.t);
eta_unit = eta(wvt.X,wvt.Y,wvt.Z,wvt.t);
p_unit = p(wvt.X,wvt.Y,wvt.Z,wvt.t);
u_error = max(max(abs(wvt.u(:)-u_unit(:))/wvt.uMax),1e-15);
v_error = max(max(abs(wvt.v(:)-v_unit(:))/wvt.uMax),1e-15);
w_error = max(max(abs(wvt.w(:)-w_unit(:))),1e-15);
eta_error = max(max(abs(wvt.eta(:)-eta_unit(:))/max(abs(eta_unit(:)))),1e-15);
p_error = max(max(abs(wvt.p(:)-p_unit(:))/max(abs(p_unit(:)))),1e-15);
end

function [totalTests,totalErrors] = recordAndReportErrorsFromSolution(totalTests,totalErrors, wvt, soln)
[u_error,v_error,w_error,eta_error,p_error] = errorsFromSolution(wvt,soln);
max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(eta_error))) round((log10(p_error)))],[],'includenan');

totalTests = totalTests + 1;
if isnan(max_error) || max_error > -10
    totalErrors = totalErrors + 1;
    fprintf('\nFound at large error at (k,l,j)=(%d,%d,%d):\n',soln.kMode,soln.lMode,soln.jMode);
    fprintf('The model solution for (u,v,w,eta,p) matches the analytical solution to 1 part in (10^%d, 10^%d, 10^%d, 10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))), round((log10(w_error))), round((log10(eta_error))), round((log10(p_error))),wvt.t);
end
end

function [u_error,v_error,w_error,eta_error,p_error] = errorsFromSolution(wvt,soln)
u_unit = soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t);
v_unit = soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t);
w_unit = soln.w(wvt.X,wvt.Y,wvt.Z,wvt.t);
eta_unit = soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t);
p_unit = soln.p(wvt.X,wvt.Y,wvt.Z,wvt.t);
u_error = max(max(abs(wvt.u(:)-u_unit(:))/wvt.uMax),1e-15);
v_error = max(max(abs(wvt.v(:)-v_unit(:))/wvt.uMax),1e-15);
w_error = max(max(abs(wvt.w(:)-w_unit(:))),1e-15);
eta_error = max(max(abs(wvt.eta(:)-eta_unit(:))/max(abs(eta_unit(:)))),1e-15);
p_error = max(max(abs(wvt.p(:)-p_unit(:))/max(abs(p_unit(:)))),1e-15);
end
