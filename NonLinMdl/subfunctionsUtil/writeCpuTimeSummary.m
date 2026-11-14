function [S, fileTxt] = writeCpuTimeSummary(tictoc_LPVMPC, strWindType, DT, strSuffix, outDir)
%WRITECPUTIMESUMMARY CPU time statistics of the qLPV-MPC and thesis text
%   [S, fileTxt] = writeCpuTimeSummary(tictoc_LPVMPC, strWindType, DT)
%   evaluates the CPU time per controller cycle logged by the qLPV-MPC
%   (timeseries tictoc_LPVMPC of runCompareCtrl / getSimulationOutputTable):
%   number of cycles, mean and maximum CPU time, and the number of cycles
%   that exceed 3 ms, 10 ms and the sampling time DT (default 0.008 s).
%   The statistics and a replacement of the real-time paragraph of the
%   thesis (Section qLMPC, Figure fig:CtrlCmpMdlSweep) are written to the
%   text file cputime_MPC<wind case><strSuffix>.txt in outDir (default:
%   current folder, as the results_table*.tex of runCompareCtrl) and to the
%   command window.
%   strWindType: 'Sweep', 'Step', 'EOG...' or the mean wind speed of the
%   turbulent wind (e.g. 18).

if nargin < 3 || isempty(DT), DT = 0.008; end
if nargin < 4, strSuffix = ''; end
if nargin < 5 || isempty(outDir), outDir = pwd; end

x = double(tictoc_LPVMPC.Data(:));
S.nCycles = numel(x);
S.duration = tictoc_LPVMPC.Time(end) - tictoc_LPVMPC.Time(1);
S.mean = mean(x);
S.max = max(x);
S.n3ms = sum(x > 3e-3);
S.n10ms = sum(x > 10e-3);
S.nDT = sum(x > DT);

if isa(strWindType, 'double')
    strCase = sprintf('NTW%02d', strWindType);
    strWind = sprintf('the turbulent wind field with a mean wind speed of %g\\,m/s', strWindType);
elseif strcmp(strWindType, 'Sweep')
    strCase = strWindType;
    strWind = 'the stepped wind sweep';
else
    strCase = strWindType;
    strWind = sprintf('the wind case %s', strWindType);
end

txt = sprintf(['The QP is solved using a linear complementarity problem formulation with Lemke''s algorithm ' ...
    '\\cite{cisneros2020velocity}, warm-started with the solution of the previous sample, which enables fast ' ...
    'computation, as real-time capability is a critical requirement for online optimization-based control. ' ...
    'Of the %d controller cycles at a sampling time of %g\\,ms needed for the %g\\,s simulation time of %s, ' ...
    '%s and %s. The mean CPU time was approximately %s\\,ms, as shown in the subplot captions of ' ...
    'Figure~\\ref{fig:CtrlCmpMdlSweep}. These measurements were obtained in MATLAB R2023a on an Intel Core ' ...
    'i7-11850H (2.5\\,GHz, 8 cores) running Windows 10. Note that for deployment on a physical turbine, a ' ...
    'fallback strategy would be required to handle the rare cases where the solver exceeds the sampling time. ' ...
    'This was not implemented in this proof-of-concept, but several options exist, such as blending the last ' ...
    'available qLMPC solution with the output of the PI feedback controller.'], ...
    S.nCycles, 1000*DT, round(S.duration), strWind, exceedText(S.n3ms, 3, true), exceedText(S.n10ms, 10, false), ...
    num2str(round(1000*S.mean, 2)));

fileTxt = fullfile(outDir, ['cputime_MPC', strCase, strSuffix, '.txt']);
fid = fopen(fileTxt, 'w');
closeF = onCleanup(@() fclose(fid));
for f = [1, fid]
    fprintf(f, 'CPU time qLPV-MPC, wind case %s (%s)\n', strCase, char(datetime('now', 'Format', 'yyyy-MM-dd')));
    fprintf(f, 'Controller cycles:        %d (%g s, sampling time %g ms)\n', S.nCycles, S.duration, 1000*DT);
    fprintf(f, 'Mean CPU time:            %.3f ms\n', 1000*S.mean);
    fprintf(f, 'Max. CPU time:            %.3f ms\n', 1000*S.max);
    fprintf(f, 'Cycles > 3 ms:            %d\n', S.n3ms);
    fprintf(f, 'Cycles > 10 ms:           %d\n', S.n10ms);
    fprintf(f, 'Cycles > sampling time:   %d\n', S.nDT);
    fprintf(f, '\nReplacement of the real-time paragraph (LaTeX):\n\n%s\n', txt);
end
fprintf('CPU time summary written to %s\n', fileTxt);
end

function s = exceedText(n, lim, first)
% "only eight cycles exceeded 3\,ms", "only one exceeded 10\,ms", "none exceeded ..."
words = {'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve'};
if n == 0
    s = sprintf('none exceeded %d\\,ms', lim);
    if first, s = sprintf('none of them exceeded %d\\,ms', lim); end
    return
end
if n <= numel(words), nStr = words{n}; else, nStr = sprintf('%d', n); end
if first
    if n == 1, unit = 'cycle'; else, unit = 'cycles'; end
    s = sprintf('only %s %s exceeded %d\\,ms', nStr, unit, lim);
else
    s = sprintf('only %s exceeded %d\\,ms', nStr, lim);
end
end
