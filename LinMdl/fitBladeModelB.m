function P = fitBladeModelB(wecs, lam, pit, Ct, rho, fastMat)
% fitBladeModelB: quasi-steady fit of the variant B blade model to FAST.
%
% Variant B: rigid blade, rotation theta about the root, modal mass mb at
% the fixed reference radius r_ref, thrust applied at the load point
% r_b(lambda, beta), stiffness Kb(lambda, beta). Tip OoP deflection R*theta.
% Static: theta = F r_b / (N Kb r_ref^2)  ->  tip = R r_b F/(N Kb r_ref^2).
%
% Steady states of the FAST wind steps (4 ... 25 m/s):
% - aerodynamic thrust F_aero = RotThrust - m_r g sin(tilt) (RotThrust is the
%   shaft force and contains the rotor weight component, 94 kN)
% - r_b: root moment arm sum(RootMOoP)/F_aero above rated; below rated
%   constant at the rated value (the low-wind values contain the
%   non-aerodynamic moment of precone, centrifugal and gravity loads)
% - Kb such that the static tip deflection with the thrust of the model
%   (Ct table at the FAST operating point) equals the FAST tip deflection
% - r_ref such that Kb at rated (11 m/s) equals the current Kb (0.699 Hz)
%
% Tables: below rated over lambda (beta = fine pitch), above rated over
% beta; f(lambda, beta) = fLam(lambda) + fBeta(beta) - fRated.
if nargin < 6, fastMat = 'OutDataStepFASTTool.mat'; end
% cache: the fit depends on the FAST data, the thrust table and the
% blade/wind parameters only
persistent cacheKey cacheP
key = {fastMat, wecs.Kb, wecs.mb, wecs.kWindEff, rho, size(Ct), sum(Ct(:)), sum(lam(:)), sum(pit(:))};
if isequal(key, cacheKey)
    P = cacheP;
    return
end
S = load(fastMat, 'OutTable'); T = S.OutTable;
T = T(T.Time >= 100, :); T.Time = T.Time - T.Time(1);
R = wecs.Rr; N = wecs.N;
Fg = wecs.mr*9.81*sind(5);
n = 22; d = zeros(n, 6);
for k = 1:n
    i = T.Time <= k*100-10 & T.Time >= k*100-30;
    V = mean(T.Wind1VelX(i)); w = mean(T.RotSpeed(i))*2*pi/60; b = mean(T.BlPitch1(i))*pi/180;
    Fa = mean(T.RotThrust(i))*1e3 - Fg;
    Mr = mean(T.RootMOoP1(i) + T.RootMOoP2(i) + T.RootMOoP3(i))*1e3;
    tip = mean((T.OoPDefl1(i) + T.OoPDefl2(i) + T.OoPDefl3(i))/3);
    l = w*R/(wecs.kWindEff*V);
    Fm = 0.5*rho*pi*R^2*interp2(lam, pit, Ct, l, b)*(wecs.kWindEff*V)^2;
    d(k, :) = [V, l, b, Mr/Fa, tip, Fm];
end
P.data = array2table(d, 'VariableNames', {'V', 'lambda', 'beta', 'rbFAST', 'tipFAST', 'Fmdl'});
iRated = find(d(:, 1) > 10.5 & d(:, 1) < 11.5, 1); % 11 m/s
rb = d(:, 4);
rb(1:iRated) = rb(iRated);
G = d(:, 5)./d(:, 6); % tip per model thrust (m/N)
P.rRef = sqrt(R*rb(iRated)/(N*wecs.Kb*G(iRated)));
Kb = R*rb./(N*P.rRef^2*G);
P.data.rb = rb; P.data.Kb = Kb; P.data.fb = sqrt(Kb/wecs.mb)/2/pi;

% Tables (sorted, unique breakpoints)
[P.lamBp, iL] = unique(d(1:iRated, 2));
P.rbLam = rb(iL); P.KbLam = Kb(iL);
[P.betaBp, iB] = unique(d(iRated:end, 3));
iB = iB + iRated - 1;
P.rbBeta = rb(iB); P.KbBeta = Kb(iB);
P.rbRated = rb(iRated); P.KbRated = Kb(iRated);
cacheKey = key; cacheP = P;
end
