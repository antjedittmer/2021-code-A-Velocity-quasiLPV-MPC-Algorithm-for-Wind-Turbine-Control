function [rr, rb, Kb, mb, Bb, gSched] = bladeParamsMdl2(wecs, lam0, beta, Ft)
% bladeParamsMdl2 returns the blade parameters of Mdl2 at the operating
% point: Simulink models Mdl_BianchiOL_BladeB and
% Mdl_BianchiOL_BladeB_WindFiltered (MATLAB Function block BladeSchedule of
% the WECS subsystem) and the linear model of compareLinearModels.
%
% Original blade model (thesis, wecs without the fields blade*): modal
% mass, stiffness and thrust load point at rr = rb = wecs.rb.
% Variant B (fields blade* of wecs, initModel5MWNREL with bladeModelVersion
% B1 or B2, getBladeModelVersion; fit: fitBladeModelB): modal mass at the
% reference radius rr = wecs.bladeRRef, load point rb(lambda0, beta),
% stiffness Kb(lambda0, beta), B2: mb = Kb/wnb^2. lambda0 = Rr*OmegaR/(Ve -
% ydot_fa) is the tip speed ratio without the blade velocity.
%
% gSched: derivative of the generalized blade force rb*Ft - N*Kb*rr^2*delta
% at the equilibrium delta = Ft*rb/(N*Kb*rr^2) with respect to
% [lambda0, beta] (zero for the original model); Ft: thrust (N).
if ~isfield(wecs, 'bladeRRef')
    rr = wecs.rb; rb = wecs.rb; Kb = wecs.Kb; mb = wecs.mb; Bb = wecs.Bb;
    gSched = [0, 0];
    return
end
rr = wecs.bladeRRef;
[rb, Kb] = schedB(wecs, lam0, beta);
if wecs.bladeMassScaled > 0
    mb = Kb/wecs.wnb^2;
else
    mb = wecs.mb;
end
Bb = 2*wecs.zetab*sqrt(Kb*mb);
delta = Ft*rb/(wecs.N*Kb*rr^2);
dl = 1e-3; db = 1e-5;
[rbp, Kbp] = schedB(wecs, lam0 + dl, beta); [rbm, Kbm] = schedB(wecs, lam0 - dl, beta);
gL = ((rbp - rbm)*Ft - wecs.N*(Kbp - Kbm)*rr^2*delta)/(2*dl);
[rbp, Kbp] = schedB(wecs, lam0, beta + db); [rbm, Kbm] = schedB(wecs, lam0, beta - db);
gB = ((rbp - rbm)*Ft - wecs.N*(Kbp - Kbm)*rr^2*delta)/(2*db);
gSched = [gL, gB];
end

function [rb, Kb] = schedB(wecs, lam0, beta)
% tables f(lambda, beta) = fLam(lambda) + fBeta(beta) - fRated (bladeScheduleB)
l = min(max(lam0, wecs.bladeLamBp(1)), wecs.bladeLamBp(end));
b = min(max(beta, wecs.bladeBetaBp(1)), wecs.bladeBetaBp(end));
rb = interp1(wecs.bladeLamBp, wecs.bladeRbLam, l) + interp1(wecs.bladeBetaBp, wecs.bladeRbBeta, b) - wecs.bladeRbRated;
Kb = interp1(wecs.bladeLamBp, wecs.bladeKbLam, l) + interp1(wecs.bladeBetaBp, wecs.bladeKbBeta, b) - wecs.bladeKbRated;
end
