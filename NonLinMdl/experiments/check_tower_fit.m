function R = check_tower_fit(Vlist)
% Frequency response (magnitude and phase) at the tower mode: FAST
% linearization vs. Mdl1, Mdl2 (compareLinearModels) and the reduced tower
% model of the qLPV-MPC (ABqL_WT_tower, without actuators)
if nargin < 1, Vlist = [12 15 18 22]; end
repo = fileparts(fileparts(fileparts(mfilename('fullpath')))); % repository root
addpath(fullfile(repo, 'LinMdl'), fullfile(repo, 'dataIn'), fullfile(repo, 'NonLinMdl', 'subfunctionsMPC'));
idxV = Vlist - 3;
[~, ~, sysCmp] = compareLinearModels(idxV, fullfile(tempdir, 'bodeTmp'), 0, 0, 0, 'off', 2);

% Tables and parameters of the controller (DD_Mdl1)
dd = Simulink.data.dictionary.open(fullfile(repo, 'NonLinMdl', 'WECSMdl', 'DD_Mdl1.sldd'));
s = getSection(dd, 'Design Data');
Cq = getValue(getEntry(s, 'Cq')); Ct = getValue(getEntry(s, 'Ct'));
lam = getValue(getEntry(s, 'lambda')); bet = getValue(getEntry(s, 'beta'));
wecs = getValue(getEntry(s, 'wecs'));
load(fullfile(repo, 'dataIn', 'NREL5MW_linearised_4to25.mat'), 'Lin');

f = [0.1 0.2 0.25 0.3 0.32 0.35 0.4 0.5 0.7];
w = 2*pi*f;
R = struct([]);
for k = 1:numel(idxV)
    S = sysCmp{idxV(k)};
    clear ABqL_WT_tower % persistent tables
    [A, B] = ABqL_WT_tower([Lin.RSpeed(idxV(k)), Lin.Pitch(idxV(k)), Lin.V(idxV(k))], lam(:)', bet(:)', Cq, Ct, wecs);
    % without actuators: inputs Tg (kNm), beta (rad), V; states 1-3, 6-7
    xs = [1 2 3 6 7];
    Bm = [A(xs, 4)*1000, A(xs, 5), B(xs, 1)];
    Am = A(xs, xs);
    Cm = [0 1 0 0 0; Am(5, :)]; Dm = [0 0 0; Bm(5, :)];
    S.MPC = ss(Am, Bm, Cm, Dm);
    clear ABqL_WT_flap
    [A, B] = ABqL_WT_flap([Lin.RSpeed(idxV(k)), Lin.Pitch(idxV(k)), Lin.V(idxV(k))], lam(:)', bet(:)', Cq, Ct, wecs);
    xs = [1 2 3 6 7 8 9];
    Am = A(xs, xs); Bm = [A(xs, 4)*1000, A(xs, 5), B(xs, 1)];
    S.FLAP = ss(Am, Bm, [0 1 0 0 0 0 0; Am(6, :)], [0 0 0; Bm(6, :)]);
    names = {'FAST', 'Mdl1', 'Mdl2', 'MPC', 'FLAP'};
    fprintf('\nV = %g m/s: pitch -> tower fore-aft acc. / rotor speed\n', S.V);
    fprintf('%-6s', 'f[Hz]'); fprintf('%8.2f', f); fprintf('\n');
    for io = [2 1]
        for n = 1:numel(names)
            H = squeeze(freqresp(S.(names{n})(io, 2), w));
            if n == 1, H0 = H; end
            outLbl = {'w ', 'a '}; fprintf('%-4s %s|H|  ', names{n}, outLbl{io}); fprintf('%8.3g', abs(H)); fprintf('\n');
            fprintf('%-4s    phase', names{n}); fprintf('%8.0f', rad2deg(angle(H ./ H0))); fprintf('   (rel. to FAST)\n');
        end
    end
    % tower mode: frequency and damping
    for n = 1:numel(names)
        [wn, z] = damp(S.(names{n}));
        sel = wn/(2*pi) > 0.2 & wn/(2*pi) < 0.45;
        fprintf('%-4s tower-range modes: %s Hz, damping %s %%\n', names{n}, mat2str(unique(round(wn(sel)'/(2*pi), 3))), mat2str(unique(round(100*z(sel)', 1))));
    end
    R(k).S = S; %#ok<AGROW>
end
close(dd);
end
