function plotParetoVerif2(rateBar)
%PLOTPARETOVERIF2 Pareto and bar plot of the evaluation criteria (18 m/s, Verif. 2)
%   plotParetoVerif2 plots the table figDirPareto/ParetoVerif2.mat of
%   runParetoVerif2 (no simulations):
%   - Pareto plot in three panels, (a) qLMPC on Mdl2 and FAST with the IPC
%     and AALC pale, (b) IPC and AALC with the qLMPC of (a) pale, (c) qLMPC
%     + IPC on FAST with the qLMPC on FAST and the IPC pale: blade DEL (blue), tower DEL (green) and power std
%     (purple) over the actuator power C_Act for the qLMPC on Mdl2 (dotted,
%     hexagrams), the qLMPC designs flapNotch (solid, diamonds),
%     flapNotchIPC (solid, triangles), flapNotchSmooth12 and
%     flapNotchSmooth12IPC (dash-dotted, down/right triangles) on FAST, the
%     IPC (dashed, circles) at 4, 8 and 13 deg/s and the AALC (squares,
%     8 deg/s). Marker fill: 4 deg/s white, 8 deg/s grey, 13 deg/s black.
%     The qLMPC on Mdl2 has no loads and no pitch moment: tower fore-aft
%     acceleration as proxy for the tower DEL, C_Act with the estimated
%     pitch moment (C_ActEst of runParetoVerif2).
%   - bar plot of the criteria blade DEL (blue), tower DEL (green), C_P
%     (grey), C_dP (purple) and C_Act (red) at the rate limit rateBar (default
%     8 deg/s), grouped by controller (as fcnPlotBar), one panel per Pareto
%     panel with its bright controllers.
%   Figures: figDirPareto/afig01_IPCMPCpareto and afig01_IPCMPCbar (.eps
%   and .png, names of the thesis figures in figures/figures_WESC25).

if nargin < 1 || isempty(rateBar)
    rateBar = 8;
end
mainDir = fileparts(mfilename('fullpath'));
outDir = fullfile(mainDir, 'figDirPareto');
S = load(fullfile(outDir, 'ParetoVerif2.mat'), 'ParetoTable');
T = S.ParetoTable;

colB = [0, 0.4470, 0.7410];      % blade DEL
colT = [0.4660, 0.6740, 0.1880]; % tower DEL
colP = [0.4940, 0.1840, 0.5560]; % power std
lw = 1.5;
ms = 7;

% Controllers (qLMPC designs, IPC, AALC): table name, legend label, line
% style, marker
ctrl = {"qLMPC", 'qLMPC (Mdl2, proxies)', ':', 'h';
    "qLMPC flapNotch", 'qLMPC (FAST)', '-', 'd';
    "qLMPC flapNotchIPC", 'qLMPC + IPC (FAST)', '-', '^';
    "qLMPC flapNotchSmooth12", 'qLMPC smooth (FAST)', '-.', 'v';
    "qLMPC flapNotchSmooth12IPC", 'qLMPC smooth + IPC (FAST)', '-.', '>';
    "IPC", 'IPC', '--', 'o';
    "AALC", 'AALC', 'none', 's'};
% Marker fill by pitch rate limit: 4 deg/s white, 8 deg/s grey, 13 deg/s black
rateFill = [4, 8, 13];
fillCol = {'w', [0.6 0.6 0.6], 'k'};

%% Pareto plot: (a) qLMPC with IPC and AALC pale, (b) IPC and AALC with
% the qLMPC pale, (c) qLMPC + IPC with the qLMPC on FAST and the IPC pale
qLMPC = ["qLMPC", "qLMPC flapNotch", "qLMPC flapNotchSmooth12"];
panels = {qLMPC, ["IPC", "AALC"], ["qLMPC flapNotchIPC", "qLMPC flapNotchSmooth12IPC"]};
panelsPale = {["IPC", "AALC"], qLMPC, ["qLMPC flapNotch", "qLMPC flapNotchSmooth12", "IPC"]};
titles = {'(a) qLMPC (pale: IPC and AALC)', '(b) IPC and AALC (pale: qLMPC)', ...
    '(c) qLMPC + IPC (pale: qLMPC and IPC alone)'};
nP = numel(panels);
yMax = max(1.1, ceil(10*max([T.C_DEL_b; T.C_DEL_t; T.C_dP; T.C_TwrAcc], [], 'omitnan'))/10);
xMax = ceil(max([T.C_Act; T.C_ActEst], [], 'omitnan'));
figure(1); clf;
set(gcf, 'Name', 'Pareto Verif. 2, 18 m/s', 'NumberTitle', 'off', 'Position', [50, 20, 560, 966]); % 2.3 times the default figure height (420)
tl = tiledlayout(nP, 1, 'TileSpacing', 'compact');
for idxP = 1:nP
    ax = nexttile(tl);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    for k = 1:size(ctrl, 1)
        if any(ctrl{k, 1} == panelsPale{idxP})
            plotController(ax, T, ctrl(k, :), {colB, colT, colP}, rateFill, fillCol, lw, ms, 0.7);
        end
    end
    for k = 1:size(ctrl, 1)
        if any(ctrl{k, 1} == panels{idxP})
            plotController(ax, T, ctrl(k, :), {colB, colT, colP}, rateFill, fillCol, lw, ms, 0);
        end
    end
    title(ax, titles{idxP}, 'FontWeight', 'normal');
    ylabel(ax, 'Evaluation criteria (-)');
    xlim(ax, [1, xMax]);
    ylim(ax, [0.3, yMax]);
    set(ax, 'FontSize', 11);
    % Arrow to the lower left corner: lower criteria and C_Act are better
    quiver(ax, 1.9, 0.62, -0.65, -0.24, 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
        'MaxHeadSize', 0.6, 'HandleVisibility', 'off');
    text(ax, 1.62, 0.62, 'better', 'Color', [0.3 0.3 0.3], 'FontSize', 10, ...
        'FontAngle', 'italic', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
xlabel(ax, 'C_{Act} (-) (Mdl2: estimated pitch moment)');
% Legend: controllers in black, criteria in colour, rate limits as marker fill
hL = gobjects(0);
for k = 1:size(ctrl, 1)
    hL(end+1) = plot(ax, NaN, NaN, 'LineStyle', ctrl{k, 3}, 'Marker', ctrl{k, 4}, 'Color', 'k', ...
        'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', ctrl{k, 2}); %#ok<AGROW>
end
hL(end+1) = plot(ax, NaN, NaN, '-', 'Color', colB, 'LineWidth', 3, 'DisplayName', 'C_{DEL,b}');
hL(end+1) = plot(ax, NaN, NaN, '-', 'Color', colT, 'LineWidth', 3, 'DisplayName', 'C_{DEL,t} (Mdl2: tower acc.)');
hL(end+1) = plot(ax, NaN, NaN, '-', 'Color', colP, 'LineWidth', 3, 'DisplayName', 'C_{\delta P}');
for i = 1:numel(rateFill)
    hL(end+1) = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', fillCol{i}, ...
        'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', sprintf('%g\\circ/s', rateFill(i))); %#ok<AGROW>
end
lgd = legend(hL, 'NumColumns', 3, 'FontSize', 10, 'Box', 'off');
lgd.Layout.Tile = 'north';
saveFig(gcf, fullfile(outDir, 'afig01_IPCMPCpareto'));

%% Bar plot at rateBar, grouped by controller (colours as fcnPlotBar),
% one panel per Pareto panel with the controllers shown bright there
% Order of the thesis tables: blade DEL, tower DEL, mean power, power std,
% actuator power
crit = {'C_DEL_b', 'C_DEL_t', 'C_P', 'C_dP', 'C_Act'};
critLabel = {'C_{DEL,b}', 'C_{DEL,t}', 'C_{P}', 'C_{\delta P}', 'C_{Act}'};
cl1 = lines(7);
critCol = [cl1(1, :); cl1(5, :); 0.5, 0.5, 0.5; cl1(4, :); 1, 0, 0];
B = T(T.RateLimit_degps == rateBar, :);
% Mdl2: C_Act with the estimated pitch moment, tower acceleration as proxy
% for the tower DEL
isMdl2 = B.Controller == "qLMPC";
B.C_Act(isMdl2) = B.C_ActEst(isMdl2);
B.C_DEL_t(isMdl2) = B.C_TwrAcc(isMdl2);
yMaxBar = 1.15*max(B{:, crit}, [], 'all', 'omitnan');
titlesBar = {'(a) qLMPC', '(b) IPC and AALC', '(c) qLMPC + IPC'};

figure(2); clf;
set(gcf, 'Name', sprintf('Criteria %g deg/s, Verif. 2, 18 m/s', rateBar), 'NumberTitle', 'off', ...
    'Position', [650, 20, 560, 966]);
tl = tiledlayout(nP, 1, 'TileSpacing', 'compact');
for idxP = 1:nP
    ax = nexttile(tl);
    names = panels{idxP};
    [~, pos] = ismember(names, B.Controller);
    names = names(pos > 0);
    K = B{pos(pos > 0), crit};
    [~, idxC] = ismember(names, [ctrl{:, 1}]);
    tickLabel = ctrl(idxC, 2);
    tickLabel(names == "qLMPC") = {'qLMPC (Mdl2)^*'};
    b = bar(ax, K);
    for k = 1:numel(crit)
        b(k).FaceColor = critCol(k, :);
        ok = ~isnan(b(k).YData);
        text(ax, b(k).XEndPoints(ok), b(k).YEndPoints(ok), compose('%.2f', b(k).YData(ok)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
    yline(ax, 1, 'k--', 'HandleVisibility', 'off');
    grid(ax, 'on');
    ylim(ax, [0, yMaxBar]);
    ylabel(ax, 'Evaluation criteria (-)');
    title(ax, titlesBar{idxP}, 'FontWeight', 'normal');
    set(ax, 'FontSize', 11);
    ax.XTickLabel = tickLabel;
    ax.XAxis.FontSize = 10;
    ax.XTickLabelRotation = 0;
end
lgd = legend(b, critLabel, 'NumColumns', numel(crit), 'FontSize', 10, 'Box', 'off');
lgd.Layout.Tile = 'north';
saveFig(gcf, fullfile(outDir, 'afig01_IPCMPCbar'));
end

function plotController(ax, T, c, colCrit, rateFill, fillCol, lw, ms, pale)
% plotController criteria of the controller c = {table name, label, line
% style, marker} over C_Act (Mdl2: proxies over the pitch rate ratio), one
% marker per rate limit, colours mixed with white by the factor pale
R = sortrows(T(T.Controller == c{1}, :), 'RateLimit_degps');
if isempty(R), return, end
if c{1} == "qLMPC" % Mdl2: proxies
    x = R.C_ActEst;
    ys = {R.C_TwrAcc, R.C_dP};
    cols = colCrit(2:3);
else
    x = R.C_Act;
    ys = {R.C_DEL_b, R.C_DEL_t, R.C_dP};
    cols = colCrit;
end
mix = @(col) (1 - pale)*validatecolor(col) + pale;
for j = 1:numel(ys)
    plot(ax, x, ys{j}, 'LineStyle', c{3}, 'Color', mix(cols{j}), 'LineWidth', lw, ...
        'HandleVisibility', 'off');
    for i = 1:numel(x)
        plot(ax, x(i), ys{j}(i), 'LineStyle', 'none', 'Marker', c{4}, 'Color', mix(cols{j}), ...
            'MarkerFaceColor', mix(fillCol{rateFill == R.RateLimit_degps(i)}), 'LineWidth', lw, ...
            'MarkerSize', ms, 'HandleVisibility', 'off');
    end
end
end

function saveFig(h, fName)
print(h, fName, '-depsc');
print(h, fName, '-dpng', '-r150');
fprintf('Saved %s.eps/.png\n', fName);
end
