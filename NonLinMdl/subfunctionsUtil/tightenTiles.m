function tightenTiles(ax, margin)
%TIGHTENTILES Same compact spacing for stacked time series tiles
%   tightenTiles(ax) prepares the axes ax of a tiled layout with
%   'TileSpacing' 'compact' (top to bottom, shared time axis):
%   - x tick labels only on the bottom axes
%   - no y ticks at or within margin*(y range) of the y limits (default
%     margin 0.07, about half the height of a tick label in a tile):
%     their labels extend beyond the axes box, and the tiled layout then
%     inserts a gap between the tiles, so that figures with
%     y limits on tick values (e.g. [43 44]) had a wider spacing than others;
%     one tile with labels at its limits widens the gaps of all tiles. If
%     fewer than two ticks remain, half the tick spacing is used.
%   Call it after the final y limits are set; it can be called again after
%   the y limits change (the y ticks are recomputed).

if nargin < 2, margin = 0.07; end
ax = ax(isgraphics(ax, 'axes'));
set(ax(1:end-1), 'XTickLabel', []);
for k = 1:numel(ax)
    ax(k).YTickMode = 'auto';
    yl = ax(k).YLim;
    t = ax(k).YTick;
    d = margin*diff(yl);
    tIn = t(t > yl(1) + d & t < yl(2) - d);
    if numel(tIn) < 2 && numel(t) >= 2
        % e.g. limits [43 44] with the ticks 43, 43.5, 44: half the tick
        % spacing (43.25, 43.5, 43.75)
        step = (t(2) - t(1))/2;
        tFine = round((ceil(yl(1)/step):floor(yl(2)/step))*step, 10);
        tIn = tFine(tFine > yl(1) + d & tFine < yl(2) - d);
    end
    if numel(tIn) >= 2
        ax(k).YTick = tIn;
    end
end
end
