function printFigure(h, fName, formats)
%PRINTFIGURE Export a figure as vector .eps and 300 dpi .png with a fixed size
%   printFigure(h, fName) saves figure h (handle or number) as fName.eps
%   and fName.png with the size set by setFigureSize (default: the
%   standard figure, 14.82 x 11.11 cm). The output is laid out for this
%   size, not for the window, which Windows shrinks to the screen. The eps
%   is a vector graphic (-painters; otherwise MATLAB embeds a bitmap for
%   plots with many points) with the full figure as bounding box (-loose),
%   so that figures of the same size have the same size in LaTeX.
%   printFigure(h, fName, formats) saves only the formats in the cell
%   array formats ('eps', 'png').

if nargin < 3 || isempty(formats)
    formats = {'eps', 'png'};
end
if isnumeric(h)
    h = figure(h);
end
if ~strcmp(get(h, 'PaperPositionMode'), 'manual') || ~strcmp(get(h, 'PaperUnits'), 'centimeters')
    setFigureSize(h, [1 1]);
end
if any(strcmp(formats, 'eps'))
    print(h, fName, '-depsc', '-painters', '-loose');
end
if any(strcmp(formats, 'png'))
    print(h, fName, '-dpng', '-r300');
end
end
