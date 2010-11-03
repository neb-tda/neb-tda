function PlotZeroCells(dataSet, boolSaveFigure, basis1, basis2)
% Plots the data and its 0-cells. If basis1 and basis2 are not both
% specified, then a planar projection is computed using principal component
% analysis.
%   
% Input:
%   dataSet - name of the data set
%   boolSaveFigure - if 1, then the figure is saved
%   basis1 - first basis vector
%   basis2 - second basis vector
%	
% Usage:
%   PlotZeroCells('GeneExpressions', 1)
%   PlotZeroCells('RangeImagePatches', 1, 1, 5)
%   PlotZeroCells('OpticalFlowPatches', 1, 1, 2)
%   PlotZeroCells('OpticalImagePatches', 1, 1, 2)
%   PlotZeroCells('SocialNetwork', 1)

data = importdata([dataSet, '.txt']);
zeroCells = importdata([dataSet, '_0cells_ordered.txt']);

if nargin < 4
    [coeff, score] = princomp(data);
    m = mean(data);
    data = score;
    zeroCells = (zeroCells - ones(size(zeroCells, 1), 1) * m) * coeff;
    basis1 = 1;
    basis2 = 2;
end

h = figure;
hold on;

markerSize = 5;
if strcmp(dataSet, 'GeneExpressions')
    markerSize = 15;
end
plot(data(:, basis1), data(:, basis2), '.', 'MarkerSize', markerSize, 'MarkerEdgeColor', [.5, .5, .5]);

plot(zeroCells(:, basis1), zeroCells(:, basis2), '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'k');

axis equal;
axis off;
if boolSaveFigure == 1
    saveas(h, [dataSet, '_0cells.fig'])
end