function PlotOneCells(dataSet, numBandPoints, boolSaveFigure, basis1, basis2)
% Plots the data, its 0-cells, and its 1-cells. If basis1 and basis2 are 
% not both specified, then a planar projection is computed using principal 
% component analysis.
%   
% Input:
%   dataSet - name of the data set
%   numBandPoints - number of points in each band
%   boolSaveFigure - if 1, then the figure is saved
%   basis1 - first basis vector
%   basis2 - second basis vector
%	
% Usage:
%   PlotOneCells('GeneExpressions', 11, 1)
%   PlotOneCells('RangeImagePatches', 11, 1, 1, 5)
%   PlotOneCells('OpticalFlowPatches', 11, 1, 1, 2)
%   PlotOneCells('OpticalImagePatches', 11, 1, 1, 2)
%   PlotOneCells('SocialNetwork', 11, 1)

data = importdata([dataSet, '.txt']);
zeroCells = importdata([dataSet, '_0cells_ordered.txt.']);
oneCells = importdata([dataSet, '_1cells_ordered.txt']);

numLines = size(oneCells, 1);
if mod(numLines, numBandPoints) ~= 0
    error('The number of points in the file must be a multiple of the number of points in each band');
end
numBands = numLines / numBandPoints;

if nargin < 4
    [coeff, score] = princomp(data);
    m = mean(data);
    data = score;
    zeroCells = (zeroCells - ones(size(zeroCells, 1), 1) * m) * coeff;
    oneCells = (oneCells - ones(numLines, 1) * m) * coeff;
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

for i = 1 : numBands
    oneCell = oneCells((i - 1) * numBandPoints + 1 : i * numBandPoints, :);
    plot(oneCell(:, basis1), oneCell(:, basis2), 'Color', [0, 0, 1], 'LineStyle', '-', 'LineWidth', 3);
end

plot(zeroCells(:, basis1), zeroCells(:, basis2), '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'k');

axis equal;
axis off;
if boolSaveFigure == 1
    saveas(h, [dataSet, '_1cells.fig'])
end