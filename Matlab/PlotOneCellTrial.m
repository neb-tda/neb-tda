function PlotOneCellTrial(dataSet, numBandPoints, boolSaveFigure, basis1, basis2)
% Plots the force norms and path of the 1-cell for a single 1-cell trial.
% The initial 1-cell is in red and the final 1-cell is in blue. The lines 
% which fade from yellow to green trace the paths of the intermediate
% nodes. If basis1 and basis2 are not both specified, then a planar 
% projection is computed using principal component analysis.
%   
% Input:
%   dataSet - name of the data set
%   numBandPoints - number of points in each band
%   boolSaveFigure - if 1, then the figure is saved
%   basis1 - first basis vector
%   basis2 - second basis vector
%	
% Usage:
%   PlotOneCellTrial('GeneExpressions', 11, 1)
%   PlotOneCellTrial('RangeImagePatches', 11, 1, 1, 5)
%   PlotOneCellTrial('OpticalFlowPatches', 11, 1, 1, 2)
%   PlotOneCellTrial('OpticalImagePatches', 11, 1, 1, 2)
%   PlotOneCellTrial('SocialNetwork', 11, 1)

forces = importdata([dataSet, '_1forces.txt']);
h1 = figure;
plot(forces)
h2 = legend('gradient force norm', 'spring force norm', 'angle force norm', 'total force norm', 4, 'Location', 'NorthEast');
set(h2, 'Interpreter', 'none')

data = importdata([dataSet, '.txt']);
zeroCells = importdata([dataSet, '_0cells_ordered.txt.']);
points = importdata([dataSet, '_1points.txt.']);

numLines = size(points, 1);
if mod(numLines, numBandPoints) ~= 0
    error('The number of points in the file must be a multiple of the number of points in each band');
end
numSteps = numLines / numBandPoints;

if nargin < 4
    [coeff, score] = princomp(data);
    m = mean(data);
    data = score;
    zeroCells = (zeroCells - ones(size(zeroCells, 1), 1) * m) * coeff;
    points = (points - ones(numLines, 1) * m) * coeff;
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

initialBandColor = [1,0,0];
finalBandColor = [0,0,1];
initialNodesColor = [1,1,0];
finalNodesColor = [0,1,0];
for i = 2 : numSteps - 1
  band = points((i - 1) * numBandPoints + 1 : i * numBandPoints, :);
  plot(band(:, basis1), band(:, basis2), 'Color', ((numSteps - 1 - i) * initialNodesColor + (i - 2) * finalNodesColor) / (numSteps - 3), 'LineStyle', 'none', 'Marker', '.');
end
band = points(1 : numBandPoints, :);
plot(band(:, basis1), band(:, basis2), 'Color', initialBandColor, 'LineStyle', '-', 'LineWidth', 3);
band = points((numSteps - 1) * numBandPoints + 1 : numSteps * numBandPoints, :);
plot(band(:, basis1), band(:, basis2), 'Color', finalBandColor, 'LineStyle', '-', 'LineWidth', 3); 

plot(zeroCells(:, basis1), zeroCells(:, basis2), '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'k');

axis equal;
axis off;
if boolSaveFigure == 1
    saveas(h1, [dataSet, '_1forces.fig'])
    saveas(h, [dataSet, '_1points.fig'])
end