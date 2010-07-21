function PlotZeroCellTrial(dataSet, boolSaveFigure, basis1, basis2)
% Plots the force norms and path of the 0-cell for a single 0-cell trial.
% The color of the 0-cell fades from red to blue as it moves from its 
% initial to final location. If basis1 and basis2 are not both specified, 
% then a planar projection is computed using principal component analysis.
%   
% Input:
%   dataSet - name of the data set
%   boolSaveFigure - if 1, then the figure is saved
%   basis1 - first basis vector
%   basis2 - second basis vector
%	
% Usage:
%   PlotZeroCellTrial('GeneExpressions', 1)
%   PlotZeroCellTrial('RangeImagePatches', 1, 1, 5)
%   PlotZeroCellTrial('OpticalFlowPatches', 1, 1, 2)
%   PlotZeroCellTrial('OpticalImagePatches', 1, 1, 2)
%   PlotZeroCellTrial('SocialNetwork', 1)

forces = importdata([dataSet, '_0forces.txt']);
h1 = figure;
plot(forces);
h2 = legend('gradient force norm', 1, 'Location', 'NorthEast');
set(h2, 'Interpreter', 'none');

data = importdata([dataSet, '.txt']);
points = importdata([dataSet, '_0points.txt.']);

if nargin < 4
    [coeff, score] = princomp(data);
    m = mean(data);
    data = score;
    points = (points - ones(size(points, 1), 1) * m) * coeff;
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

numSteps = size(points, 1);
startColor = [1,0,0];
endColor = [0,0,1];
for i = 1 : numSteps
    plot(points(i, basis1), points(i, basis2), '.', 'MarkerSize', 5, 'MarkerEdgeColor', ((numSteps - i) * startColor + (i - 1) * endColor) / (numSteps - 1));
end
    
axis equal;
axis off;
if boolSaveFigure == 1
    saveas(h1, [dataSet, '_0forces.fig'])
    saveas(h, [dataSet, '_0points.fig'])
end