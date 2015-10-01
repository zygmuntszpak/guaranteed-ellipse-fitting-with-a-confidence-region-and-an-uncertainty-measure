% function [normalizedPoints T] = normalize_data_isotropically( 
%                                                     pointOnEllipseMatrix)
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: February 2013
%
% Description: This procedure takes as input a matrix of two-dimensional 
%              coordinates and normalizes the coordinates so that they
%              lie inside a unit box. 
%
% Parameters : dataPts               - an nPoints x 2 matrix of 
%                                      coordinates
%
% Return     : normalizedPts         - an nPoints x 2 matrix of 
%                                      coordinates which are constrained
%                                      to lie inside a unit box. 
%
%            : T                     - a 3x3 affine transformation matrix 
%                                      T that was used to transform the 
%                                      (homogenous coordinates) of the data 
%                                      points so that they lie inside a 
%                                      unit box.
%              
%
%
% Example of Usage:
%
%
%
% Citation: W. Chojnacki and M. Brookes, "On the Consistency of the
%           Normalized Eight-Point Algorithm", J Math Imaging Vis (2007)
%           28: 19-27
%
%
% Last Modified:
%
%

function [normalizedPts, T] = normalize_data_isotropically( ...
                                                      dataPts)

[nPoints, ~] = size(dataPts);
% homogenous representation of data points resulting in a 3 x nPoints
% matrix, where the first row contains all the x-coordinates, the second
% row contains all the y-coordinates and the last row contains the
% homogenous coordinate 1.
points  = [dataPts ones(nPoints,1)]';

meanX = mean(points(1,:));
meanY = mean(points(2,:));

% isotropic scaling factor
s = sqrt(  (1/(2*nPoints))*sum((points(1,:) - meanX).^2 + ...
                                                (points(2,:)-meanY).^2)  );
T = [s^-1    0      -s^-1*meanX;
     0       s^-1   -s^-1*meanY;
     0       0                 1];
 
 
 
 normalizedPts = T*points; 
 % remove homogenous coordinate
 normalizedPts = normalizedPts';
 normalizedPts(:,3) = [];

end

