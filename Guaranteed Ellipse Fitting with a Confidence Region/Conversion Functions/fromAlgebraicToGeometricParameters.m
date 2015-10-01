% function geometricEllipseParameters = fromAlgebraicToGeometricParameters(
%                                            ...algebraicEllipseParameters)
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: February 2013
%
% Description: This procedure takes as input a vector representing the 
%              algebraic parameters of an ellipse equation and produces 
%              equivalent geometric parameters so that the ellipse can be
%              expressed in parametric form.
%              
%
% Parameters : algebraicEllipseParameters  - a length-6 vector containing
%                                            an the ellipse parameters 
%                                            theta = [a b c d e f]
%                                            associated with the ellipse
%                                            equation,
%
%                            a*x^2+ b * x y + c * y^2 + d * x + e*y + f = 0
%           
%                          with the additional result that b^2 - 4 a c < 0.
%
%
%
% Return     : geometricEllipseParameters - a length-5 vector of 
%                                           geometrically meaningful 
%                                           ellipse parameters,
%
%                              majAxis: half-length of major axis
%                              minAxis: half-length of minor axis
%                              xCenter: x-coordinate of ellipse centroid
%                              yCenter: y-coordinate of ellipse centroid
%                              tilt   : rotation of ellipse axis in radians 
%                                       as measured counter-clockwise from
%                                       the postive x-axis 
%                                          
%              
%
% Example of Usage:
%
%
%
%
% Last Modified:
%
%
function geometricEllipseParameters = fromAlgebraicToGeometricParameters(...
                                                algebraicEllipseParameters)


a = algebraicEllipseParameters(1);
b = algebraicEllipseParameters(2);
c = algebraicEllipseParameters(3);
d = algebraicEllipseParameters(4);
e = algebraicEllipseParameters(5);
f = algebraicEllipseParameters(6);

 
delta = b^2 - 4*a*c;
lambdaPlus = 0.5*(a + c - (b^2 + (a - c)^2)^0.5);
lambdaMinus = 0.5*(a + c + (b^2 + (a - c)^2)^0.5);

psi = b*d*e - a*e^2 - b^2*f + c*(4*a*f - d^2);
Vplus = (psi/(lambdaPlus*delta))^0.5;
Vminus = (psi/(lambdaMinus*delta))^0.5;

% major semi-axis
axisA = max(Vplus,Vminus);
% minor semi-axis
axisB = min(Vplus,Vminus);

% determine x-coordinate of ellipse centroid
xCenter = (2*c*d - b*e)/(delta);
yCenter = (2*a*e - b*d)/(delta);

% angle between x-axis and major axis 
tau = 0;
% determine tilt of ellipse in radians
if (Vplus >= Vminus)
  if(b == 0 && a < c)
      tau = 0;
  elseif (b == 0 && a >= c)
      tau = 0.5*pi;
  elseif (b < 0 && a < c)
      tau = 0.5*acot((a - c)/b);
  elseif (b < 0 && a == c) 
      tau = pi/4;
  elseif (b < 0 && a > c)
      tau = 0.5*acot((a - c)/b) + pi/2;
  elseif (b > 0 && a < c)
      tau = 0.5*acot((a - c)/b) + pi;
  elseif (b > 0 && a == c)
      tau = pi*(3/4);
  elseif (b > 0 && a > c)
      tau = 0.5*acot((a - c)/b) + pi/2;
  end
elseif (Vplus < Vminus)
  if(b == 0 && a < c)
       tau = pi/2;
  elseif (b == 0 && a >= c)
           tau = 0;
  elseif (b < 0 && a < c)
      tau = 0.5*acot((a - c)/b) + pi/2;
  elseif (b < 0 && a == c)
       tau = pi*(3/4);
  elseif (b < 0 && a > c)
      tau = 0.5*acot((a - c)/b) + pi;
  elseif (b > 0 && a < c)
      tau = 0.5*acot((a - c)/b) + pi/2;
  elseif (b > 0 && a == c)
      tau = pi/4;
  elseif (b > 0 && a > c)
      tau = 0.5*acot((a - c)/b);
  end
end
 

% notation in paper
geometricEllipseParameters = [axisA axisB xCenter yCenter tau]';
    
end

