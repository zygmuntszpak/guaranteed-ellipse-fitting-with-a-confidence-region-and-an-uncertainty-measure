
% function estimatedParameters = fast_guaranteed_ellipse_estimate(...
%                                        dataPts, covList)
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: March 2014
%
% Description: This procedure takes as input a matrix of two-dimensional 
%              coordinates and estimates a best fit ellipse using the
%              sampson distance between a data point and the ellipse
%              equations as the error measure. The Sampson distance is
%              an excellent approximation to the orthogonal distance for
%              small noise levels. The Sampson distance is often also
%              referred to as the approximate maximum likelihood (AML).
%              The user can specify a list of covariance matrices for the
%              data points. If the user does not specify a list of
%              covariance matrices then isotropic homogeneous Gaussian
%              noise is assumed. 
%
% Parameters : initialParameters    - initial parameters use to seed the
%                                     iterative estimation process
%
%              dataPts                - an nPoints x 2 matrix of 
%                                      coordinates
%
%              covList               - a list of N 2x2 covariance matrices 
%                                      representing the uncertainty of the 
%                                      coordinates of each data point.
%                                      if this parameter is not specified 
%                                      then  default isotropic  (diagonal)
%                                      and homogeneous (same noise level 
%                                      for each data point) covariance
%                                      matrices are assumed.
%
%
% Return     : a length-6 vector containing an estimate of the ellipse
%              parameters theta = [a b c d e f] associated with the ellipse
%              equation 
%
%                   a*x^2+ b * x y + c * y^2 + d * x + e*y + f = 0
%           
%              with the additional result that b^2 - 4 a c < 0.
%
% Example of Usage:
%
%
%
%
% Last Modified:
% 18/03/2014
%
function [estimatedParameters, iterations] = ...
    fast_guaranteed_ellipse_estimate(dataPts, covList)

nPts = length(dataPts);

% Check to see if the user passed in their own list of covariance matrices
if (~exist('covList', 'var'))
   % Generate a list of diagonal covariance matrices   
   covList = mat2cell(repmat(eye(2),1,nPts),2,2.*(ones(1,nPts)));
end

% estimate an initial ellipse using the direct ellipse fit method
initialEllipseParameters   = compute_directellipse_estimates(dataPts);

% scale and translate data points so that they lie inside a unit box
[normalizedPoints, T] = normalize_data_isotropically(dataPts);


% transfer initialParameters to normalized coordinate system
% the formula appears in the paper Z.Szpak, W. Chojnacki and A. van den
% Hengel, "A comparison of ellipse fitting methods and implications for
% multiple view geometry", Digital Image Computing Techniques and
% Applications, Dec 2012, pp 1--8
initialEllipseParameters = initialEllipseParameters / ...
                                            norm(initialEllipseParameters);
E = diag([1,2^-1,1,2^-1,2^-1,1]);
% permutation matrix for interchanging 3rd and 4th
% entries of a length-6 vector
P34 = kron(diag([0,1,0]), [0 1;1 0]) + kron(diag([1,0,1]), [1 0; 0 1]);
% 9 x 6 duplication matrix
D3 = [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 1 0 0 0 0;
      0 0 0 1 0 0   
      0 0 0 0 1 0;
      0 0 1 0 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1];    
  
initialEllipseParametersNormalizedSpace = E \ P34*pinv(D3)*...
                         inv(kron(T,T))'*D3*P34*E*initialEllipseParameters;
                     
initialEllipseParametersNormalizedSpace = ...
                    initialEllipseParametersNormalizedSpace / ...
                             norm(initialEllipseParametersNormalizedSpace);

                         
% Becase the data points are now in a new normalised coordinate system,
% the data covariance matrices also need to be tranformed into the 
% new normalised coordinate system. The transformation of the covariance
% matrices into the new coordinate system can be achieved by embedding the
% covariance matrices in a 3x3 matrix (by padding the 2x2 covariance
% matrices by zeros) and by  multiply the covariance matrices by the 
% matrix T from the left and T' from the right. 
normalised_CovList = cell(1,nPts);
for iPts = 1 : nPts
    covX_i = zeros(3,3);
    covX_i(1:2,1:2) = covList{iPts};
    covX_i = T*covX_i*T';    
    % the upper-left 2x2 matrix now represents the covariance of the 
    % coordinates of the data point in the normalised coordinate system
    normalised_CovList{iPts} = covX_i(1:2,1:2);
end
                        
                         
                         
                         
% To guarantee an ellipse we utilise a special parameterisation which
% by definition excludes the possiblity of a hyperbola. In theory 
% a parabola could be estimated, but this is very unlikely because
% an equality constraint is difficult to satisfy when there is noisy data.
% As an extra guard to ensure that a parabolic fit is not possible we
% terminate our algorithm when the discriminant of the conic equation
% approaches zero. 


% convert our original parameterisation to one that excludes hyperbolas
% NB, it is assumed that the initialParameters that were passed into the
% function do not represent a hyperbola or parabola.
para = initialEllipseParametersNormalizedSpace;
% p = para(2)/(2*para(1));
% q = 1 / sqrt(para(3)/para(1) - (para(2)/(2*para(1)))^2);
% r = para(4) / para(1);
% s = para(5) / para(1);
% t = para(6) / para(1);
p = para(2)/(2*para(1));
q = (para(3)/para(1) - (para(2)/(2*para(1)))^2)^(1/2);
r = para(4) / para(1);
s = para(5) / para(1);
t = para(6) / para(1);

latentParameters  = [p q r s t];


[ellipseParametersFinal, iterations] = ...
          fastGuaranteedEllipseFit(latentParameters,...
                                    normalizedPoints',normalised_CovList);
    
           
ellipseParametersFinal = ellipseParametersFinal /...
                                              norm(ellipseParametersFinal);
                                     
% convert final ellipse parameters back to the original coordinate system
estimatedParameters =  E \ P34*pinv(D3)*...
                    kron(T,T)'*D3*P34*E*ellipseParametersFinal;


estimatedParameters = estimatedParameters / norm(estimatedParameters);
estimatedParameters = estimatedParameters / sign(estimatedParameters(end));

   
end

 





