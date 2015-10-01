% function  [covarianceMatrix, covarianceMatrixNormalisedSpace]  = ...
%         compute_covariance_of_sampson_estimate(...
%                   algebraicEllipseParameters, dataPts, covList)
%
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: February 2013
%
% Description: This procedure takes as input an estimated of ellipse 
%              parameters [a b c d e f] associated with the equation
%
%                        a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%              where it is assumed that the parameters were estimated
%              using the approximate maximum likelihood (Sampson distance).
%              
%               The function also takes as input the data points on which
%               the ellipse parameters were estimated, and optionally
%               a list of covariance matrices for the data points. If the
%               list of data point covariance matrices is not specified,
%               then it is assumed that the data points were corrupted
%               by isotropic homogeneous Gaussian noise, in which case
%               the noise level is estimated from the data points. 
%
%               This function return a covariance matrix for the estimated
%               algebraic parameters in the original data space. It also
%               returns the covariance matrix in a normalised coordinate
%               system because the covariance matrix in the original 
%               coordinate system usually consists of tiny numbers.
%               Both covariance matrices encode the same information,
%               but the covariance matrix in the normalised coordinate
%               system produces values that lie in a more sensible range
%               and is numerically more useful for subsequent analysis. 
%               
%
%              
%
% Parameters :  algebraicEllipseParameters -  an estimate of the algebraic
%                                             ellipse parameters based on
%                                             the Sampson cost function in 
%                                             the original (unnormalised) 
%                                             data space.

%               dataPts                   -  an nPoints x 2 matrix of 
%                                            coordinates
%
%               covList                   - a list of N 2x2 covariance 
%                                           matrices representing the
%                                           uncertainty of the coordinates
%                                           of each data point.
%                                           if this parameter is not 
%                                           specified then  default 
%                                           isotropic  (diagonal)and 
%                                           homogeneous (same noise level 
%                                           for each data point) covariance
%                                           matrices are assumed, and the
%                                           noise level is estimated from
%                                           the data points
%
%
%
% Return     : covarianceMatrix summarising the uncertainty of the 
%              algebraic parameter estimate associated with the Sampson
%              cost function (aka approximate maximum likelihood cost)
%           
%
%
% Example of Usage:
%
%
%
% Credit:
%
%
% Last Modified:
% 19/03/2014
%
function [covarianceMatrix, covarianceMatrixNormalisedSpace] = ...
                                compute_covariance_of_sampson_estimate(...
    algebraicEllipseParameters,...
    dataPts, covList)

nPts = length(dataPts);

% Check to see if the user passed in their own list of covariance matrices
if (~exist('covList', 'var'))
    % Generate a list of diagonal covariance matrices
    covList = mat2cell(repmat(eye(2),1,nPts),2,2.*(ones(1,nPts)));
    sigma_squared = estimateNoiseLevel(algebraicEllipseParameters,...
                                                 dataPts, covList);
    
    % ensure that the isotropic covariance matrices are scaled with
    % an estimate of the noise level
    for iPts = 1:nPts
     covList{iPts} = covList{iPts} * sigma_squared;
    end
end

% the algebraicEllipseParameters were computed in a hartley normalised
% coordinate system so in order to correctly characterise the uncertainty
% of the estimate, we need to know the transformation matrix T that maps
% between the original coordinate system and the hartley normalised
% coordinate system.
%
% scale and translate data points so that they lie inside a unit box
[dataPts, T] = normalize_data_isotropically(dataPts);


% transfer initialParameters to normalized coordinate system
% the formula appears in the paper Z.Szpak, W. Chojnacki and A. van den
% Hengel, "A comparison of ellipse fitting methods and implications for
% multiple view geometry", Digital Image Computing Techniques and
% Applications, Dec 2012, pp 1--8
algebraicEllipseParameters = algebraicEllipseParameters / ...
                                          norm(algebraicEllipseParameters);
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
  
algebraicEllipseParametersNormalisedSpace = E \ P34*pinv(D3)*...
                       inv(kron(T,T))'*D3*P34*E*algebraicEllipseParameters;
                     
algebraicEllipseParametersNormalisedSpace = ...
                    algebraicEllipseParametersNormalisedSpace / ...
                 norm(algebraicEllipseParametersNormalisedSpace);

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


t = algebraicEllipseParametersNormalisedSpace;
t = t / norm(t);
x = dataPts;
numberOfPoints = length(x);
M = zeros(6,6);
aml = 0;
for i = 1:numberOfPoints
    m = x(i,:);
    ux_i = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
    dux_i =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
    
    A = ux_i * ux_i';
    
    % covariance matrix of the ith data pont
    covX_i = normalised_CovList{i}; 
    
    B = dux_i * covX_i * dux_i';
    
    tAt = t' * A * t;
    tBt = t' * B * t;
    aml = aml + abs(tAt/tBt);
    M = M + A /(tBt);
end
Pt = eye(6) - (t*t')/norm(t)^2;
% compute rank-5 constrained pseudo-inverse of M
[U, D, V] = svd(M);
for i = 1:5
    D(i,i) = 1/D(i,i);
end
D(6,6) = 0;
pinvM = V*D*U';

covarianceMatrixNormalisedSpace = Pt*pinvM*Pt;

% transform covariance matrix from normalised coordinate system
% back to the original coordinate system
t = algebraicEllipseParameters / norm(algebraicEllipseParameters);                            
F = E \ P34*pinv(D3)*kron(T,T)'*D3*P34*E;                
P = eye(6) - norm(t)^2*(t*t'); 
covarianceMatrix = ...
          norm(t)^-2* P*F*covarianceMatrixNormalisedSpace*F'*P;

end

function sigma_squared = estimateNoiseLevel(algebraicEllipseParameters,...
                                                 dataPts, covList)
    t = algebraicEllipseParameters;
    t = t / norm(t);
    x = dataPts;
    numberOfPoints = length(x);
    M = zeros(6,6);
    aml = 0;
    for i = 1:numberOfPoints
        m = x(i,:);
        ux_i = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
        dux_i =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
    
        A = ux_i * ux_i';
        
        % covariance matrix of the ith data pont
        covX_i = covList{i}; 
    
        B = dux_i * covX_i * dux_i';
    
        tAt = t' * A * t;
        tBt = t' * B * t;
        aml = aml + abs(tAt/tBt);
        M = M + A /(tBt);
    end
    
    % estimate of the noise level based on average residual
    sigma_squared = aml / (numberOfPoints-5);
    
end


