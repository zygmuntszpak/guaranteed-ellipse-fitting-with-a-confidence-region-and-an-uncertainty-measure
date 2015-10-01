% function geometricCovarianceMatrix = ...
%           compute_covariance_of_geometric_parameters(...
%                             algebraicEllipseParameters, dataPts, covList)
%
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: February 2013
%
%              This procedure takes as input an estimated of ellipse 
%              parameters [a b c d e f] associated with the equation
%
%                        a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%              where it is assumed that the parameters were estimated
%              using the approximate maximum likelihood (Sampson distance).
%              
%              The function also takes as input the data points on which
%              the ellipse parameters were estimated, and optionally
%              a list of covariance matrices for the data points. If the
%              list of data point covariance matrices is not specified,
%              then it is assumed that the data points were corrupted
%              by isotropic homogeneous Gaussian noise, in which case
%              the noise level is estimated from the data points. 
%
%              This function return a covariance matrix for geometric
%              ellipse parameters in the original data space.
%
%              Note than when the ellipse as actually a circle the 
%              covariance of the orientation becomes undefined.
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
% Return     : covarianceMatrix summarising the uncertainty of  
%              geometric ellipse parameters associated with the Sampson
%              cost function (aka approximate maximum likelihood cost)
%
%              the diagonal of the covariance matrix containts the
%              uncertainty associated with the major axis, minor axis,
%              x-centre, y-centre and orientation (in radians)
%              respectively.
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
%
%
function geometricCovarianceMatrix = ...
      compute_covariance_of_geometric_parameters(...
                              algebraicEllipseParameters, dataPts, covList)

  if (~exist('covList', 'var'))
    [~, thetaCovarianceMatrixNormalisedSpace] = ...
                                compute_covariance_of_sampson_estimate(...
    algebraicEllipseParameters,...
    dataPts);
  else
    [~, thetaCovarianceMatrixNormalisedSpace] = ...
                                compute_covariance_of_sampson_estimate(...
    algebraicEllipseParameters,...
    dataPts, covList);
  end



% the algebraicEllipseParameters were computed in a hartley normalised
% coordinate system so in order to correctly characterise the uncertainty
% of the estimate, we need to know the transformation matrix T that maps
% between the original coordinate system and the hartley normalised
% coordinate system.
%
% scale and translate data points so that they lie inside a unit box
[~, T] = normalize_data_isotropically(dataPts);

% extract isotropic scaling factor from matrix T which we will require
% in the transformation of the geometric parameter covariance matrix
% from a normalised coordinate system to the original data space
s = T(1,1)^-1;


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

                         

                               
a = algebraicEllipseParametersNormalisedSpace(1);
b = algebraicEllipseParametersNormalisedSpace(2);
c = algebraicEllipseParametersNormalisedSpace(3);
d = algebraicEllipseParametersNormalisedSpace(4);
e = algebraicEllipseParametersNormalisedSpace(5);
f = algebraicEllipseParametersNormalisedSpace(6);  



% various computations needed to build the jacobian matrix
% of the transfromation from theta (algebraic parameters) to
% eta (geometric parameters)
delta = b^2 - 4*a*c;
lambdaPlus = 0.5*(a + c - (b^2 + (a - c)^2)^0.5);
lambdaMinus = 0.5*(a + c + (b^2 + (a - c)^2)^0.5);

psi = b*d*e - a*e^2 - b^2*f + c*(4*a*f - d^2);
Vplus = (psi/(lambdaPlus*delta))^0.5;
Vminus = (psi/(lambdaMinus*delta))^0.5;

dXcenter = derivativeXcenter(a,b,c,d,e,delta);
dYcenter = derivativeYcenter(a,b,c,d,e,delta);

dTau = derivativeTau(a,b,c);

dVplus = derivativeVplus(a,b,c,d,e,f,psi,lambdaPlus,delta);
dVminus = derivativeVminus(a,b,c,d,e,f,psi,lambdaMinus,delta);

A = max(Vplus,Vminus);
B = min(Vplus,Vminus);

if (A == Vplus)
    dA = dVplus;
    dB = dVminus;
else
    dA = dVminus;
    dB = dVplus;
end

% jacobian matrix of the transformation from theta to eta (geometric)
% parameters
etaDtheta = [dA';dB';dXcenter';dYcenter';dTau'];

% propogate uncertainty from the algebraic parameter space (theta)
% to the geometric parameter space (eta) in a normalised coordinate
% system for maximum numerical accuracy
etaCovarianceMatrixNormalisedSpace =...
                  etaDtheta*thetaCovarianceMatrixNormalisedSpace*etaDtheta';

% apply denormalisation step to determine geometric parameter covariance
% matrix in the original data space
denormalisationMatrix = diag([s,s,s,s,1]);
etaCovarianceMatrix = ...
        denormalisationMatrix*etaCovarianceMatrixNormalisedSpace*...
                                                denormalisationMatrix';


geometricCovarianceMatrix = etaCovarianceMatrix;
end

function gradXcenter = derivativeXcenter(a,b,c,d,e,delta)
gradXcenter = zeros(6,1);
gradXcenter(1) =  (4*c*(2*c*d-b*e))/delta^2;
gradXcenter(2) =  (b^2*e+4*a*c*e-4*b*c*d)/delta^2;
gradXcenter(3) =  (2*b*(b*d-2*a*e))/delta^2;
gradXcenter(4) =  (2*c)/delta;
gradXcenter(5) =   -b/delta;  
gradXcenter(6) =   0; 
end

function gradYcenter = derivativeYcenter(a,b,c,d,e,delta)
gradYcenter = zeros(6,1);
gradYcenter(1) =  (2*b*(b*e-2*c*d))/delta^2;
gradYcenter(2) =  (b^2*d+4*a*c*d-4*a*b*e)/delta^2;
gradYcenter(3) =  (4*a*(2*a*e-b*d))/delta^2;
gradYcenter(4) =  -b / delta;
gradYcenter(5) =  2*a/delta; 
gradYcenter(6) =  0; 
end

function gradTau = derivativeTau(a,b,c)
% add condition for when a==c and b = 0
gradTau = zeros(6,1);
gradTau(1) = -b/(2*(b^2+(a-c)^2)); 
gradTau(2) =  (a-c)/(2*(b^2+(a-c)^2));
gradTau(3) =  b/(2*(b^2+(a-c)^2));
gradTau(4) =  0;
gradTau(5) =  0;  
gradTau(6) =  0;   
end

function gradVplus = derivativeVplus(a,b,c,d,e,f,psi,lambdaPlus,delta)
gradVplus = zeros(6,1);
gradVplus(1) =  computeDerivativeVplusA(a,b,c,d,e,f,psi,lambdaPlus,delta);
gradVplus(2) =  computeDerivativeVplusB(a,b,c,d,e,f,psi,lambdaPlus,delta);
gradVplus(3) =  computeDerivativeVplusC(a,b,c,d,e,f,psi,lambdaPlus,delta);
gradVplus(4) =  computeDerivativeVplusD(a,b,c,d,e,f,psi,lambdaPlus,delta);
gradVplus(5) =  computeDerivativeVplusE(a,b,c,d,e,f,psi,lambdaPlus,delta);
gradVplus(6) =  computeDerivativeVplusF(a,b,c,d,e,f,psi,lambdaPlus,delta);
end

function dVplusA = computeDerivativeVplusA(a,b,c,~,e,f,psi,...
                                                         lambdaPlus,delta)
part1 = 1/(2*lambdaPlus*delta);
part2 = (psi/(lambdaPlus*delta))^(-0.5);
part3 = 4*c*f - e^2 + 4*delta^(-1)*c*psi;
part4 = (psi/(2*lambdaPlus))*(1 + ((c - a)/((a - c)^2 + b^2)^0.5));
dVplusA = part1*part2*(part3 - part4);
end

function dVplusB = computeDerivativeVplusB(a,b,c,d,e,f,psi,...
                                                          lambdaPlus,delta)
part1 = 1/(2*lambdaPlus*delta);
part2 = (psi/(lambdaPlus*delta))^(-0.5);
part3 = d*e - 2*b*f - 2*delta^(-1)*b*psi;
part4 = (b*psi)/((2*lambdaPlus)*((a - c)^2 + b^2)^0.5);
dVplusB =  part1*part2*(part3 + part4);                                                     
end

function dVplusC = computeDerivativeVplusC(a,b,c,d,~,f,psi,...
                                                        lambdaPlus,delta)
part1 = 1/(2*lambdaPlus*delta);
part2 = (psi/(lambdaPlus*delta))^(-0.5);
part3 = 4*a*f - d^2 + 4*delta^(-1)*a*psi;
part4 = (psi/(2*lambdaPlus))*(1 + ((a - c)/((a - c)^2 + b^2)^0.5));
dVplusC =part1*part2*(part3 - part4);
end

function dVplusD = computeDerivativeVplusD(~,b,c,d,e,~,psi,...
                                                         lambdaPlus,delta)
part1 = (psi/(lambdaPlus*delta))^(0.5);
part2 = (b*e - 2*c*d)/(2*psi);
dVplusD = part1*part2;                                                     
end

function dVplusE = computeDerivativeVplusE(a,b,~,d,e,~,psi,...
                                                         lambdaPlus,delta)
                                                    
part1 = (psi/(lambdaPlus*delta))^(0.5);
part2 = (b*d - 2*a*e)/(2*psi);
dVplusE = part1*part2;
end

function dVplusF = computeDerivativeVplusF(~,~,~,~,~,~,psi,...
                                                         lambdaPlus,delta)
part1 = -(2*lambdaPlus)^-1;
part2 = (psi/(lambdaPlus*delta))^(-0.5);
dVplusF = part1*part2;
end

function gradVminus = derivativeVminus(a,b,c,d,e,f,psi,lambdaMinus,delta)
gradVminus = zeros(6,1);
gradVminus(1) = computeDerivativeVminusA(a,b,c,d,e,f,psi,lambdaMinus,delta);
gradVminus(2) = computeDerivativeVminusB(a,b,c,d,e,f,psi,lambdaMinus,delta);
gradVminus(3) = computeDerivativeVminusC(a,b,c,d,e,f,psi,lambdaMinus,delta);
gradVminus(4) = computeDerivativeVminusD(a,b,c,d,e,f,psi,lambdaMinus,delta);
gradVminus(5) = computeDerivativeVminusE(a,b,c,d,e,f,psi,lambdaMinus,delta);
gradVminus(6) = computeDerivativeVminusF(a,b,c,d,e,f,psi,lambdaMinus,delta);
end

function dVminusA = computeDerivativeVminusA(a,b,c,~,e,f,psi,...
                                                         lambdaMinus,delta)
part1 = 1/(2*lambdaMinus*delta);
part2 = (psi/(lambdaMinus*delta))^(-0.5);
part3 = 4* c*f - e^2 + 4*delta^(-1)*c*psi;
part4 = (psi/(2*lambdaMinus))*(1 - ((c - a)/((a - c)^2 + b^2)^0.5));
dVminusA = part1*part2*(part3 - part4);
end

function dVminusB = computeDerivativeVminusB(a,b,c,d,e,f,psi,...
                                                         lambdaMinus,delta)
part1 = 1/(2*lambdaMinus*delta);
part2 = (psi/(lambdaMinus*delta))^(-0.5);
part3 = d*e - 2*b*f - 2*delta^(-1)*b*psi;
part4 = (b*psi)/((2*lambdaMinus)*((a - c)^2 + b^2)^0.5);
dVminusB = part1*part2*(part3 - part4);
end

function dVminusC = computeDerivativeVminusC(a,b,c,d,~,f,psi,...
                                                         lambdaMinus,delta)
part1 = 1/(2*lambdaMinus*delta);
part2 = (psi/(lambdaMinus*delta))^(-0.5);
part3 = 4*a*f - d^2 + 4*delta^(-1)*a*psi;
part4 = (psi/(2*lambdaMinus))*(1 - ((a - c)/((a - c)^2 + b^2)^0.5));
dVminusC =  part1*part2*(part3 - part4);
end

function dVminusD = computeDerivativeVminusD(~,b,c,d,e,~,psi,...
                                                         lambdaMinus,delta)
part1 = (psi/(lambdaMinus*delta))^(0.5);
part2 = (b*e - 2*c*d)/(2*psi);
dVminusD =  part1*part2;
end

function dVminusE = computeDerivativeVminusE(a,b,~,d,e,~,psi,...
                                                         lambdaMinus,delta)
part1 = (psi/(lambdaMinus*delta))^(0.5);
part2 = (b*d - 2*a*e)/(2*psi);
dVminusE = part1*part2;
end

function dVminusF = computeDerivativeVminusF(~,~,~,~,~,~,psi,...
                                                         lambdaMinus,delta)
part1 = -(2*lambdaMinus)^-1;
part2 = (psi/(lambdaMinus*delta))^(-0.5);
dVminusF = part1*part2;
end