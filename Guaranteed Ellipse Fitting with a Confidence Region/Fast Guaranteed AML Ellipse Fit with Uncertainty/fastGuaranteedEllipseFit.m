%   Function: fastGuaranteedEllipseFit
%
%   This function implements the ellipse fitting algorithm described in
%   Z.Szpak, W. Chojnacki and A. van den Hengel
%   "Guaranteed Ellipse Fitting with an Uncertainty Measure for Centre, 
%    Axes, and Orientation"
%
%
%   Parameters:
%
%      latentParameters    - an initial seed for latent parameters
%                            [p q r s t] which through a transformation
%                            are related to parameters  [a b c d e f] 
%                            associated with the conic equation  
%                            
%                             a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%      dataPts             - a 2xN matrix where N is the number of data
%                            points
%
%      covList             - a list of N 2x2 covariance matrices 
%                            representing the uncertainty of the 
%                            coordinates of each data point.
%
%
%   Returns: 
%
%     a length-6 vector [a b c d e f] representing the parameters of the
%     equation
%    
%     a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%     with the additional result that b^2 - 4 a c < 0.
%
%   See Also: 
%
%    compute_guaranteedellipse_estimates
%    levenbergMarquardtStep
%    lineSearchStep
%
%  Zygmunt L. Szpak (c) 2014
%  Last modified 18/3/2014 
function [theta, iterations] = fastGuaranteedEllipseFit(latentParameters,...
                                                           dataPts,covList)

eta = latentParameters;
% convert latent variables into length-6 vector (called t) representing
% the equation of an ellipse
t = [1, 2*eta(1), eta(1)^2 +  abs(eta(2))^2, eta(3),eta(4), eta(5)]';
t = t / norm(t); 

% various variable initialisations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primary loop variable
keep_going = true;
% in some case a LevenbergMarquardtStep does not decrease the cost
% function and so the parameters (eta) are not updated
struct.eta_updated = false;
% damping parameter in LevenbergMarquadtStep
struct.lambda = 0.01;
% loop counter (matlab arrays start at index 1, not index 0)
struct.k = 1;
% used to modify the tradeoff between gradient descent and hessian based
% descent in LevenbergMarquadtStep
struct.damping_multiplier = 15;
% used to modify the tradeoff between gradient descent and hessian based
% descent in LevenbergMarquadtStep
struct.damping_divisor = 1.2;
% number of data points
struct.numberOfPoints = length(dataPts);
% data points that we are going to fit an ellipse to
struct.data_points = dataPts;
% a list of 2x2 covariance matrices representing the uncertainty
% in the coordinates of the data points
struct.covList = covList;

% various parameters that determine stopping criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum loop iterations 
maxIter = 200; 
% step-size tolerance
struct.tolDelta = 1e-7;
% cost tolerance
struct.tolCost = 1e-7;
% parameter tolerance
struct.tolEta = 1e-7;
% gradient tolerance
struct.tolGrad = 1e-7;
% barrier tolerance (prevent ellipse from converging on parabola)
struct.tolBar = 15.5;
% minimum allowable magnitude of conic determinant (prevent ellipse from 
% convering on degenerate parabola (eg. two parallel lines) 
struct.tolDet = 1e-5;

Fprim = [0 0 2; 0 -1 0; 2 0 0];
F = [Fprim zeros(3,3); zeros(3,3) zeros(3,3)];
I = eye(6,6);

% various initial memory allocations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocate space for cost of each iteration
struct.cost = zeros(1,maxIter);
% allocate space for the latent parameters of each iteration
struct.eta = zeros(5,maxIter);
% and for the parameters representing the ellipse equation
struct.t = zeros(6,maxIter);
% allocate space for the parameter direction of each iteration
struct.delta = zeros(5,maxIter);
% make parameter vector a unit norm vector for numerical stability
%t = t / norm(t);
% store the parameters associated with the first iteration
struct.t(:,struct.k) = t;
struct.eta(:,struct.k) = eta;
% start with some random search direction (here we choose all 1's)
% we can initialise with anything we want, so long as the norm of the
% vector is not smaller than tolDeta. The initial search direction 
% is not used in any way in the algorithm. 
struct.delta(:,struct.k) = ones(5,1);
% main estimation loop
 while (keep_going && struct.k < maxIter)
    
    % allocate space for residuals 
    struct.r = zeros(struct.numberOfPoints,1);
    % allocate space for the jacobian matrix based on AML component
    struct.jacobian_matrix = zeros(struct.numberOfPoints,5);    
    % grab the current latent parameter estimates
    eta = struct.eta(:,struct.k);
    % convert latent variables into length-6 vector (called t) representing
    % the equation of an ellipse
    t = [1, 2*eta(1), eta(1)^2 +  abs(eta(2))^2, eta(3),eta(4), eta(5)]';
 
    % jacobian matrix of the transformation from eta to theta parameters
    jacob_latentParameters = ...
        [0                     0                          0 0 0;
         2                     0                          0 0 0;
         2*eta(1) 2*abs(eta(2))^(2-1)*sign(eta(2))  0 0 0;
         0                     0                          1 0 0;
         0                     0                          0 1 0;
         0                     0                          0 0 1];
    
    % we impose the additional constraint that theta will be unit norm
    % so we need to modify the jacobian matrix accordingly
    Pt = eye(6) - ((t*t')/(norm(t,2)^2));
    jacob_latentParameters = (1/norm(t,2))*Pt*jacob_latentParameters;
    % unit norm constraint
    t = t / norm(t); 
    
    % residuals computed on data points
    for i = 1:struct.numberOfPoints
        m = dataPts(:,i);
        % transformed data point
        ux_i = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
        % derivative of transformed data point
        dux_i =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
        
        % outer product
        A = ux_i * ux_i';
        
        % covariance matrix of the ith data pont
        covX_i = covList{i}; 
        
        B = dux_i * covX_i * dux_i';
        
        tBt = t' * B * t;
        tAt = t' * A * t;
        
        % AML cost for i'th data point
        struct.r(i) = sqrt(abs(tAt/tBt));
        
        % derivative AML component
        M = (A / tBt);
        Xbits = B * ((tAt) / (tBt^2));
        X = M - Xbits;
                                
        % gradient for AML cost function (row vector)
        grad = ((X*t) / sqrt((abs(tAt/tBt)+eps)))';
        % build up jacobian matrix
        struct.jacobian_matrix(i,:) = grad * jacob_latentParameters; 
    end            
    
    % approximate Hessian matrix
    struct.H =  struct.jacobian_matrix'* struct.jacobian_matrix;
       
    % sum of squares cost for the current iteration
    struct.cost(struct.k) = (struct.r'*struct.r);
    
    struct.jacob_latentParameters =  jacob_latentParameters;
      
    % use LevenbergMarquadt step to update parameters
     struct = fastLevenbergMarquardtStep(struct,2);
     
    % Preparations for various stopping criteria tests
    

    % convert latent variables into length-6 vector (called t) representing
    % the equation of an ellipse
    eta = struct.eta(:,struct.k+1);
    t = [1, 2*eta(1), eta(1)^2 +  abs(eta(2))^2, eta(3),eta(4), eta(5)]';    
    t = t / norm(t);  
    
    % First criterion checks to see if discriminant approaches zero by using
    % a barrier 
    tIt = t' * I * t;
    tFt = t' * F * t;
    barrier = (tIt/tFt);
    
    % Second criterion checks to see if the determinant of conic approaches
    % zero
    M = [t(1) , t(2)/2 ,t(4)/2 ; 
        t(2)/2 , t(3) ,t(5)/2 ; 
        t(4)/2 ,t(5)/2 , t(6)];  
    DeterminantConic = det(M);
   
    % Check for various stopping criteria to end the main loop
    if (min(norm(struct.eta(:,struct.k+1)-struct.eta(:,struct.k)), ...
            norm(struct.eta(:,struct.k+1)+struct.eta(:,struct.k))) < ...
                                      struct.tolEta && struct.eta_updated)
        keep_going = false;
    elseif (abs(struct.cost(struct.k) - struct.cost(struct.k+1)) ...
                                    < struct.tolCost && struct.eta_updated)
        keep_going = false;
    elseif (norm(struct.delta(:,struct.k+1)) <  ...
                                     struct.tolDelta && struct.eta_updated)
        keep_going = false;
    elseif   norm(grad) < struct.tolGrad;
        keep_going = false;
    elseif (log(barrier) > struct.tolBar || ...
                                    abs(DeterminantConic) < struct.tolDet)
        keep_going = false;    
    end
    
    struct.k = struct.k + 1;    
 end

 iterations = struct.k;
 theta = struct.t(:,struct.k);
 theta = theta / norm(theta);

end

