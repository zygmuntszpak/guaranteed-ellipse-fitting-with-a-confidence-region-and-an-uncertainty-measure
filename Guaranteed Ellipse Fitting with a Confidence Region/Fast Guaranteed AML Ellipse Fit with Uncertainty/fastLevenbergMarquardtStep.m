%   Function: fastLevenbergMarquardtStep
%
%   This function is used in the main loop of guaranteedEllipseFit in the
%   process of minimizing an approximate maximum likelihood cost function 
%   of an ellipse fit to data.  It computes an update for the parameters
%   representing the ellipse, using the method of Levenberg-Marquardt for
%   non-linear optimisation. 
%   See: http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
%
%   However, unlike the traditional LevenbergMarquardt step, we do not
%   add a multiple of the identity matrix to the approximate Hessian,
%   but instead a different positive semi-definite matrix. Our choice
%   particular choice of the different matrix corresponds to the 
%   gradient descent direction in the theta coordinate system, 
%   transformed to the eta coordinate system. We found empirically
%   that taking steps according to the theta coordinate system
%   instead of the eta coordinate system lead to faster convergence.
%
%   Parameters:
%
%      struct     - a data structure containing various parameters
%                   needed for the optimisation process. 	
%
%   Returns: 
%
%     the same data structure 'struct', except that relevant fields have
%     been updated
%
%   See Also: 
%
%    fastGuaranteedEllipseFit
%
%  Zygmunt L. Szpak (c) 2014
%  Last modified 18/3/2014 
function struct = fastLevenbergMarquardtStep( struct,rho )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract variables from data structure                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    jacobian_matrix = struct.jacobian_matrix;
    r = struct.r;
    lambda = struct.lambda;
    delta = struct.delta(struct.k);
    damping_multiplier = struct.damping_multiplier;
    damping_divisor = struct.damping_divisor;
    current_cost = struct.cost(struct.k);
    data_points = struct.data_points;
    covList = struct.covList;
    numberOfPoints = struct.numberOfPoints;
    H = struct.H;    
    jlp = struct.jacob_latentParameters;        
    eta = struct.eta(:,struct.k);
    
    % convert latent variables into length-6 vector (called t) representing
    % the equation of an ellipse
    t = [1, 2*eta(1), eta(1)^2 +  abs(eta(2))^rho, eta(3),eta(4), eta(5)]';
    % we impose unit norm constraint on theta
    t = t /norm(t);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute two potential updates for theta based on different         %
    % weightings of the identity matrix.                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    jacob = [jacobian_matrix]'*r;           
    DMP = (jlp'*jlp)*lambda;
    update_a = - (H+DMP)\jacob;     
       
    % In a similar fashion, the second potential search direction 
    % is computed
    
    DMP = (jlp'*jlp)*lambda/damping_divisor;
    update_b = - (H+DMP)\jacob;     
   
    % the potential new parameters are then 
    eta_potential_a = eta + update_a;
    eta_potential_b = eta + update_b;
    
    % we need to convert from eta to theta and impose unit norm constraint
    t_potential_a =  [1, 2*eta_potential_a(1), eta_potential_a(1)^2 +  ...
                      abs(eta_potential_a(2))^rho, eta_potential_a(3), ... 
                                 eta_potential_a(4), eta_potential_a(5)]';
    t_potential_a = t_potential_a/norm(t_potential_a);
    
    t_potential_b =   [1, 2*eta_potential_b(1), eta_potential_b(1)^2 + ...
                      abs(eta_potential_b(2))^rho, eta_potential_b(3), ... 
                                  eta_potential_b(4), eta_potential_b(5)]'; 
    t_potential_b = t_potential_b/norm(t_potential_b);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new residuals and costs based on these updates             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % residuals computed on data points
    cost_a = 0;
    cost_b = 0;
    for i = 1:numberOfPoints
        m = data_points(:,i);
        % transformed data point
        ux_i = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
        % derivative of transformed data point
        dux_i =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
        
        % outer product
        A = ux_i * ux_i';
        
        % covariance matrix of the ith data pont
        covX_i = covList{i}; 
        
        B = dux_i * covX_i * dux_i';
    
        t_aBt_a = t_potential_a' * B * t_potential_a;
        t_aAt_a = t_potential_a'  * A * t_potential_a;
        
        t_bBt_b = t_potential_b' * B * t_potential_b;
        t_bAt_b = t_potential_b'  * A * t_potential_b;
           
        % AML cost for i'th data point
        cost_a = cost_a +  abs(t_aAt_a/t_aBt_a) ;
        cost_b = cost_b +  abs(t_bAt_b/t_bBt_b) ;     
               
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % determine appropriate damping and if possible select an update     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (cost_a >= current_cost && cost_b >= current_cost)
        % neither update reduced the cost
        struct.eta_updated = false;
        % no change in the cost
        struct.cost(struct.k+1) = current_cost;
        % no change in parameters
        struct.eta(:,struct.k+1) = eta;
        struct.t(:,struct.k+1) = t;
        % no changes in step direction
        struct.delta(:,struct.k+1) = delta;
        % next iteration add more Identity matrix
        struct.lambda = lambda * damping_multiplier; 
    elseif (cost_b < current_cost)
        % update 'b' reduced the cost function
        struct.eta_updated = true;
        % store the new cost
        struct.cost(struct.k+1) = cost_b;
        % choose update 'b'
        struct.eta(:,struct.k+1) = eta_potential_b;
        struct.t(:,struct.k+1) = t_potential_b;
        % store the step direction
        struct.delta(:,struct.k+1) = update_b';
        % next iteration add less Identity matrix
        struct.lambda = lambda / damping_divisor;
    else
        % update 'a' reduced the cost function
        struct.eta_updated = true;
        % store the new cost
        struct.cost(struct.k+1) = cost_a;
        % choose update 'a'
        struct.eta(:,struct.k+1) = eta_potential_a;
        struct.t(:,struct.k+1) = t_potential_a;
        % store the step direction
        struct.delta(:,struct.k+1) = update_a';
        % keep the same damping for the next iteration
        struct.lambda = lambda;
    end

    % return a data structure containing all the updates
    struct;
end

