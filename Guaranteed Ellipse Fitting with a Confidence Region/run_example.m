%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the guaranteed ellipse fitting
% method described in the paper:
%
% Z. L. Szpak, W. Chojnacki, and A. van den Hengel. 
% Guaranteed ellipse fitting with a confidence region and an uncertainty
% measure for centre, axes, and orientation. 
% J. Math. Imaging Vision, 2015. 
% http://dx.doi.org/10.1007/s10851-014-0536-x
%
% The general equation of a conic is
% a m_1^2 + b m_1m_2 + c m_2^2 + d m_1 + e m_2 + f = 0
%
% and when the discriminant b^2 - 4 a c < 0, then the conic
% is an ellipse. Through a particular choice of algebraic 
% parameters our estimation method guarantees an ellipse fit.
%
% However, the parametrisation does not guarantee a non-degenerate
% ellipse. A conic is non-degenerate, if the determinant of the 
% conic 
%
%     | a   b/2 d/2 |
% D = |b/2  c   e/2 |
%     |d/2  e/2 f   |
%
% is non-zero. When D = 0 the conic is degenerate. Degenerate ellipse are,
% for example, two parallel lines or a double line. See 
% http://en.wikipedia.org/wiki/Degenerate_conic
% for more details. 
%
% To prevent degenerate ellipses, our algorithm monitors 
% the determinant D and terminates if it falls below a threshold.
% Refer to our paper for more details.
%
% Degenerate ellipses usually only occur when points are sampled
% from a very small faction of an ellipse perimiter. 
%
% Our method can take as input a 2 x 2 covariance matrix for each data
% point. If data point covariance matrices are not provided, then the 
% algorithm assumes homogeneous Gaussian noise and estimates the level of
% noise from the data. 
%
% The algorithm produces the algebraic parameters and corresponding 
% covariance matrix, as well as geometric ellipse parameters and 
% corresponding covariance matrix. 
% 
% Additionally, it draws a 95% confidence region which gives a visual
% representation of the uncertainty in the estimate. Note however, that
% when points are sampled from a small portion of the ellipse (a very
% ill-posed problem), then certain assumptions underlying the generation
% of the confidence region are violated and the confidence region may
% no longer give the correct coverage. 
%
% The ellipse estimation and covariance matrix estimation itself is really 
% fast. However, drawing the confidence region in Matlab takes a long time,
% primarily because I don't know how to vectorise and parallelise the code
% properly. There are also occasional issues with how I overlay the
% confidence region over the existing plot. The confidence region
% computation is correct, but I struggled to render and overlay it
% over the plot in MATLAB and resorted to some workaround.
%
% If you have any suggestions on how to speed up the rendering of the
% confidence region I'd love to hear from you. 
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mySeed = 10; 
rand( 'seed', mySeed );
randn( 'seed', mySeed );

% noisy data points sampled from an ellipse
data_points = [554.4255448860487, 252.25695182609456;
    563.2198518844277, 266.393133656879;
    529.4069191555953,  293.8378723522843;
    540.4542393161609, 305.7810697735044;
    529.979807146737, 316.8126082331031;
    504.904046118169, 337.4825804687354;
    481.8704189162051, 364.2838371680541;
    444.52860559075106, 363.1246963316071;
    399.2429477714536, 377.5681669353429;
    389.3711948977463, 385.5186736607949;
    349.4010115271496,386.2628649432312;
    296.5840933162226, 399.73317608469233;
    255.11250601488572, 395.1411289783602;
    231.89205825187307, 381.61074589981473;
    182.90013818617746,384.047826837435;
    144.59336259440926,391.4128662876713;
    97.60448452233265, 391.35548699175257;
    76.28734052405164,369.8323364160187;
    25.345944215923343,367.78769391856787;
    -2.5177280792243053,339.56485858969074;
    -4.3254400982999766,305.06389938514843;
    -29.074863473090797,319.05942892395166;
    -48.73339406796903, 288.8232699112716;
    -85.21858392682377,262.78052469826156;
    -77.7559980456105,255.63059999535096;
    -63.96460001891381,240.5558431574499;
    -47.50679886929984, 221.79850196104366;
    -67.29391840704433, 209.6724162607794;
    -23.86911438928088, 176.9721535710241;
    -17.624006835223113, 174.24367926036066;
    31.38415275933525, 142.62191100176753;
    46.11460906558086, 128.39144744279602;
    65.38963556066517,126.90848784795256;
    116.43356672461559,125.50674998521659;
    135.24857840307698, 109.847540360864;
    212.24483407970993, 105.53786425855019;
    227.26351500381813, 101.51465346090531;
    255.18532643477042, 107.83975623109346;
    306.56565641797596, 100.03580557868527;
    355.07952555897765, 143.5543623673256;
    376.5664452416278, 116.9383797341007;
    416.266045280784, 102.4712728975143;
    451.06342651230636, 144.36320327683816;
    475.96962498992485,  153.3519382363760;
    507.98479900085863, 186.347765271335;
    522.3993024410767, 181.1256514766796;
    544.0355186255582, 222.30892095626606;
    561.1263519611979, 205.47120976368097;
    581.8305025046524,241.94356048025762;
    546.5064585867074, 232.51223747541997];


%%%%%%%%%%%%%%%%%%%% Example with ALL data points %%%%%%%%%%%%%%%%%%%%%%%%
% An example of fitting to all the data points
fprintf('**************************************************************\n')
fprintf('* Example with ALL data points assuming homogeneous Gaussian *\n')
fprintf('* noise with the noise level automatically estimated from    *\n')
fprintf('* data points                                                *\n')
fprintf('**************************************************************\n')

fprintf('Algebraic ellipse parameters of direct ellipse fit: \n')
[theta_dir]  = compute_directellipse_estimates(data_points)
fprintf('Algebraic ellipse parameters of our method: \n')
[theta_fastguaranteed] = fast_guaranteed_ellipse_estimate(data_points)


fprintf('Geometric ellipse parameters \n')
fprintf('(majAxis, minAxis, xCenter,yCenter, orientation (radians)): \n')
geometricEllipseParameters = ...
            fromAlgebraicToGeometricParameters(theta_fastguaranteed)

fprintf('Covariance matrix of geometric parameters: \n')
geoCov =  compute_covariance_of_geometric_parameters(...
                               theta_fastguaranteed, data_points)
 
 fprintf('Standard deviation of geometric parameters: \n')                          
 stds = sqrt(diag(geoCov)) 
 
 
[S, thetaCovarianceMatrixNormalisedSpace] = ...
                                 compute_covariance_of_sampson_estimate(...
                     theta_fastguaranteed, data_points);
   
% plot the data points
x = data_points;
n = length(x);
figure('Color',[1 1 1])
plot(x(:,1),x(:,2),'k.');
hold on

% determine data range
% minX = min(min(x(:,1))) - 20;
% minY = min(min(x(:,2))) - 20;
% maxX = max(max(x(:,1))) + 20;
% maxY = max(max(x(:,2)))  + 20;

% the particular way we plot confidence regions in matlab
% requires that we work with square axes
minX = -100;
maxX = 600;
minY = -100;
maxY = 600;

% plot the direct ellipse fit
a = theta_dir(1); b = theta_dir(2); c = theta_dir(3);
d = theta_dir(4); e = theta_dir(5); f = theta_dir(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
set(h, 'Color', [0 0 0]);
set(h,'LineWidth',1.5);  
axis([minX maxX minY maxY]);

% plot the guaranteed ellipse fit
a = theta_fastguaranteed(1); b = theta_fastguaranteed(2);
c = theta_fastguaranteed(3); d = theta_fastguaranteed(4); 
e = theta_fastguaranteed(5); f = theta_fastguaranteed(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);
set(h, 'Color', [1 0 0]);
set(h,'LineWidth',1.5);  
title('All data and homogeneous Gaussian noise')
legend('DATA POINT','DIRECT ELLIPSE FIT','FAST GUARANTEED ELLIPSE FIT',...
                                                'Location','SouthOutside');
                                            

xres=700;
yres=700;
xv = linspace(minX, maxX, xres);
yv = linspace(minY, maxY, yres);

% scale and translate data points so that they lie inside a unit box
% in this instance we are actually just after the transformation matrix T
[~, T] = normalize_data_isotropically(data_points);

% convert parameters to normalised space

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
        
theta_guaranteed_covweightedNormalised = E \ P34*pinv(D3)*...
            inv(kron(T,T))'*D3*P34*E*theta_fastguaranteed;


% now we can transform our grid points to normalised space
pts = [xv; yv; ones(1,xres)];
pts_normalised_space = T*pts;
xv_normalised = pts_normalised_space(1,:);
yv_normalised = pts_normalised_space(2,:);
Z = compute_confidence_band(xv_normalised,yv_normalised,...
                         theta_guaranteed_covweightedNormalised, ...
                         thetaCovarianceMatrixNormalisedSpace,11.07);

% this next step is a matlab hack to try to overlay the confidence
% region over the existing plot. I don't know of a better way to do
% this, and sometime the code below doesn't work....
% it relies on knowledge of the different ellipse colors
hold on
[C,h] = contourf(xv,yv,Z,[0,1],'LineColor','none','LineWidth',2);
% cmap = [[234,228,241]./255; 1 1 1];
cmap = [[237,221,221]./255; 1 1 1];
colormap(cmap)

fastGuaranteedPlot = findobj(gca, 'Color', [1 0 0]);
uistack(fastGuaranteedPlot, 'top')
%truthAndDataPlot = findobj(gca, 'Color', [0 1 0]);
directEllipsePlot = findobj(gca, 'Color', [0 0 0]);
uistack(directEllipsePlot, 'top') 
dataPointPlot = findobj(gca, 'Color', [0 0 0]);
uistack(dataPointPlot, 'top')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%% Example with only a portion of data points %%%%%%%%%%%
fprintf('**************************************************************\n')
fprintf('* Example with a only a portion of data points assuming      *\n')
fprintf('* homogeneous Gaussian noise with the noise level            *\n')
fprintf('* automatically estimated from the data points               *\n')
fprintf('**************************************************************\n')
%data_points_portion = data_points(:,1:1:(end/2));
data_points_portion = data_points(1:1:(end/2),:);
% An example of fitting to a portion of the data points
fprintf('Algebraic ellipse parameters of direct ellipse fit: \n')
[theta_dir]  = compute_directellipse_estimates(data_points_portion)
fprintf('Algebraic ellipse parameters of our method: \n')
[theta_fastguaranteed] = ...
                fast_guaranteed_ellipse_estimate(data_points_portion);
            
fprintf('Geometric ellipse parameters \n')
fprintf('(majAxis, minAxis, xCenter,yCenter, orientation (radians)): \n') 
geometricEllipseParameters = ...
            fromAlgebraicToGeometricParameters(theta_fastguaranteed)            

fprintf('Covariance matrix of geometric parameters: \n')
geoCov = compute_covariance_of_geometric_parameters(...
                              theta_fastguaranteed, data_points_portion)

fprintf('Standard deviation of geometric parameters: \n')                          
stds = sqrt(diag(geoCov))      

[S, thetaCovarianceMatrixNormalisedSpace] = ...
                                compute_covariance_of_sampson_estimate(...
                    theta_fastguaranteed, data_points_portion);
                


% plot the data points
x = data_points;
xportion = data_points_portion;
n = length(x);
figure('Color',[1 1 1])
%plot(x(:,1),x(:,2),'b.');
plot(xportion(:,1),xportion(:,2),'k.');
hold on

% determine data range
% minX = min(min(x(:,1))) - 20;
% minY = min(min(x(:,2))) - 20;
% maxX = max(max(x(:,1))) + 20;
% maxY = max(max(x(:,2)))  + 20;

% the particular way we plot confidence regions in matlab
% requires that we work with square axes
minX = -100;
maxX = 600;
minY = -100;
maxY = 600;

% plot the direct ellipse fit
a = theta_dir(1); b = theta_dir(2); c = theta_dir(3);
d = theta_dir(4); e = theta_dir(5); f = theta_dir(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
set(h, 'Color', [0 0 0]);
set(h,'LineWidth',1.5); 
axis([minX maxX minY maxY]);


% plot the guaranteed ellipse fit
a = theta_fastguaranteed(1); b = theta_fastguaranteed(2); 
c = theta_fastguaranteed(3); d = theta_fastguaranteed(4);
e = theta_fastguaranteed(5); f = theta_fastguaranteed(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);
set(h, 'Color', [1 0 0]);
set(h,'LineWidth',1.5); 
title('Partial data and homogeneous Gaussian noise')
legend('DATA POINT','DIRECT ELLIPSE FIT','FAST GUARANTEED ELLIPSE FIT',...
                                                'Location','SouthOutside');
                                            
                                            
xres=700;
yres=700;
xv = linspace(minX, maxX, xres);
yv = linspace(minY, maxY, yres);

% scale and translate data points so that they lie inside a unit box
% in this instance we are actually just after the transformation matrix T
[~, T] = normalize_data_isotropically(data_points);

% convert parameters to normalised space

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
        
theta_guaranteed_covweightedNormalised = E \ P34*pinv(D3)*...
            inv(kron(T,T))'*D3*P34*E*theta_fastguaranteed;


% now we can transform our grid points to normalised space
pts = [xv; yv; ones(1,xres)];
pts_normalised_space = T*pts;
xv_normalised = pts_normalised_space(1,:);
yv_normalised = pts_normalised_space(2,:);
Z = compute_confidence_band(xv_normalised,yv_normalised,...
                         theta_guaranteed_covweightedNormalised, ...
                         thetaCovarianceMatrixNormalisedSpace,11.07);

% this next step is a matlab hack to try to overlay the confidence
% region over the existing plot. I don't know of a better way to do
% this, and sometime the code below doesn't work....
% it relies on knowledge of the different ellipse colors
hold on
[C,h] = contourf(xv,yv,Z,[0,1],'LineColor','none','LineWidth',2);
% cmap = [[234,228,241]./255; 1 1 1];
cmap = [[237,221,221]./255; 1 1 1];
colormap(cmap)

fastGuaranteedPlot = findobj(gca, 'Color', [1 0 0]);
uistack(fastGuaranteedPlot, 'top')
%truthAndDataPlot = findobj(gca, 'Color', [0 1 0]);
directEllipsePlot = findobj(gca, 'Color', [0 0 0]);
uistack(directEllipsePlot, 'top') 
dataPointPlot = findobj(gca, 'Color', [0 0 0]);
uistack(dataPointPlot, 'top')                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%% Example with data points from elongated ellipse and  %
%                       user-specified covariance matrices                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('**************************************************************\n')
fprintf('* Example with ALL data points assuming inhomogeneous        *\n')
fprintf('* Gaussian noise, where a covariance matrix is provided      *\n')
fprintf('* for each data points                                       *\n')
fprintf('**************************************************************\n')        
data_points_eccentric = [ 362.5,  162.5;
          327.993236135912,          203.769681636721;
          285.034694946748,          250.049135605403;
          238.279600940802,         296.323264968822;
          192.794594033649,          337.577549793355;
          153.508680496279,          369.341448192079;
          124.679098644223,          388.172848905829;
          109.429981008205,           392.03107746798;
          109.413806098526,         380.498034947239;
           124.63232671763,          354.823505454321;
          153.436380016758,          317.789722641156;
          192.704599880014,          273.409871513897;
          238.181665376003,          226.493197566525;
          284.939430808858,          182.123850408441;
          327.910966772446,          145.109937265579;
          362.439640573031,          119.462489957005;
          384.783734636744,          107.960807359753;
           392.52191966744,          111.851275250228;
          384.815643246731,          130.712301065668;
                     362.5,                     162.5];

% We add inhomogeneous isotropic noise to all of the data points
% and store the covariance matrices for each data point
for iPts = 1:length(data_points_eccentric)
    standardDevX = 1 + (6-1)*rand(1,1);
    standardDevY = 1 + (6-1)*rand(1,1);
    covList{iPts} = diag([standardDevX^2, standardDevY^2]);
    perturbationX = standardDevX * randn(1,1);
    perturbationY = standardDevX * randn(1,1);
    data_points_eccentric(iPts,:) = ...
        data_points_eccentric(iPts,:) + [perturbationX, perturbationY];    
end

% An example of algebraic ellipse fitting  to data points sampled
% from an eccentric ellipse 
fprintf('Algebraic ellipse parameters of direct ellipse fit: \n')
[theta_dir]  = compute_directellipse_estimates(data_points_eccentric)

% An example of Sampson distance ellipse fitting to data points
% sampled from an eccentric ellipse assumping isotropic homogeneous noise
fprintf('Algebraic ellipse parameters of our method assuming \n')
fprintf('homogeneous Gaussian noise with noise level automatically \n')
fprintf('estimated from the data: \n')
 [theta_fastguaranteed] = ...
           fast_guaranteed_ellipse_estimate(data_points_eccentric);
                        

% An example of Sampson distance ellipse fitting to data points 
% sampled from an eccentric ellipse assuming isotropic inhomogeneous noise
% with user-specified covariance matrices
fprintf('Algebraic ellipse parameters of our method using \n')
fprintf('specified data point covariance matrices: \n')
[theta_guaranteed_covweighted] = ...
          fast_guaranteed_ellipse_estimate(data_points_eccentric,...
                                                                  covList);       

fprintf('Geometric ellipse parameters \n')
fprintf('(majAxis, minAxis, xCenter,yCenter, orientation (radians)): \n') 
geometricEllipseParameters = ...
           fromAlgebraicToGeometricParameters(theta_guaranteed_covweighted)            

fprintf('Covariance matrix of geometric parameters: \n')
geoCov = compute_covariance_of_geometric_parameters(...
                         theta_guaranteed_covweighted, data_points_portion)

fprintf('Standard deviation of geometric parameters: \n')                          
stds = sqrt(diag(geoCov))                                                                  
                                                              
[S, thetaCovarianceMatrixNormalisedSpace] = ...
                                compute_covariance_of_sampson_estimate(...
              theta_guaranteed_covweighted, data_points_eccentric,covList);
                
                
% plot the data points
x = data_points_eccentric;
xportion = data_points_eccentric;
n = length(x);
figure('Color',[1 1 1])
%plot(x(:,1),x(:,2),'b.');
plot(xportion(:,1),xportion(:,2),'k.');
hold on

% determine data range
% minX = min(min(x(:,1))) - 20;
% minY = min(min(x(:,2))) - 20;
% maxX = max(max(x(:,1))) + 20;
% maxY = max(max(x(:,2)))  + 20;

% the particular way we plot confidence regions in matlab
% requires that we work with square axes
minX = 50;
maxX = 450;
minY = 50;
maxY = 450;

% plot the direct ellipse fit
a = theta_dir(1); b = theta_dir(2); c = theta_dir(3);
d = theta_dir(4); e = theta_dir(5); f = theta_dir(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
set(h, 'Color', [0 0 0]);
set(h,'LineWidth',1.5); 
axis([minX maxX minY maxY]);

% plot the guaranteed ellipse fit
a = theta_fastguaranteed(1); b = theta_fastguaranteed(2);
c = theta_fastguaranteed(3); d = theta_fastguaranteed(4); 
e = theta_fastguaranteed(5); f = theta_fastguaranteed(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
set(h, 'Color', [1 0 0]);
set(h,'LineWidth',1.5); 
axis([minX maxX minY maxY]);



% plot the guaranteed ellipse fit using covariance weights (inhomogeneous)
a = theta_guaranteed_covweighted(1); b = theta_guaranteed_covweighted(2); 
c = theta_guaranteed_covweighted(3);d = theta_guaranteed_covweighted(4);
e = theta_guaranteed_covweighted(5); f = theta_guaranteed_covweighted(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);
set(h, 'Color', [0 0 1]);

title('All data and inhomogeneous Gaussian noise')
legend('DATA POINT','DIRECT ELLIPSE FIT','FAST GUARANTEED ELLIPSE FIT',...
      'FAST GUARANTEED ELLIPSE FIT WITH COVARIANCE WEIGHTS',...
            'Location','SouthOutside');
        
        
                                                    
xres=400;
yres=400;
xv = linspace(minX, maxX, xres);
yv = linspace(minY, maxY, yres);

% scale and translate data points so that they lie inside a unit box
% in this instance we are actually just after the transformation matrix T
[~, T] = normalize_data_isotropically(data_points_eccentric);

% convert parameters to normalised space

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
        
theta_guaranteed_covweightedNormalised = E \ P34*pinv(D3)*...
            inv(kron(T,T))'*D3*P34*E*theta_guaranteed_covweighted;


% now we can transform our grid points to normalised space
pts = [xv; yv; ones(1,xres)];
pts_normalised_space = T*pts;
xv_normalised = pts_normalised_space(1,:);
yv_normalised = pts_normalised_space(2,:);
Z = compute_confidence_band(xv_normalised,yv_normalised,...
                         theta_guaranteed_covweightedNormalised, ...
                         thetaCovarianceMatrixNormalisedSpace,11.07);

% this next step is a matlab hack to try to overlay the confidence
% region over the existing plot. I don't know of a better way to do
% this, and sometimes the code below doesn't work....
% it relies on knowledge of the different ellipse colors
hold on
[C,h] = contourf(xv,yv,Z,[0,1],'LineColor','none','LineWidth',2);
% cmap = [[234,228,241]./255; 1 1 1];
cmap = [[221,221,237]./255; 1 1 1];
colormap(cmap)

fastGuaranteedCovWeightedPlot = findobj(gca, 'Color', [0 0 1]);
uistack(fastGuaranteedCovWeightedPlot, 'top')
fastGuaranteedPlot = findobj(gca, 'Color', [1 0 0]);
uistack(fastGuaranteedPlot, 'top')
directEllipsePlot = findobj(gca, 'Color', [0 0 0]);
uistack(directEllipsePlot, 'top') 
dataPointPlot = findobj(gca, 'Color', [0 0 0]);

uistack(dataPointPlot, 'top')                                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

