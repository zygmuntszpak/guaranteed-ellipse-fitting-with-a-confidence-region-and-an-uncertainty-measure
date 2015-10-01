% function   Z = compute_confidence_band(xv,yv,...
%                       tNormalisedSpace,covNormalisedSpace,criticalValue )
%
%
% Author: Zygmunt L. Szpak (zygmunt.szpak@gmail.com)
% Date: September 2014
%
% Description:  This procedure samples point over a regular grid and tests 
%               each point to see if it could plausibly have originated
%               from the ellipse parameterised by 'tNormalisedSpace'
%               with covariance matrix 'covNormalisedSpace'.
%               The test cut-off value is specified by 'criticalValue'. 
%               Points that pass the test form the confidence region for
%               the ellipse.
%               
%
%              
%
% Parameters :  xv                        -  vector specifing x-coordinates
%                                            that will be used to form
%                                            a grid of 2D coordinates.
%
%               yv                        -  vector specifing y-coordinates
%                                            that will be used to form a
%                                            grid of 2D coordinates.
%
%               criticalValue             -  a threshold on the Chi-Squared
%                                            value associated with
%                                            a chosen level of confidence
%                                            and p-value with 5 degrees
%                                            of freedom. Refer to
%                                            http://en.wikipedia.org/wiki/
%                                            Chi-squared_distribution
%                                            for more details.
%                                            
%
%
%
%
% Return     : a grid (represented as a matrix) with a value of 0
%              if a point could have originated from the ellipse
%              and a value of 1 otherwise. 
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
function Z = compute_confidence_band(xv,yv,...
    tNormalisedSpace,covNormalisedSpace,criticalValue )


[Xinterp,Yinterp] = meshgrid(xv,yv);
Z = ones(length(xv),length(yv));

for i = 1:length(xv)
    for j = 1:length(yv)
       
        x = Xinterp(i,j);
        y = Yinterp(i,j);        

        tNormalisedSpace = tNormalisedSpace /  norm(tNormalisedSpace);
       
        ux = [x^2 x*y y^2 x y 1]';
        zSquared = tNormalisedSpace'*(ux*ux')*tNormalisedSpace;
        var = ux'*covNormalisedSpace*ux;
        %val = zSquared - pval*var;        
        
        if (zSquared/var < criticalValue)
            Z(i,j) = 0;
        end
        
    end
end

end

