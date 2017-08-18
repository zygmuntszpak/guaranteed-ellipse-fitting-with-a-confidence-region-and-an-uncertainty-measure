Guaranteed Ellipse Fitting with a Confidence Region and an Uncertainty
Measure for Center, Axes, and Orientation 
Zygmunt L Szpak (c) 2014 

Version 2, 18 August 2017 
Version 1, 25 September 2014 

If you find this code useful in your research, please cite:

Z. L. Szpak, W. Chojnacki, and A. van den Hengel. 
Guaranteed ellipse fitting with a confidence region and an uncertainty
measure for centre, axes, and orientation. 
J. Math. Imaging Vision, 2015. 
http://dx.doi.org/10.1007/s10851-014-0536-x

CHANGES 

Version 2 Bug fix

The covariance matrix of the algebraic parameters associated 
with the original (unnormalised) coordinate system was incorrect due to 
an oversight in the implementation. I incorrectly unit normalised a 
parameter vector when I shouldn't have just before applying a critical 
formula. 

This issue did not affect the covariance matrix of the geometric 
parameters, nor did it affect the covariance matrix of the algebraic 
parameters in the normalised coordinate system. Hence, if you used the 
normalised algebraic covariance matrix or the covariance matrix of the 
geometric parameters for any subsequent calculations/propagations your 
results are correct. However, if for some reason you utilised the 
covariance of the algebraic parameters in the original (unnormalised) 
coordinate system then you would have had an inflated covariance matrix. 
That is, your error bars would have been greater than they needed to be. 

The issue has been resolved in the current release. 



Version 1 Initial submission 

DISCLAIMER OF WARRANTY 

This source code is provided "as is" and without warranties as to 
performance or merchantability. The author and/or distributors of this 
source code may have made statements about this source code. Any such 
statements do not constitute warranties and shall not be relied on by 
the user in deciding whether to use this source code. 

This source code is provided without any express or implied warranties 
whatsoever. Because of the diversity of conditions and hardware under 
which this source code may be used, no warranty of fitness for a 
particular purpose is offered. The user is advised to test the source 
code thoroughly before relying on it. The user must assume the entire 
risk of using the source code. 

HOW TO USE 

Unzip the folder into a directory of your choice and add the folder and 
all its sub-folders to the matlab path. Then execute run_example. 

This will run the guaranteed ellipse fit and the direct ellipse fit on 
synthetic data and produce a graph comparing the two methods. 



KNOWN ISSUES 

* None 

FUTURE RELEASE 

TBD


