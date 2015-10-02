# Guaranteed Ellipse Fitting with a Confidence Region and an Uncertainty Measure for Centre, Axes and Orientation
MATLAB implementation of Guaranteed ellipse fitting with a confidence region and an uncertainty measure for centre, axes, and orientation described in the publication:

Z. L. Szpak, W. Chojnacki, and A. van den Hengel. 
Guaranteed ellipse fitting with a confidence region and an uncertainty measure for centre, axes, and orientation. 
J. Math. Imaging Vision, 52(2):173-199, 2015. 

In this work we utilise the Sampson Distance to measure the quality of fit between a candidate ellipse and data points, hence the method is still very fast and simple to implement. Our method also produces a measure of uncertainty in the form of a covariance matrix for the estimated parameters. In addition, we provide formulae for converting the algebraic ellipse parameters into geometric ellipse parameters such as the semi-major and semi-minor axes, the centre of the ellipse and the orientation of the ellipse (in radians). Moreover, we propagate the measure of uncertainty through to the geometric parameters and produce a covariance matrix for the geometric parameters. Finally, we produce a planar 95% confidence region for the ellipse estimate, to assist in interpreting the uncertainty contained in the covariance matrices. 

If you have covariance matrices associated with your data points you can specify them and our method will utilise them in the estimation process. If you do not submit data point covariance matrices, our method assumes homogeneous Gaussian noise and estimates the level of noise automatically from the data.
