===================================================================================
Timothy D Barfoot and Paul T Furgale, 
Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
DOI: 10.1109/TRO.2014.2298059
tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
===================================================================================

These files provide an implementation of the major equations from the
paper mentioned above. For formulas with both series and closed form
representations, both implementations are provided.


Experiment Scripts
------------------		
compound_exp1.m		generates plots for first pose compound experiment in paper (figure 2)
compound_exp2.m		generates plots for second pose compound experiment in paper (figure 3)
fusion_exp1.m		generates plots for first pose fusion experiment in paper (figure 5)
fusion_exp2.m		generates plots for second pose fusion experiment in paper (figure 6)
measure_exp1.m		generates plots for passing uncertain pose through stereo camera model in paper (figures 7 and 8)
measure_exp2.m		(not in paper) generates 3D plot showing effect of uncertain pose applied to constant vector


Math Implementations
--------------------
bernoullinumber.m  	provides Bernoulli number sequence
cam.m  				stereo camera model
camHess.m  			Hessian of stereo camera model
camInv.m  			inverse of stereo camera model
camJac.m			Jacobian of stereo camera model
compound.m			compound two uncertain pose transformations into one
curlyhat.m			turns a 6x1 column into associated Adjoint
hat.m				turns a 3x1 into skew-symmetric matrix or 6x1 into special 4x4 matrix
plotcov.m			plots a covariance matrix
point2fs.m			turns a 4x1 homogeneous point into a special 4x6 matrix
point2sf.m			turns a 4x1 homogeneous point into a special 6x4 matrix
rot2vec.m			turns a 3x3 rotation matrix into a 3x1 rotation vector
rotValidate.m		checks if a 3x3 rotation matrix is valid or not
tran2vec.m			turns a 4x4 transformation matrix into a 6x1 column
tranAd.m			turns a 4x4 transformation matrix into a 6x6 Adjoint matrix
tranValidate.m		checks if a 4x4 transformation matrix is valid or not
vec2Q.m				computes the 3x3 Q matrix from a 6x1 column
vec2jac.m			computes the 3x3 and 6x6 Jacobian matrices for SO(3) and SE(3) in closed form
vec2jacInv.m		computes the inverse of the 3x3 and 6x6 Jacobian matrices for SO(3) and SE(3) in closed form
vec2jacInvSeries.m	approximates the inverse of the 3x3 and 6x6 Jacobian matrices for SO(3) and SE(3) in series form
vec2jacSeries.m		approximates the 3x3 and 6x6 Jacobian matrices for SO(3) and SE(3) in series form
vec2rot.m			computes a 3x3 rotation matrix from a 3x1 rotation vector in closed form
vec2rotSeries.m		approximates the 3x3 rotation matrix from a 3x1 rotation vector in series form
vec2tran.m			computes a 4x4 transformation matrix from a 6x1 column in closed form
vec2tranSeries.m	approximates a 4x4 transformation matrix from a 6x1 column in series form
