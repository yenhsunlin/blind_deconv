Thie repository includes the Bayesian blind-deconvolution scripts. One is in principle capable of using such script to re-focus the picture deteoriated by the motion blur.

Note that the motion blur kernel should be spatial invariant and contains no rotation. The scripts in this repository are experimental and case sensitive. Different pictures and blur kernels may require specific fine-tuning on the hyper parameters.

The de-blurred result is also not promising and while imposing hyper-Laplacian prior (0<p<1) as the regularizer (penalty term), the scipy optimizer module could fail to find the correct result. 
