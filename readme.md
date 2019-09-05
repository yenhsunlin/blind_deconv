<script type="text/javascript" async
src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js? 
config=TeX-MML-AM_CHTML"
</script>

# Bayesian image process: Blind-deconvolution
Thie repository includes the Bayesian blind-deconvolution scripts. One is in principle capable of using such script to re-focus the picture deteoriated by the motion blur without knowing the blur kerner in advance.

Note that the motion blur kernel should be spatial invariant and contains no rotation. The scripts in this repository are experimental and case sensitive. Different pictures and blur kernels may require specific fine-tuning on the hyper parameters.

The de-blurred result is also not promising and while imposing hyper-Laplacian norm $ L_p $ ($ 0<p<1 $) as the regularizer (penalty term), the `fsolve` in `scipy.optimize` module could fail to find the correct result. 

If you found any glitch while running these scripts, please contact via <a href='mailto:yenhsun@gate.sinica.edu.tw'>`yenhsun@gate.sinica.edu.tw`</a> or <a href='mailto:yenhsun@phys.ncku.edu.tw'>`yenhsun@phys.ncku.edu.tw`</a>, much thanks!

More about me @ <a href='https://sites.google.com/view/yenhsun' title='Google site'>*Google site*</a>.
