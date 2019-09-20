# Bayesian blind-deconvolution (*experimental*)
Thie repository includes the *experimental* Bayesian blind-deconvolution scripts. One is in principle capable of using such scripts to re-focus the picture deteoriated by the motion blur without knowing the blur kerner in advance.

Note that the motion blur kernel should be spatial invariant and contains no rotation. The scripts in this repository are experimental and case sensitive. Different pictures and blur kernels may require fine-tuning on the hyper parameters *ad hoc*.

## Required packages
To run these scripts, the followings third-party packages are required:
- `cv2`
- `numpy`
- `pypher`
- `scipy`

Both `cv2` and `pypher` can be installed via `pip`.

## Known issues
The de-blurred result is not promising in general and while imposing hyper-Laplacian norm *L<sub>p</sub>* (0<*p*<1) as the regularizer (penalty term), the `fsolve` in `scipy.optimize` module could fail to find the correct result. 

## Important notes!
This blind-deconvolution scripts is just for fun and part of the core contents are transcribed from J. Kotera et al.'s work, see `DOI: 10.1007/978-3-642-40246-3_8` for further details. If you have other purpose for commercial using or deploying into your own projects instead of fun-playing or academic research, pleas DO contact the original authors for permission! 

If you found any glitch while running these scripts, please contact via <a href='mailto:yenhsun@gate.sinica.edu.tw'>`yenhsun@gate.sinica.edu.tw`</a> or <a href='mailto:yenhsun@phys.ncku.edu.tw'>`yenhsun@phys.ncku.edu.tw`</a>, much thanks!

More about me @ <a href='https://sites.google.com/view/yenhsun' title='Google site'>*Google site*</a>.
