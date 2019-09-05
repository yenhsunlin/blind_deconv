import sys
sys.path.append('./')
import numpy as np
from copy import copy,deepcopy
import cv2
from deconv import *

def bayes_deblur(im, psf_size, maxiter=5, maxit_u = 5, maxit_h = 5, Lp = 0, psf=None, params=None):
    """
    Blind-deconvolution algorithm based on alternating MAP estimation with
    heavy-tailed priors.
    This python version transcribes from J. Kotera et al.'s work, see
    DOI: 10.1007/978-3-642-40246-3_8 for further details.
    At this moment, this transcription supports Lp = 0 and 1 only.
    
    Input
    -----
    im: array, blurred image color or grey
    psf_size: int, a guess on the psf_size, must be an odd number to have an
              anchor point
    maxiter: int, number for iterations
    psf: array, ground truth PSF function for later comparison, default is None
    params: dict, changing the default parameters
    
    Return
    ------
    Deblur image, estimated PSF
    """
    
    # Define parameters
    gamma = 3e2
    #Lp = 0
    tol = 1e-3
    iter_u = maxit_u
    iter_h = maxit_h
    chop_size = 1024
    
    # Check the image is color or grey
    if im.ndim == 3 and im.shape[-1] == 3:
        color = True
        # convert color image into grey-scale for blind-deconvolution
        G = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    elif im.ndim == 2:
        color = False
        # copy a new image
        G = im.copy()
    else:
        raise ValueError('The inpupt image must have dimension 3 for color or dimension 2 for grey.')
    
    # Normalize the image intensity between 0 and 1 and get the original min and max
    # intensity (should between 0 and 255)
    G, norm_m, norm_v = normalize_im(G)
    # If the blurred image is colorful, we have to prepare additional stuff for later use
    if color:
        norm_m = []
        norm_v = []
        colorG = np.zeros_like(im, dtype=np.float32)
        for i in range(3):
            colorG[:,:,i], m, v = normalize_im(im[:,:,i])
            norm_m.append(m)
            norm_v.append(v)
    else: pass
    
    # Chop image if any side of its W, H is greater than 1024
    if G.shape > (chop_size, chop_size):
        G = chop_im(im, chop_size)
    else:
        pass
    
    # Initializing PSF function as a delta function
    psf_guest = np.zeros((psf_size,psf_size))
    psf_guest[psf_size//2,psf_size//2] = 1
    
    # Starting PSF estimation
    latentG, psf_est = PSFest(G, psf_guest, maxiter, iter_u, iter_h, gamma, Lp, tol)
    
    # Regularize estimated PSF for non-blind decovolution process
    psf_est[psf_est < 0] = 0
    psf_est = psf_est/np.sum(psf_est)
    
    # Starting deconvolution process through estimated PSF
    if color:
        deconvG = np.zeros_like(colorG)
        for ch in range(3):
            deconvG[:,:,ch] = deconv(colorG[:,:,ch], psf_est, maxiter, gamma, tol)
            deconvG[:,:,ch] = deconvG[:,:,ch]*norm_v[ch] + norm_m[ch]       
    else:
        deconvG = deconv(G, psf_est, maxiter, gamma, tol)
        deconvG = deconvG*norm_v + norm_m
    
    return np.uint8(deconvG), psf_est
