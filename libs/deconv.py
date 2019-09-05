import numpy as np
from pypher.pypher import psf2otf
from copy import copy, deepcopy
from scipy.optimize import fsolve
from scipy.signal import correlate2d, convolve2d


# abbreviations
fft2 = np.fft.fft2
ifft2 = np.fft.ifft2
conj = np.conj
real = np.real
conv2 = convolve2d


def LnormPrior(x, nx, Lp, alpha, beta):
    
    def aLn(DU, normDU, u_star, k):
        V = np.zeros_like(DU)
        DUp = DU[normDU > u_star]
        normDUp = normDU[normDU > u_star]
        V[normDU > u_star] = DUp*(normDUp - k)/normDUp
        return V
    
    if Lp == 1:
        v_star = 0
        u_star = alpha/beta
    elif Lp == 0:
        v_star = np.sqrt(2*alpha/beta)
        u_star = np.sqrt(2*alpha/beta)
    elif 0 < Lp <1:
        v1 = lambda v : -v+alpha/beta*v**(Lp-1)*(1-Lp)*Lp
        v2 = lambda v : -0.5*v**2+alpha/beta*v**Lp*(1-Lp)
        lbound = fsolve(v1,10)[0]
        v_star = fsolve(v2, [lbound, 10])[1]
        u_star = v_star + alpha/beta*Lp*v_star**(Lp-1)    
    return aLn(x, nx, u_star, u_star - v_star)


def PSFest(I, psf,  max_it, max_it_u, max_it_h, gamma, Lp, ccreltol):
    Lp = Lp
    ccreltol = ccreltol
    gamma = gamma
    
    U = np.zeros_like(I)
    H = psf
    # Initializing params for min-U steps
    # FFT of U, the latent image
    FU = fft2(U)
    # FFT of x and y directives
    FDx = fft2(np.array([[1,-1]]), I.shape) 
    FDy = fft2(np.array([[1],[-1]]), I.shape) 
    DTD = conj(FDx)*FDx + conj(FDy)*FDy
    # FFT of x and y directives of U
    FUx = np.zeros_like(I)
    FUy = np.zeros_like(I)
    
    # Auxiliary vars for image gradients and blurs
    # initialize to zeros
    Vx = np.zeros_like(I)
    Vy = np.zeros_like(I)
    Vh = np.zeros_like(I)
    # Bergman iterations
    Bx = np.zeros_like(I)
    By = np.zeros_like(I)
    Bh = np.zeros_like(I)
    
    FeGu = fft2(I)
    FeGx = FDx*FeGu # FFT of dU/dx
    FeGy = FDy*FeGu # FFT of dU/dy
    
    for i in range(max_it):
        
        # U-estimation
        FHS = fft2(H, I.shape)#psf2otf(psf, I.shape)
        FHTH = conj(FHS)*FHS
        FGs = conj(FHS)*FeGu
        
        beta_u = 3*1e-2*gamma      
        alpha_u = 3*1e-4*gamma  
        
        for i in range(max_it_u):
            FUp = FU.copy()
            b = FGs + beta_u/gamma*(conj(FDx)*fft2(Vx + Bx) + conj(FDy)*fft2(Vy + By))
            FU = b/(FHTH + beta_u/gamma*DTD)
            
            FUx = FDx*FU
            FUy = FDy*FU
            xD = real(ifft2(FUx))
            yD = real(ifft2(FUy))
            xDm = xD - Bx
            yDm = yD - By
            nDm = np.sqrt(xDm**2 + yDm**2)
            Vy = LnormPrior(yDm, nDm, Lp, alpha_u, beta_u)
            Vx = LnormPrior(xDm, nDm, Lp, alpha_u, beta_u)
            # update Bregman vars
            Bx = Bx - xD  + Vx
            By = By - yD  + Vy
            #E = np.sqrt(Vy**2 + Vx**2)
            
            relcon = np.sqrt(np.sum(np.abs(FUp - FU)**2)/np.sum(np.abs(FU)**2))
            if relcon < ccreltol:
                break
            else: pass
            # update latent image
        U = real(ifft2(FU))
        
        # h-estimation
        FUD = FeGx*conj(FUx) + FeGy*conj(FUy)
        FUTU = conj(FUx)*FUx + conj(FUy)*FUy
        FH = fft2(H, I.shape) #fft2(psf,I.shape) #psf2otf(psf, I.shape)
        # forcing windows outside kernel size are zero
        mask = np.zeros_like(I)
        mask[:H.shape[0],:H.shape[1]] = 1
        
        beta_h = 1e4*gamma
        alpha_h = 5*1e0*gamma
        
        for i in range(max_it_h):
            FHp = FH.copy()
            b = beta_h/gamma*fft2(Vh + Bh) + FUD
            FH = b/(FUTU + beta_h/gamma)
            relcon = np.sqrt(np.sum(np.abs(FHp - FH)**2)/np.sum(np.abs(FH)**2))
            
            hI = real(ifft2(FH)) #real(ifft2(FH)) #otf2psf(FH, I.shape)
            hIm = hI - Bh
            nIm = np.abs(hIm)
            Vh = LnormPrior(hIm, nIm, 1, alpha_h, beta_h)
            Vh[Vh < 0] = 0
            Vh = Vh*mask
            Bh = Bh - hI  + Vh
            H = hI[:H.shape[0],:H.shape[1]]
            if relcon < ccreltol:
                break
            else: pass
        
        # Increasing gamma at every step
        gamma = gamma*1.8
    
    return U, H


def deconv(I, psf, max_it_u, gamma, ccreltoltol):
    gamma = gamma
    beta_u = 3*1e-2*gamma      
    alpha_u = 3*1e-4*gamma
    ccreltol = ccreltoltol
    #Lp = Lp
    
    FU = 0
    FDx = fft2(np.array([[1,-1]]), I.shape) 
    FDy = fft2(np.array([[1],[-1]]), I.shape) 
    FH = psf2otf(psf, I.shape) #fft2(psf, I.shape)
    FHTH = conj(FH)*FH
    
    FGu = fft2(I)
    FGs = conj(FH)*FGu
    
    DTD = conj(FDx)*FDx + conj(FDy)*FDy
    
    # Bregman vars
    Bx = np.zeros_like(I)
    By = np.zeros_like(I)
    Vx = np.zeros_like(I)
    Vy = np.zeros_like(I)
    
    for i in range(max_it_u):
        FUp = FU
        b = FGs + beta_u/gamma*(conj(FDx)*fft2(Vx + Bx) + conj(FDy)*fft2(Vy + By))
        FU = b/(FHTH + beta_u/gamma*DTD)
        
        xD = real(ifft2(FDx*FU))
        yD = real(ifft2(FDy*FU))
        xDm = xD - Bx
        yDm = yD - By
        nDm = np.sqrt(xDm**2 + yDm**2)
        Vy =  LnormPrior(yDm, nDm, 1, alpha_u, beta_u)
        Vx =  LnormPrior(xDm, nDm, 1, alpha_u, beta_u)
        
        Bx = Bx + Vx - xD
        By = By + Vy - yD
        
        relcon = np.sqrt(np.sum(np.abs(FUp - FU)**2)/np.sum(np.abs(FU)**2))
        if relcon < ccreltol:
            break
        else: pass
            
    U = real(ifft2(FU))
    U = normalize_im(U)[0]
    return U
        

# Auxiliary functions
def normalize_im(im):
        """
        Constraint the intensity values between 0 and 1
        
        Input
        -----
        im: image-array
        
        Return
        ------
        Image array, Original min intensity, Original max intensity
        """
        newim = np.float32(im)
        #newIm = np.zeros_like(newim)
        maxIm = np.max(newim)
        minIm = np.min(newim)
        newIm = (newim - minIm)/(maxIm - minIm)
        return newIm, minIm, maxIm-minIm


def chop_im(im, chop_size):
    """
    Return image with both sides <= 1024
    """
    row, col = im.shape
    if row > chop_size and col > chop_size:
        G = im[row//2 - chop_size//2, row//2 + 1 + chop_size//2, col//2 - chop_size//2, col//2 + 1 + chop_size//2]
    elif row > chop_size:
        G = im[row//2 - chop_size//2, row//2 + 1 + chop_size//2,:]
    elif col > chop_size:
        G = im[:, col//2 - chop_size//2, col//2 + 1 + chop_size//2]
    return G