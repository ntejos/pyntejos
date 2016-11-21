import numpy as np

"""This module is meant to deal with the Soneira-Peebles fractal
distribution, with applications to astrophysics

Terminology
-----------

gamma : slope of the correlation
L : Number of fractal levels 
R : radius of the 0-th level sphere (this define global scale)
eta : Number of spheres in a given level
l : shrinking factor of sphere radii from one level to the other R_i = R_(i-1)/l
M : number of dimensions

r0 : correlation length (will depend on the scale R)

The approach will be to chose the number of levels.
"""


def get_eta(l, gamma, M):
    log_l = np.log10(l)
    return 10 ** (log_l * (M - gamma))


def get_lambda(eta, gamma, M):
    log_eta = np.log10(eta)
    return 10 ** (log_eta / (M - gamma))

class Sphere(object):
    def __init__(self, pos, r, M=1):
        self.pos = pos
        self.radius = r
        self.dimension = M
    
    def get_random_pos(N): 
        """Return N random position within the sphere"""
        if self.dimension == 1:
            return np.random.uniform(pos-r,pos+r,N)

# lets try a simple case first
eta = 2 # Number of spheres per level
L = 3 # Number of levels
l = 1.5 # shrinking factor
R = 1 #initial radius

# total number of particles, including intermediate levels
Npart = np.sum(eta ** np.arange(0, L+1))
pos = [(0,0,0)] * Npart

Ngal = eta ** L  # Number of galaxies in the final level
gal = [(0,0,0)] * Ngal # positions of Ngal galaxies

ng = 1
ic = 0
ig = 0

#set the original random seed
np.random.seed(30)

for i in range(L+1):
    print i
    Ri = R / (l ** i)
    Nc_prev = eta ** (i - 1)
    Nc = eta ** (i + 1)

    if i == 0:
        pos[i][0] = np.random.uniform(0,1,1) # pos in the x-axis
        # pos[i][1] = np.random.uniform(0,1,1) # pos in the y-axis
        # pos[i][3] = np.random.uniform(0,1,1) # pos in the z-axis

    for j in range(0, Nc_prev):
        xc = pos[ic][0]       
        # yc = pos[ic][1]       
        # zc = pos[ic][2]

        for k in range(1, eta+1):
            xg = np.random.uniform(xc - Ri, xc + Ri, 1)
            # ylen = np.sqrt(Ri**2 -  )   
            



def soneira_peebles_dist(L, gamma, d=1):
    """Gives the positions of N points with two-point correlation
    function following a slope of -gamma, following the model of Soneira
    and Peebles, currently only implemented for dimension (d=1)

    Parameters
    ----------
    N : int
        Number of points to retrieve
    gamma : float
        Slope of the correlation function ~r**(-gamma)
    d : int, optional
        Number of dimensions.

    Returns
    -------
    N_corr : array
        The relative positions in a dominion (0, 1) of N points
        following the given correlation power

    """

    eta = 2 # Number of spheres per level
    N = int(eta)**int(L) # Final number of points
    l = 1.5 # shrinking factor

    N_pos = np.sum(eta ** np.arange(0,))

    int(eta)**int(L) + 1 # Number of levels including intermediate



"""
; This code computes the Soneira-Peebles algorithm
; to generate a spatial distribution of objects with a 
; given set of parameters

; Last Update: 13-Oct-2008

pro soneira_peebles,Param,Gal,NGal

OutFile = '../out/test'


; Debug Parameters
    Plot  =  0      ; It will plot the result for each level
    Debug =  0
    Stop  =  0
; General parameters

    


;   R       = 25.   ; Size of the zero-level sphere
;   eta     = 5.    ; Number of sub-spheres 
;   L       = 3.    ; Number of levels
;   Lambda  = 3     ; Radius = R/Lambda^{L-1}

;   Seed = 31235

    R       = Param[0]
    eta     = Param[1]
    L       = Param[2]
    Lambda  = Param[3]  
    Seed    = Param[4]  
    XSize   = Param[5]
    YSize   = Param[6]



; 2D Example

;   XSize = 50.
;   YSize = 50.

    NPart = total(eta^(indgen(L+1)+1) + 1

    Pos = replicate({x:0.,y:0.},NPart)
    
    NGal = eta^(L+1)
    Gal  = replicate({x:0.,y:0.},NGal)

    if Debug then begin
        print,'NPart = ',NPart
        print,'NGal = ',NGal
    endif
    

    if Plot then begin
;       window,0,retain=2
        loadct,39
        plotsym,0,/fill
        plot,indgen(2),/nodata,xr=[0,XSize],[0,YSize],/isotropic,/xs,/ys
    endif   
        
    ng = 1l
    ic = 0l
    ig = 0l
    for i = 0,L do begin

        R0 = R/(Lambda^i)           ; Radius of each sphere at the i-level
        Nc_prev = eta^(i-1)             ; Number of spheres at the i-1 level
        Nc = eta^(i+1)              ; Number of spheres at the i-level
    
        if i eq 0 then begin
            Pos[0].x = randomu(seed)*XSize
            Pos[0].y = randomu(seed)*YSize
            if Debug then print,'Zero-sphere coordinates: ',Pos[0].x,Pos[0].y
            continue
        endif
        
        for j = 0l,Nc_prev-1 do begin   

            xc = Pos[ic].x
            yc = Pos[ic].y

            for k = 1,eta do begin
                xg = randomu(Seed)*(2*R0) + (xc-R0)
                ylen = sqrt(R0^2 - (xg - xc)^2)
                yg = randomu(Seed)*(2*ylen) + (yc -ylen)

                Pos[ng].x = xg
                Pos[ng].y = yg

                if i eq L then begin    ;Here the galaxy positions are recorde
d
                    Gal[ig].x = xg
                    Gal[ig].y = yg
                    ig++
                endif

                if Plot then begin
;                   oplot,[xg,xg],[yg,yg],ps=3
                    polyfill,circle(xg,yg,r0),/fill,color=i*40+30
                    plots,circle(xg,yg,r0)
;                   xyouts,xg,yg,strn(ng),charsiz=.7
                endif
                ng++
                if Debug then begin
                    print,'ng = ',ng
                    print,'ic = ',ic
                    print,'i, j, k = ',i,j,k
                    print,'R0 = ',r0
    ;               stop,''

                endif

            endfor
if Stop then stop,''    
        ic++
        endfor
    R0_prev = R0
    endfor

    return

end     
"""