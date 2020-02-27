'''

Create a r-theta grid, using Chebyshev points for r and Gauss-Legendre points for theta.
Define the Chebyshev and Gauss-Legendre weights for integration using quadrature rules 
Define the Chebyshev differentiation matrix for taking derivatives on the r-grid.
Define the Legendre polynomials and their derivatives for spherical harmonic expansions
and taking derivatives on the theta-grid.

'''

# Author: Laura Kulowski
# Last updated: Nov 2019 

# SciPy imports 
import numpy as np

class Grid():
    
    def __init__(self, nr, nt, r_inner, r_outer):
        
        '''
        : param nr: number of radial grid points (chebyshev)
        : param nt: number of theta grid points (gauss-legendre)
        : param r_inner: inner boundary of the radial grid (0 < r_inner < 1)
        : param r_outer: outer boundary of the radial grid (0 < r_outer <= 1)
        '''
        
        self.nr = nr
        self.nt = nt
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.leg_deg = self.nt
        self.xc = None
        self.xcs = None
        self.x = None
        self.w = None
        self.wr = None
        self.theta = None
        
    def cgl_points(self):
        
        '''
        chebyshev-gauss-legendre (cgl) points in interval [-1, 1]
        : return self.xc: chebyshev-gauss-legendre points
        
        '''
                
        js = np.linspace(0, self.nr, self.nr + 1)
        self.xc = np.cos(js * np.pi / self.nr)  
        self.xc = self.xc[::-1]
        
        return self.xc

    def cheby_shift(self):
        
        '''
        chebyshev-gauss-legendre points shifted into interval [a, b]
        : return self.xcs: shifted chebyshev-gauss-legendre points
        '''
        
        js = np.linspace(0, self.nr, self.nr + 1)
        self.xcs = 0.5 * (self.r_inner + self.r_outer) + 0.5 * (self.r_outer - self.r_inner) * np.cos(js * np.pi / self.nr)
        self.xcs = self.xcs[::-1]
        
        return self.xcs

    def gauss_leg_xw(self):
        
        '''
        gauss-legendre points and weights
        : return [self.x, self.w]: [gauss-legendre points, gauss-legendre weights]
        '''
        
        [self.x, self.w] = np.polynomial.legendre.leggauss(self.nt)
        return self.x, self.w

    def cheby_wr(self):

        '''
        chebyshev weights for quadrature
        : return self.wr: chebyshev quadrature weights
        '''
        
        w_r = np.pi/self.nr * np.ones(self.nr + 1)
        w_r[0] = np.pi / (2. * self.nr)
        w_r[-1] = np.pi / (2. * self.nr)
        self.wr = w_r
        
        return self.wr
        

    def theta_lat_grids(self):

        '''
        theta grid points in radials and colatitude, useful when making plots
        : return self.theta, lats_gl: theta points, colatitude points
        '''
        
        [xl, w] = self.gauss_leg_xw()
        self.theta = np.arccos(xl)
        theta_lat = np.pi/2 - self.theta
        lats_gl = theta_lat * 180./np.pi

        return self.theta, lats_gl

    def T(self, n, x):
        ''' cheby polys, abs(x) <= 1 '''
        
        return np.cos(n * np.arccos(x))

    def cheby_diff_matrix(self):
        
        '''
        chebyshev differentiation matrix for taking derivatives on the r-grid
        : return DN: chebyshev differentiation matrix
        '''
        
        self.rcs = self.cheby_shift()
        DN = np.zeros([self.nr+1, self.nr+1])

        # i = row, j = column
        # off-diagonal terms
        for i in range(self.nr + 1):
            for j in range(self.nr + 1):

                # off-diagonal entries  
                if i != j:

                    if i == 0 or i == self.nr:
                        ci = 2.
                    else:
                        ci = 1.
                    if j == 0 or j == self.nr:
                        cj = 2.
                    else:
                        cj = 1. 

                    d = (ci / cj) * (-1.)**(i + j) / (self.rcs[i] - self.rcs[j])
                    DN[i, j] = d

        # diagonal terms
        for i in range(self.nr + 1):
            DN[i, i] = -1. * sum(DN[i, :])

        return DN

    def leg_P(self):

        '''
        legendre polynomials up to degree self.leg_deg
        : return P_matrix: matrix of legendre polynomials evaluated at chebyshev points (xl, or cos(theta));
        :                  shape (maximum degree, nt), that is, the first row is the legendre polynomial of
        :                  degree zero evaluated at each chebyshev point on the grid 
        '''
        
        xl, w = self.gauss_leg_xw() 
        
        P_matrix = np.zeros([self.leg_deg + 1, np.size(xl)])

        for ii in range(self.leg_deg + 1):
            if ii == 0:
                P_matrix[0, :] = np.ones(np.size(xl))
            if ii == 1:
                P_matrix[1, :] = xl
            if ii > 1:
                P_matrix[ii, :] = (2. * ii - 1.)/ii * xl * P_matrix[ii-1, :] - (ii - 1.)/ii * P_matrix[ii-2, :]
    
        return P_matrix
    
    def leg_P_deriv(self):

        '''
        legendre polynomial derivatives
        : return dP_dz: legendre derivatives with respect to chebyshev points, xl [-1, 1], evaluated at
        :               each chebyshev point on the grid; shape (maximum degree, number of theta grid points)
        : return dP_dtheta: legendre derivatives with respect to theta points [0, pi], evaluated at each
        :                chebyshev point on the grid; shape (maximum degree, nt)
        '''
        
        xl, w = self.gauss_leg_xw()
        theta, lats_gl = self.theta_lat_grids()
        
        P_matrix = self.leg_P()
        dP_dz = np.zeros([self.leg_deg + 1, self.nt])

        for ii in range(self.leg_deg + 1):
            if ii == 0:
                dP_dz[ii, :] = np.zeros(self.nt)
            else:
                dP_dz[ii, :] = ii/(xl**2 - 1.) * (xl*P_matrix[ii, :] - P_matrix[ii-1, :])

        dP_dtheta = -np.sin(theta) * dP_dz 

        return dP_dz, dP_dtheta
