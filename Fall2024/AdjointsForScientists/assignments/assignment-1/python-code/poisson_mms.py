import numpy as np
import mfem.ser as mfem
import matplotlib.pyplot as plt

def calc_gamma(x: np.float64, y: np.float64):
    return (np.pi/(np.e - 1.0)) * np.exp(x)

def calc_uexact(x: np.float64, y: np.float64):
    val = (np.pi/(np.e - 1.0)) * (np.exp(x) - 1.0)
    uex = np.exp(y) * np.sin(val)
    return uex

def calc_grad_uexact(x: np.float64, y: np.float64):
    val = (np.pi/(np.e - 1.0)) * (np.exp(x) - 1.0)
    dudx = np.exp(y) * np.cos(val) * calc_gamma(x, y)
    dudy = calc_uexact(x, y)
    return np.array([dudx, dudy])

def calc_beta(x: np.float64, y: np.float64, bOption: np.int64=1):
    if bOption == 1:
        beta = 1.0
    else:
        beta = np.pi**2.0 * (np.exp(x)-1.0) * (np.e - np.exp(x))/ (np.e - 1.0)**2.0
    return beta

def calc_load(x: np.float64, y: np.float64):
    u = calc_uexact(x, y)
    G = calc_gamma(x, y)
    Du= calc_grad_uexact(x, y)
    return -u * G * (1.0 - G**2.0) - 2.0 * G * Du[0]


class ForcingCoefficient(mfem.PyCoefficient):
    def __init__(self, force, dim):
        self.dim = dim
        self.force = force
        super(ForcingCoefficient, self).__init__()
    
    def EvalValue(self, p):
        assert p.size==self.dim, "MMS solution is not 2D"
        x = p[0]
        y = p[1]
        return self.force(x, y)


class PoissonMMS:
    def __init__(self, thermalOpt: dict):
        # set-up mesh
        self.nx = thermalOpt['mesh']['nx']
        self.ny = thermalOpt['mesh']['ny']
        self.xd = thermalOpt['mesh']['xd']
        self.yd = thermalOpt['mesh']['yd']
        
        # static-condensation
        self.sc = thermalOpt['fe-disc']['static-condensation']
        
        self.mesh = mfem.Mesh(mfem.Mesh.MakeCartesian2D(self.nx, self.ny, mfem.Element.QUADRILATERAL, sx=self.xd, sy=self.yd))
        self.mesh.EnsureNodes()
        
        # set-up mfem fe_col and fespace
        self.dim   = 2
        self.order = np.int64(thermalOpt['fe-disc']['order'])
        self.h1  = thermalOpt['fe-disc']['h1']
        if self.h1:
            print('H1 continuous Finite Elements are used')
            self.fe_col = mfem.H1_FECollection(self.order, self.dim)
        else:
            print('L2 discontinuous Finite Elements are used')
            self.fe_col = mfem.DG_FECollection(self.order, self.dim)
        self.fespace= mfem.FiniteElementSpace(self.mesh, self.fe_col)
        
        size = self.fespace.GetTrueVSize()
        print("Number of finite element unknowns: " + str(size))
        
        # set-up solution vector
        self.u_gf = mfem.GridFunction(self.fespace)
        self.u    = np.zeros( self.fespace.GetTrueVSize(), dtype=np.float64 )
        
        # set-up dirchlet boundary condition coefficient
        self.dbcCoef = ForcingCoefficient(calc_uexact, self.dim)
        
        # set-up forcing function
        self.force = ForcingCoefficient(calc_load, self.dim)
        
        # set-up gamma coefficient
        self.gamma = ForcingCoefficient(calc_gamma, self.dim)
        
        self.maxiter = thermalOpt['linear-solve']['max-iter']
        self.flux_bdr_marker = np.array(thermalOpt['objective']['flux-bdr-marker'], dtype=np.int64)
        self.bOption = thermalOpt['objective']['bOption']
        self.beta = lambda x, y: calc_beta(x, y, self.bOption)
        
        # self.kappa = (self.order+1.0)**2.0
        self.kappa = thermalOpt['fe-disc']['kappa']
        self.sigma = -1.0
    
    def solveForState(self):
        # 1. Marking the dirchlet boundary arrays
        dbc_bdr = mfem.intArray(self.mesh.bdr_attributes.Max())
        dbc_bdr.Assign(1)
        
        # 2. Assign grid-function to 0.0 everywhere
        self.u_gf.Assign(0.0)
        
        # 3. Create bilinear form
        a = mfem.BilinearForm(self.fespace)    
        # 3.1 Add Domain Diffusion integrator
        a.AddDomainIntegrator(mfem.DiffusionIntegrator(self.gamma))
        if self.h1 !=1 :
            # 3.2 Add the interfacial portion of the Laplace operator if L2 elements are used
            a.AddInteriorFaceIntegrator(mfem.DGDiffusionIntegrator(self.gamma,
                                                               self.sigma, self.kappa))
        # 3.3 Counteract the n.Grad(u) term on the Dirichlet portion of the boundary
        a.AddBdrFaceIntegrator(mfem.DGDiffusionIntegrator(self.gamma, self.sigma, self.kappa),
                               dbc_bdr)
        # 3.4 Assemble bilinear form
        a.Assemble()
        
        # 4. Create linear form
        b = mfem.LinearForm(self.fespace)
        # 4.1 Add the forcing linear form integrator
        b.AddDomainIntegrator(mfem.DomainLFIntegrator(self.force))
        # 4.2 Add the desired value for the Dirichlet boundary
        b.AddBdrFaceIntegrator(mfem.DGDirichletLFIntegrator(self.dbcCoef, self.gamma,
                                                            self.sigma, self.kappa), dbc_bdr)
        b.Assemble()
        

        # 5. Solving Linear System to find the grid-function
        A = mfem.OperatorPtr()
        B = mfem.Vector()
        X = mfem.Vector()

        ess_tdof_list = mfem.intArray()
        a.FormLinearSystem(ess_tdof_list, self.u_gf, b, A, X, B)

        print("Size of linear system: " + str(A.Height()))
        print("Size of grid function: " + str(self.u_gf.Size()))
        
        # Solve the linear system A X = B.
        AA = mfem.OperatorHandle2SparseMatrix(A)
        
        # numpy solving    
        Adense = np.array( AA.ToDenseMatrix().GetDataArray() )
        bdense = np.array( B.GetDataArray() )
        self.u = np.linalg.solve(Adense, bdense)
        self.u_gf.Assign(self.u)
        
        self.Lstar = Adense.transpose()
        
    
    def calcL2Error(self):
        err = self.u_gf.ComputeL2Error(self.dbcCoef)
        return err
    
    def compute_objective(self):
        flx_bdr = mfem.intArray(self.mesh.bdr_attributes.Max())
        flx_bdr[0] = 1
        
        # Set-up coefficients        
        beta = ForcingCoefficient(self.beta, self.dim)
        corr = mfem.ProductCoefficient(beta, self.dbcCoef)
        
        lf = mfem.LinearForm(self.fespace)
        lf.AddBdrFaceIntegrator(mfem.DGDirichletLFIntegrator(beta, self.gamma, -self.sigma, -self.kappa), flx_bdr)
        lf.Assemble()
        gh = lf.GetDataArray()
        
        ones = np.ones_like(self.u)
        lf2 = mfem.LinearForm(self.fespace)
        lf2.AddBdrFaceIntegrator(mfem.DGDirichletLFIntegrator(corr, self.gamma, 0.0, self.kappa), flx_bdr)
        lf2.Assemble()
        
        cval = np.dot(ones, lf2.GetDataArray())
        
        J = np.dot(gh, self.u) + cval
        return J, gh
    
    def solveForAdjoints(self):
        dJdu = self.compute_objective()[1]
        self.psi = np.linalg.solve(self.Lstar, -dJdu)
    
    def plotAdjoints(self, fsave: bool=False):
        
        nodeVec = mfem.Vector(self.mesh.GetNV() * self.mesh.Dimension())
        self.mesh.GetNodes(nodeVec)
        nodes = ( np.array( nodeVec.GetDataArray(), dtype=np.float64 ) )
        nodes = nodes.reshape([ self.mesh.GetNV(), self.mesh.Dimension() ])
        
        psi_nodeVec = mfem.Vector()
        psi_gf = mfem.GridFunction(self.fespace, mfem.Vector(self.psi), 0)
        psi_gf.GetNodalValues( psi_nodeVec, 1 )
        psi_nodes= np.array( psi_nodeVec.GetDataArray(), dtype=np.float64 )
        
        x = np.linspace(0.0, self.xd, self.nx+1)
        y = np.linspace(0.0, self.yd, self.ny+1)
        X, Y = np.meshgrid(x, y)
        U = psi_nodes.reshape([self.nx+1, self.ny+1])
        
        fig1, ax2 = plt.subplots(layout='constrained')
        CS = ax2.contourf(X, Y, U, 15, cmap=plt.cm.bone)
        CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='r')

        ax2.set_title(r'$\psi$ adjoint var contours')
        ax2.set_xlabel(r'$x$')
        ax2.set_ylabel(r'$y$')

        # Make a colorbar for the ContourSet returned by the contourf call.
        cbar = fig1.colorbar(CS)
        cbar.ax.set_ylabel(r'$\psi$')
        # Add the contour line levels to the colorbar
        cbar.add_lines(CS2)
        
        if fsave:
            fname = 'adj_bOption'+str(self.bOption)+'.pdf'
            plt.savefig(fname, transparent=True, bbox_inches='tight')
        plt.show()
