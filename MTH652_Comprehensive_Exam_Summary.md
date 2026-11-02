# MTH 652 Advanced Numerical Analysis: Comprehensive Exam Summary
Based on Lectures 1-14 (Pazner, Winter 2025), the lecture notebooks, and Homework 1-4. Proofs are skipped; key ideas and results are kept. Notation: $(u,v)=\int_\Omega uv\,dx$, $a(u,v)=\int_\Omega \nabla u\cdot\nabla v\,dx$, $M$ = mass matrix, $A$ or $K$ = stiffness matrix, $V_h$ = finite element space.
# Part I: Parabolic Problems and Time Discretization (Lectures 1-3)
## 1. Classification of second-order PDEs (L1)
- Operator $Lu=\sum_{i,j} a_{ij}\,\partial_{x_i}\partial_{x_j}u$ with symmetric coefficient matrix $A=(a_{ij})$.
- Elliptic: $A$ definite. Parabolic: $A$ semidefinite (one zero eigenvalue). Hyperbolic: $A$ indefinite and nonsingular.
- For $n=2$, $Lu=au_{xx}+bu_{xy}+cu_{yy}$ and $\Delta=b^2-4ac=-4\det A$. Elliptic iff $\Delta<0$, parabolic iff $\Delta=0$, hyperbolic iff $\Delta>0$.
- Prototypes: Laplace/Poisson (elliptic), heat (parabolic, infinite speed of propagation), wave (hyperbolic, finite speed of propagation).
## 2. Heat equation and semi-discretization (L1)
- Model problem: $u_t-\Delta u=f$ in $\Omega\times[0,T]$, $u(x,0)=u_0$, $u=0$ on $\partial\Omega$.
- Stability for $f=0$: $\|u(t)\|_{L^2}\le\|u_0\|_{L^2}$ and $\|u_t(t)\|_{L^2}\le \frac{C}{t}\|u_0\|_{L^2}$.
- Method of lines (semi-discretization): discretize in space first (FEM), then apply an ODE solver in time.
- Semi-discrete weak form: find $u_h(t)\in V_h$ with $(\partial_t u_h,v_h)+a(u_h,v_h)=(f,v_h)$ for all $v_h\in V_h$.
- Matrix form with $u_h=\sum_j u_j(t)\phi_j$: $$M\frac{du}{dt}+Au=f,\qquad M_{ij}=(\phi_i,\phi_j),\quad A_{ij}=a(\phi_i,\phi_j).$$
- Semi-discrete energy stability (take $v_h=u_h$, $f=0$): $\frac12\frac{d}{dt}\|u_h\|^2=-\|\nabla u_h\|^2\le0$.
- With $M=E^TE$ and $\eta=Eu$: $\eta'=-\bar A\eta$ where $\bar A=E^{-T}AE^{-1}$ is SPD. Solution is $\eta(t)=e^{-\bar A t}\eta_0$ and $\|e^{-\bar A t}\|\le1$. Also $\|\bar A\|\lesssim h^{-2}$.
- Matrix exponential: $e^{B}=\sum_{k\ge0}B^k/k!$.
- Semi-discrete error estimate (Thm, duality argument): $$\max_{t\le T}\|u-u_h\|_{L^2}\le C\Big(1+\log\tfrac{T}{h^2}\Big)\max_{t\le T}h^2\|u(t)\|_{H^2}.$$
- Proof ideas: elliptic (Ritz) projection $\tilde u_h$ with $a(\tilde u_h,v_h)=a(u,v_h)$. Split error as $u-u_h=(u-\tilde u_h)+(\tilde u_h-u_h)$. Use a backward-in-time dual problem and the bound $\int_0^{t_f}\|\partial_t w_h\|\,dt\le C(1+\log\frac{t_f}{h^2})\|w_f\|$.
- Qualitative behavior (1D, $u_0=\sum a_j\sin(jx)$): $u=\sum a_j e^{-j^2t}\sin(jx)$. High frequencies decay fastest, so the solution smooths in time. This motivates fine meshes and small steps early (initial transients), and coarse meshes and large steps later.
## 3. Time integration for the model ODE $u'=Bu$, $B$ negative definite (L2)
Exact solution $u(t)=e^{Bt}u_0$ with $\|u\|$ decreasing. We want the discrete analogue $\|u_{i+1}\|\le\|u_i\|$.
- Forward Euler (explicit, 1st order): $u_{i+1}=(I+\Delta tB)u_i$.
- Forward Euler stability needs $|1+\Delta t\lambda|<1$ for all eigenvalues $\lambda$ of $B$, i.e. $\Delta t<2/|\lambda|_{\max}$. This is conditional stability.
- For the heat equation $\lambda_{\max}\sim h^{-2}$, so $\Delta t\le Ch^2$ (CFL condition, Courant-Friedrichs-Lewy 1928). Each uniform refinement of the mesh costs 4 times as many time steps. Need $\Delta t\to0$ like $h^2$ for convergence.
- Backward Euler (implicit, 1st order): $(I-\Delta tB)u_{i+1}=u_i$. Eigenvalues of the amplification matrix are $\frac{1}{1-\Delta t\lambda}$, which have modulus $<1$ always. This is unconditionally stable.
- The notebook `forward_euler_backward_euler` uses $u'=-u$, exact solution $e^{-t}$. Forward Euler blows up if $\Delta t>2$. Backward Euler is stable for any $\Delta t$.
- FEM versions (with $M$): $$\text{FE: } u_{i+1}=u_i-\Delta t\,M^{-1}Au_i,\qquad \text{BE: } (M+\Delta tA)u_{i+1}=Mu_i.$$
- Variational forms: FE is $(u_{i+1},v)=(u_i,v)-\Delta t\,a(u_i,v)$. BE is $(u_{i+1},v)=(u_i,v)-\Delta t\,a(u_{i+1},v)$.
- FE needs a mass-matrix solve ($\kappa(M)=O(1)$, easy for CG). BE needs a solve with $M+\Delta tA$, with $\kappa(M+\Delta tA)\approx\Delta t/h^2$. For large $\Delta t$ this is ill-conditioned, so a good preconditioner is needed.
- Theta methods: $\frac{u_{i+1}-u_i}{\Delta t}=\theta Bu_i+(1-\theta)Bu_{i+1}$ (this course's convention: $\theta=1$ is FE, $\theta=0$ is BE).
- Update: $u_{i+1}=(I-(1-\theta)\Delta tB)^{-1}(I+\theta\Delta tB)u_i$. Unconditionally stable iff $\theta\le\frac12$.
- Crank-Nicolson ($\theta=\frac12$): unconditionally stable, 2nd order (trapezoidal rule analogue). All other $\theta$ are 1st order.
## 4. Runge-Kutta methods (L3)
- For $u'=f(t,u)$, an $s$-stage RK method is $$k_j=f\Big(t_i+c_j\Delta t,\ u_i+\Delta t\sum_{\ell}a_{j\ell}k_\ell\Big),\qquad u_{i+1}=u_i+\Delta t\sum_j b_jk_j.$$
- Butcher tableau: $c$ (nodes/abscissas), $A=(a_{ij})$ (RK matrix), $b$ (weights).
- $A$ strictly lower triangular gives an explicit method. $A$ lower triangular gives diagonally implicit (DIRK), with $N\times N$ solves. Otherwise fully implicit, with a coupled $sN\times sN$ system.
- Explicit midpoint (2nd order): $k_1=f(t_i,u_i)$, $k_2=f(t_i+\frac{\Delta t}{2},u_i+\frac{\Delta t}{2}k_1)$, $u_{i+1}=u_i+\Delta t\,k_2$.
- RK4 (4th order, 4 stages): $u_{i+1}=u_i+\frac{\Delta t}{6}(k_1+2k_2+2k_3+k_4)$ with $k_2,k_3$ at the midpoint and $k_4$ at $t_i+\Delta t$.
- Even RK4 is explicit and so conditionally stable. Implicit methods are preferred for parabolic (stiff) problems.
- Link to quadrature: implicit RK methods are generated by quadrature rules (CN is the trapezoidal rule).
## 5. Discontinuous Galerkin in time (L3)
- Space $W_h=\{v:[0,T]\to V_h,\ v|_{I_i}\in P_q(I_i)\}$, with partition $0=t_0<\dots<t_K=T$. No continuity across $t_i$.
- Traces $v^\pm(t_i)$ and jump $[\![v(t_i)]\!]=v^+(t_i)-v^-(t_i)$.
- Upwind numerical flux: $\hat u_h(t_i)=u_h^-(t_i)$ (information flows forward in time), and $\hat u_h(t_0)=u_0$.
- Interval-by-interval problem: find $u_h\in P_q(I_i)$ such that for all $v_h\in P_q(I_i)$ $$-\int_{I_i}(u_h,\partial_tv_h)\,dt+(u_h^-(t_i),v_h^-(t_i))+\int_{I_i}a(u_h,v_h)\,dt=\int_{I_i}(f,v_h)\,dt+(u_h^-(t_{i-1}),v_h^+(t_{i-1})).$$
- $q=0$ gives a modified backward Euler: $(u^i,v)+\Delta t_i\,a(u^i,v)=(u^{i-1},v)+\int_{I_i}(f,v)\,dt$. The only change is that $f(t_i)$ is replaced by its time average.
- $q=1$ with trapezoidal quadrature gives a block system equivalent to the fully implicit RK method Lobatto IIIC, with tableau $\begin{pmatrix}0&\frac12&-\frac12\\1&\frac12&\frac12\\&\frac12&\frac12\end{pmatrix}$.
- In general DG in time of order $q$ corresponds to fully implicit RK.
## 6. Adaptive time step and mesh control (L3)
- Error estimate for the $q=0$ scheme: $$\|u(t_n)-u_h^n\|\le C\Big(\max_{i\le n}\Delta t_i\|\partial_tu\|_{\infty,I_i}+\max_t h^2\|u(t)\|_{H^2}\Big),$$ where $C$ contains a slowly growing $\sqrt{1+\log(t_n/\Delta t_n)}$. The first term is the temporal error and the second is the spatial error.
- Computable proxy: $\Delta t_i\|\partial_tu\|_{\infty,I_i}\approx\|u_h^i-u_h^{i-1}\|$. Choose $\Delta t_i$ by trial and error so that $\|u_h^i-u_h^{i-1}\|\approx\delta/C$. Accept the provisional step if the tolerance is met, otherwise rescale and retry.
- Spatial error: use $\|u\|_{H^2}\lesssim\|u_t\|_{L^2}$ (from the PDE with $f=0$ and elliptic regularity). Choose the mesh size $h_i$ such that $h_i^2\|u_h^i-u_h^{i-1}\|\le\delta\Delta t_i/C$.
## 7. Homework 1 highlights (heat equation, FE/BE/theta)
- Test problem: $f=0$, $u_0=\sin(\pi x)\sin(\pi y)$ on $[0,1]^2$. Exact solution $u(x,y,t)=e^{-2\pi^2t}\sin(\pi x)\sin(\pi y)$.
- Forward Euler is stable only under the CFL restriction $\Delta t\lesssim h^2$. The code returns NaN if the norm grows by more than 10 percent.
- Expected convergence under refinement (halving both $h$ and $\Delta t$): FE and BE are first order in time, CN is second order in time. Note FE is stable only if $\Delta t\lesssim h^2$, so halving $\Delta t$ and $h$ together eventually becomes unstable.
# Part II: Linear Elasticity (Lectures 4-5)
## 8. Kinematics and equilibrium (L4)
- Deformation $\Phi(x)=x+u(x)$ with displacement $u$. Deformation gradient $F=D\Phi=I+Du$.
- Rigid motion $\Phi(x)=x_1+Q(x-x_0)$ with $Q^TQ=I$. Displacement is not deformation (rigid motions have $u\ne0$ but no deformation).
- Green-Lagrange strain $E=\frac12(F^TF-I)=\frac12(Du+Du^T+Du^TDu)$. It vanishes for rigid motions.
- Linearization for small $\|Du\|$: infinitesimal strain $$\varepsilon(u)=\nabla^su=\tfrac12(Du+Du^T).$$ Diagonal entries are normal strains (stretching), off-diagonal entries are shear strains.
- Kernel of $\nabla^s$: translations plus infinitesimal rotations, $w(x)=c+Zx$ with $Z$ antisymmetric. These are the rigid body modes $RM$. In 2D: $w=c+b(y,-x)^T$, dimension 3. In 3D: dimension 6. Planar rotation by $\theta$ gives $\varepsilon\ne0$ at second order, which vanishes only for linearized rotations $\begin{pmatrix}1&-\theta\\\theta&1\end{pmatrix}$.
- Stress tensor $\sigma$ (symmetric by Cauchy's theorem). Traction on a surface with normal $n$ is $\sigma n$.
- Force balance on arbitrary subsets $\omega$ gives Cauchy's equilibrium equation $$-\nabla\cdot\sigma=f.$$ In $\mathbb R^3$ this is 3 equations in 6 stress unknowns, so a constitutive law is needed.
- Hooke's law (isotropic, linear): $$\sigma=2\mu\,\varepsilon(u)+\lambda\,\mathrm{tr}(\varepsilon)I=2\mu\nabla^su+\lambda(\nabla\cdot u)I.$$
- Lame parameters: $\mu=\frac{E}{2(1+\nu)}$, $\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}$ ($E$ = Young's modulus, $\nu$ = Poisson ratio). Nearly incompressible means $\nu\to\frac12$, so $\lambda\to\infty$.
- Boundary conditions: displacement (Dirichlet) $u=g_D$ on $\Gamma_D$, traction (Neumann) $\sigma n=g_N$ on $\Gamma_T=\partial\Omega\setminus\Gamma_D$.
## 9. Weak formulation (L4)
- Space $V=\{v\in[H^1(\Omega)]^d: v|_{\Gamma_D}=0\}$.
- Find $u$ with $u|_{\Gamma_D}=g$ and $a(u,v)=\int_\Omega f\cdot v\,dx+\int_{\Gamma_T}t\cdot v\,ds$ where $$a(u,v)=\int_\Omega\big(2\mu\,\nabla^su:\nabla^sv+\lambda(\nabla\cdot u)(\nabla\cdot v)\big)dx.$$
- It is symmetric because $\nabla^su:\nabla v=\nabla^su:\nabla^sv$ (symmetric:antisymmetric contraction is zero).
## 10. Well-posedness (L5)
- Continuity of $a$ is easy (Cauchy-Schwarz). Coercivity is the hard part: need $\|\nabla^sv\|_{L^2}\gtrsim\|v\|_{H^1}$ (Korn-type inequality).
- Lemma (right inverse of div): for $p\in L^2(\Omega)$ there exists $v\in H^1$ with $\nabla\cdot v=p$ and $\|v\|_{H^1}\le C\|p\|_{L^2}$. If $\int p=0$, $v$ can be chosen in $H^1_0$. Idea: solve $-\Delta w=p$ and set $v=-\nabla w$. For zero trace use the Neumann problem plus a $\mathrm{curl}\,\psi$ correction.
- Define $\widehat H^1=\{v:\int v=0,\ \int\mathrm{rot}\,v=0\}$ (removes rigid motions), with $\mathrm{rot}\,v=-\partial_yv_1+\partial_xv_2$.
- Korn on $\widehat H^1$: $\|\nabla^sv\|_{L^2}\ge C\|v\|_{H^1}$ for $v\in\widehat H^1$.
- Korn's second inequality (with lower-order term): $\|\nabla^sv\|_{L^2}+\|v\|_{L^2}\ge\alpha\|v\|_{H^1}$ for all $v\in H^1$.
- Decomposition: $H^1=\widehat H^1\oplus RM$. Proof of Korn II uses decomposition, Korn on $\widehat H^1$, compactness, and finite dimensionality of $RM$.
- If $|\Gamma_D|>0$, the only rigid motion vanishing on $\Gamma_D$ is zero, so $\|\nabla^sv\|\ge C\|v\|_{H^1}$ on $V$, so $a$ is coercive and the problem has a unique solution.
- Pure traction problem: solvable iff the compatibility condition $\int_\Omega f\cdot v+\int_{\Gamma_T}t\cdot v=0$ holds for all $v\in RM$. Then the solution is unique in $\widehat H^1$ (unique up to rigid motions in general).
- FEM: $V_h=[V_h]^d$ (vector-valued Lagrange). Galerkin gives $\|u-u_h\|_{H^1}\le Ch\|u\|_{H^2}$ for P1, but $C$ depends on $\mu,\lambda$ and blows up as $\lambda\to\infty$.
## 11. Volumetric locking (L5, HW2)
- Stability with $\mu=1$, pure displacement BC: $\|u\|_{H^2}+\lambda\|\nabla\cdot u\|_{H^1}\le C\|f\|_{L^2}$ with $C$ independent of $\lambda$. So $\nabla\cdot u\to0$ as $\lambda\to\infty$.
- Heuristic explanation: the energy forces $\nabla\cdot u_h\approx0$. A continuous P1 vector field with zero boundary values and zero divergence must be $0$, so $u_h\to0$ while $u\not\to0$. The scheme "locks".
- Remedies: higher polynomial degree (nontrivial divergence-free fields exist, analysis is hard), mixed methods (discretize $\sigma_h$ and $u_h$ together, giving a non-SPD saddle system with better stress approximations), or reduced/selective integration.
## 12. Homework 2 highlights
- Nullspace of $\nabla^s$ is the set $\{c+Zx: Z\ \text{antisymmetric}\}$ (infinitesimal rigid motions, dimension $\binom{d+1}{2}$: 3 in 2D, 6 in 3D). Related to rotations: $Q=e^{Z}\approx I+Z$ for small rotations.
- In 2D, if $w\in RM$ vanishes on two distinct points of $\Gamma_D$ then $w\equiv0$ (so any $\Gamma_D$ with nonempty open set suffices for uniqueness).
- Decomposition $H^1=\widehat H^1+RM$: take $b=-\frac{1}{2|\Omega|}\int\mathrm{rot}\,v$, $c=\frac{1}{|\Omega|}\int(v-b(y,-x)^T)$, $w=c+b(y,-x)^T$, and $v-w\in\widehat H^1$.
- Divergence-free P1 (lecture claim): the only divergence-free field in the continuous piecewise-linear space with zero boundary values is $v_h=0$. This is the mechanism behind locking. Without boundary conditions the space can contain nonzero fields (for example constants and infinitesimal rotations).
- Manufactured solution $u=(-\cos\pi x\sin\pi y,\ \sin\pi x\cos\pi y)$ is divergence-free, $\nabla\cdot u=0$, and the corresponding right-hand side is $f=2\mu\pi^2\,u$, independent of $\lambda$.
- Locking experiment ($\lambda\in\{1,10^2,10^5\}$): for small $\lambda$ the $L^2$ error converges at the optimal rate $O(h^2)$. For large $\lambda$ the convergence rate degrades or stalls (locking). CG iterations grow with $n_x$ (condition number $\sim h^{-2}$) and also with $\lambda$ (the stiffness matrix becomes ill-conditioned because of the penalty term).
# Part III: Non-Symmetric Problems and Krylov Methods (Lectures 6-8)
## 13. Symmetric theory recap and Lax-Milgram (L6)
- Variational problem: find $u\in V$ with $a(u,v)=F(v)$ for all $v\in V$.
- Conditions: (i) symmetry, (ii) continuity $a(u,v)\le\gamma\|u\|\|v\|$, (iii) coercivity $a(u,u)\ge\alpha\|u\|^2$. With (i)-(iii) and $F\in V'$, Riesz representation gives existence and uniqueness.
- Lax-Milgram: (ii) and (iii) suffice (drop symmetry). Unique $u\in V$ with $\|u\|_V\le\|F\|_{V'}/\alpha$. Proof idea: Riesz map $\tau$, operator $T(v)=v-\rho(\tau(Av)-\tau(F))$ is a contraction for $0<\rho<2\alpha/\gamma^2$.
- Cea's Lemma (non-symmetric): $\|u-u_h\|_V\le\frac{\gamma}{\alpha}\min_{v_h\in V_h}\|u-v_h\|_V$. Galerkin orthogonality $a(u-u_h,v_h)=0$ still holds. In the symmetric case the constant improves to $\sqrt{\gamma/\alpha}$ (error is minimal in the energy norm).
## 14. Convection-diffusion-reaction (L6)
- Equation: $-\mu\Delta u+\beta\cdot\nabla u+\varepsilon u=f$, $u=0$ on $\partial\Omega$. Weak form: $\mu(\nabla u,\nabla v)+(\beta\cdot\nabla u,v)+\varepsilon(u,v)=(f,v)$.
- The advection term is non-symmetric (in fact skew-symmetric when $\nabla\cdot\beta=0$): $(\beta\cdot\nabla u,u)=0$. Hence $a(u,u)\gtrsim\mu\|\nabla u\|^2+\varepsilon\|u\|^2$ and coercivity holds for constant or divergence-free $\beta$. Well-posed by Lax-Milgram.
- Convection-dominated regime ($\mu\ll1$): the coercivity constant $\alpha\sim\mu$ is small, so the Cea bound $\gamma/\alpha$ blows up and standard Galerkin FEM can fail (oscillations). This motivates stabilization and DG methods.
- Linear algebra: the matrix $A_{ij}=a(\phi_j,\phi_i)$ is not symmetric, so we cannot use CG or Cholesky. Use LU for direct solves and GMRES for iterative solves.
## 15. Krylov subspaces and CG recap (L7)
- Krylov subspace: $\mathcal K_m(A,u)=\mathrm{span}\{u,Au,\dots,A^{m-1}u\}$, nested and with dimension at most $m$.
- Which method: SPD gives CG, symmetric indefinite gives MINRES, general non-symmetric gives GMRES.
- CG: iterate $x^{(i)}$ minimizes $\|x-x^{(i)}\|_A$ over $x^{(0)}+\mathcal K_i(A,r^{(0)})$. Residuals and search directions lie in $\mathcal K_{i+1}$.
- Naive least squares in $\mathcal K_m$ (basis $K_m$) stagnates: columns $A^ku$ align with the dominant eigenvector (power method, $A^kx/\lambda_1^k\to v_1$), so $K_m$ becomes nearly singular.
## 16. Arnoldi and GMRES (L7, L8)
- Arnoldi iteration builds an orthonormal basis $q_1,\dots,q_m$ of $\mathcal K_m(A,u)$ via Gram-Schmidt: $$Aq_k=\sum_{i=1}^{k+1}H_{ik}q_i,\qquad AQ_m=Q_{m+1}H_m,$$ where $H_m\in\mathbb R^{(m+1)\times m}$ is upper Hessenberg ($H_{ij}=0$ for $i>j+1$).
- GMRES: $x_m=\arg\min_{y\in\mathcal K_m(A,b)}\|b-Ay\|$. Writing $x_m=Q_my_m$ and using $q_1=b/\|b\|$, $$\min_{y_m\in\mathbb R^m}\ \big\|\,\|b\|e_1-H_my_m\big\|.$$ This is a small $(m+1)\times m$ least-squares problem.
- Converges in at most $N$ iterations for any invertible $A$ (exact arithmetic).
- Numerical stability: classical Gram-Schmidt loses orthogonality. Modified Gram-Schmidt (update $z$ inside the loop) is better. Householder-based Arnoldi is the most stable.
- Householder reflector $P=I-2ww^T$ with $\|w\|=1$ is symmetric, orthogonal, $P=P^T=P^{-1}$. To zero out $x$ below the first entry choose $\alpha=-\mathrm{sign}(x_1)\|x\|$ and $$w=\frac{x+\mathrm{sign}(x_1)\|x\|e_1}{\|x+\mathrm{sign}(x_1)\|x\|e_1\|}.$$ This avoids cancellation. $P$ is never formed explicitly: $PK=K-2w(w^TK)$.
- Least-squares via Givens rotations $G(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix}$ gives a cheap incremental QR of $H_m$ (one new rotation per iteration).
- Cost: must store all Krylov vectors, so memory $O(m N)$ and work $O(m^2N+m\,\mathrm{nnz})$ after $m$ steps. Worst case $O(N^2)$ memory and $O(N^3)$ work. Goal: $n_{it}\ll N$, ideally $O(1)$, via preconditioning.
- Restarted GMRES(M): memory $O(MN)$, work $O(n_{it}MN)$, but restarting can slow or stall convergence, so $n_{it}$ increases.
- No simple convergence bound in terms of extremal eigenvalues for general $A$ (non-normal). Preconditioning: left ($MA$) or right ($AM$) with $M\approx A^{-1}$ (for example incomplete LU).
## 17. Lanczos and MINRES (L8)
- For symmetric $A$, $H_m'=Q_m^TAQ_m$ is symmetric and Hessenberg, so tridiagonal. The three-term recurrence is $$Aq_m=H_{m-1,m}q_{m-1}+H_{mm}q_m+H_{m+1,m}q_{m+1}.$$ This is the Lanczos iteration.
- MINRES is GMRES with Lanczos: same iterates in exact arithmetic, but per iteration $O(N+\mathrm{nnz})$ work, $O(N)$ memory, and no restarts.
## 18. Homework 3 highlights
- Householder Arnoldi: with $P_k=I-2w_kw_k^T$, $q_j=P_1P_2\cdots P_je_j$. The first $j$ columns of $P_1\cdots P_m$ equal those of $P_1\cdots P_j$. The last column of $H_{k-1}$, namely $h_{k-1}$, is the first $k$ entries of $P_kP_{k-1}\cdots P_1Aq_{k-1}$.
- Cost: each step applies up to $k$ reflectors, each $O(N)$, plus one matvec, so total $O(m^2N+m\,\mathrm{nnz}(A))$. Memory $O(mN)$ for the $w$ vectors, since the $P_j$ are never formed.
- Givens rotation to zero $b$: $\cos\theta=a/\sqrt{a^2+b^2}$, $\sin\theta=-b/\sqrt{a^2+b^2}$ (so that $G(\theta)(a,b)^T=(\sqrt{a^2+b^2},0)^T$). Apply successive rotations to the subdiagonal entries of the Hessenberg $H$ to get $H=QR$ with $R$ upper triangular. The GMRES least-squares problem $\min\|\beta e_1-H_my\|$ is solved by applying the same rotations to $\beta e_1$ and back-substituting with $R$. One new rotation per iteration gives $O(1)$ update.
- Symmetric $A$ with Lanczos: per-iteration cost $O(N+\mathrm{nnz})$ (orthogonalize against 2 previous vectors), memory $O(N)$ for a short recurrence, compared with $O(mN)$ per iteration and $O(mN)$ memory for GMRES. Restarting is unnecessary. This is the MINRES vs GMRES advantage.
- GMRES(m) vs CG on the FEM stiffness matrix ($nx=20$): CG is cheaper (short recurrence). GMRES with restart 10 converges slowly or stagnates, while restart 50 is better but costs more memory per iteration. Convergence of GMRES is not governed by the spectrum the way CG is.
# Part IV: Mixed Finite Elements and Saddle-Point Systems (Lectures 9-12)
## 19. Stokes equations (L9)
- $-\Delta u+\nabla p=f$, $\nabla\cdot u=0$, $u=0$ on $\partial\Omega$. Valid at low Reynolds number, with incompressible fluid.
- Two unknowns and two equations means different spaces for $u_h$ and $p_h$ (mixed FEM). Stable pairing is subtle. Other examples: Darcy flow, mixed elasticity (avoids locking).
## 20. Abstract saddle-point theory (L9, L10)
- Problem: $V,P$ Banach/Hilbert spaces, $a:V\times V\to\mathbb R$, $b:V\times P\to\mathbb R$, $F\in V'$, $G\in P'$. Find $(u,p)\in V\times P$ with $$a(u,v)+b(v,p)=F(v)\ \ \forall v,\qquad b(u,q)=G(q)\ \ \forall q.$$ Operator form: $\begin{pmatrix}A&B'\\B&0\end{pmatrix}\begin{pmatrix}u\\p\end{pmatrix}=\begin{pmatrix}F\\G\end{pmatrix}$.
- Kernel space: $Z=\{v\in V: b(v,q)=0\ \forall q\in P\}=\ker B$.
- Finite-dimensional intuition: $u\in\ker B$. Coercivity of $A$ on $Z$ gives a unique $u$ (Lax-Milgram on $Z$). Then $F-Au\in Z^\perp=\mathrm{range}(B^T)$, so $p$ exists. Uniqueness of $p$ needs $B^T$ injective.
- Inf-sup (LBB, Ladyzhenskaya-Babuska-Brezzi) condition: $$\inf_{0\ne q\in P}\sup_{0\ne v\in V}\frac{b(v,q)}{\|v\|_V\|q\|_P}\ge\beta>0.$$
- Equivalent forms: $B'$ is an isomorphism of $P$ onto $Z^\circ$ with $\|B'q\|\ge\beta\|q\|$. Also $B$ is an isomorphism of $Z^\perp$ onto $P'$ with $\|Bv\|\ge\beta\|v\|$ (so $B$ is surjective). Proof tool: Closed Range Theorem.
- Main theorem: if $a,b$ are continuous, $a$ is coercive on $Z$ (constant $\alpha$), and $b$ satisfies inf-sup (constant $\beta$), then there is a unique solution and $$\|u\|_V+\|p\|_P\le C(\|F\|_{V'}+\|G\|_{P'}),$$ with explicit bounds like $\|u\|\lesssim\frac1\alpha\|F\|+\frac{C_a}{\alpha\beta}\|G\|$ and $\|p\|\lesssim\frac{C_a}{\alpha\beta}\|F\|+\frac{C_a^2}{\alpha\beta^2}\|G\|$.
- Idea for $G\ne0$: pick $u_\xi\in Z^\perp$ with $Bu_\xi=G$ (exists by surjectivity of $B$, with $\|u_\xi\|\le\|G\|/\beta$), then solve the homogeneous problem for $u-u_\xi$.
## 21. Stokes is well-posed (L10)
- $a(u,v)=(\nabla u,\nabla v)$, $b(v,p)=-(\nabla\cdot v,p)$, $V=[H^1_0]^d$, $P=L^2_0(\Omega)$ (zero mean, since only $\nabla p$ appears).
- $a$ is coercive on all of $V$ (Poincare). The inf-sup condition follows from the lemma: for $q\in L^2_0$ there is $v\in H^1_0$ with $\nabla\cdot v=q$, $\|v\|_{H^1}\le C\|q\|$, giving $\beta=1/C$.
## 22. Discrete inf-sup and stable pairs (L10, L11)
- Choose $V_h\subset V$, $P_h\subset P$. Continuity and coercivity pass to subspaces for Stokes. The discrete inf-sup condition is the nontrivial requirement and must hold uniformly in $h$.
- Discrete inf-sup implies $B_h$ surjective, so $\dim V_h\ge\dim P_h$ is necessary.
- P1-P0 (continuous linear velocity, piecewise constant pressure) fails. Angle counting in 2D: $2\pi n_V<\pi n_T$ implies $2n_V<n_T$, so $B$ has more rows than columns and cannot be surjective. The saddle matrix is singular.
- Fortin operator: $\Pi_h:V\to V_h$ with (i) $b(\Pi_hv,q_h)=b(v,q_h)$ for all $q_h\in P_h$ and (ii) $\|\Pi_hv\|_V\le C_f\|v\|_V$. For Stokes, (i) says the divergence commutes with projections: $\nabla\cdot\Pi_h=Q_h\nabla\cdot$ where $Q_h$ is the $L^2$ projection onto $P_h$.
- Theorem: a Fortin operator implies the discrete inf-sup condition with constant $\beta/C_f$.
- Construction for P2-P0: need $\int_e v=\int_e\Pi_hv$ on every edge $e$. Use quadratic edge bubbles $b_e=6\phi_i\phi_j/|e|$ (so $\int_eb_e=1$) and set $\Pi_hv=\sum_e\big(\int_ev\big)b_e$. Stability by Cauchy-Schwarz, inverse and trace inequalities.
- Other stable pairs: Taylor-Hood $P_k$-$P_{k-1}$ for $k\ge2$, with $P_2$-$P_1$ the lowest order.
- Error estimate (quasi-optimality for mixed problems): $$\|u-u_h\|_V\le C_1\inf_{v_h}\|u-v_h\|_V+C_2\inf_{q_h}\|p-q_h\|_P,\qquad \|p-p_h\|_P\le C_3\inf_{v_h}\|u-v_h\|_V+C_4\inf_{q_h}\|p-q_h\|_P,$$ with constants depending on $C_a,C_b,\alpha,\beta$. Galerkin orthogonality: $a(u-u_h,v_h)+b(v_h,p-p_h)=0$ and $b(u-u_h,q_h)=0$.
## 23. Darcy flow and H(div) (L11)
- Darcy's law: $u=-\frac k\mu\nabla p$. Incompressibility: $\nabla\cdot u=0$. Together: $-\nabla\cdot(\frac k\mu\nabla p)=0$, a Poisson problem for $p$.
- Standard FEM does not give a locally conservative flux. Mixed methods do (flux balance on each control volume).
- First-order system: $u+\nabla p=f$, $\nabla\cdot u=g$, $p=0$ on $\partial\Omega$. Weak form with $a(u,v)=\int u\cdot v$, $b(v,p)=-\int(\nabla\cdot v)p$, $V=H(\mathrm{div},\Omega)=\{w\in L^2:\nabla\cdot w\in L^2\}$, $P=L^2(\Omega)$, norm $\|v\|^2_{H(\mathrm{div})}=\|v\|^2_{L^2}+\|\nabla\cdot v\|^2_{L^2}$.
- $H(\mathrm{div})$ functions have well-defined normal traces (normal continuity across element interfaces).
- Well-posedness: $Z=\{v:\nabla\cdot v=0\}$, on which $a(v,v)=\|v\|^2_{H(\mathrm{div})}$ (coercive with constant 1). Inf-sup holds via the same divergence right-inverse lemma.
- Discrete spaces: $H(\mathrm{div})$-conforming elements such as Raviart-Thomas.
## 24. Linear algebra of saddle-point systems (L11)
- System: $$\mathcal A=\begin{pmatrix}A&B^T\\B&0\end{pmatrix}\begin{pmatrix}u\\p\end{pmatrix}=\begin{pmatrix}f\\g\end{pmatrix}.$$ It is the stationarity condition of the Lagrangian $L(u,p)=\frac12u^TAu-f^Tu+(Bu-g)^Tp$ for the constrained minimization $\min\frac12u^TAu-f^Tu$ subject to $Bu=g$. Here $p$ is the Lagrange multiplier.
- Assumptions: $A$ is SPD ($n\times n$), $B$ has full row rank ($m\times n$). Then $\mathcal A$ is symmetric and invertible.
- Block factorization: $$\begin{pmatrix}A&B^T\\B&0\end{pmatrix}=\begin{pmatrix}I&0\\BA^{-1}&I\end{pmatrix}\begin{pmatrix}A&0\\0&S\end{pmatrix}\begin{pmatrix}I&A^{-1}B^T\\0&I\end{pmatrix},\qquad S=-BA^{-1}B^T.$$
- The Schur complement $S$ is negative definite. By Sylvester's Law of Inertia (congruent symmetric matrices have the same numbers of positive/zero/negative eigenvalues), $\mathcal A$ has $n$ positive and $m$ negative eigenvalues. So it is indefinite and CG cannot be used.
- Use MINRES (residual minimization over $\mathcal K_m(\mathcal A,b)$ with a Lanczos three-term recurrence, no restarts).
- SPD case of MINRES bound: $\|r_m\|\le2\big(\frac{\sqrt\kappa-1}{\sqrt\kappa+1}\big)^m\|r_0\|$. General $\kappa=|\lambda|_{\max}/|\lambda|_{\min}$.
## 25. MINRES convergence (L12)
- Polynomial characterization: $\|r_k\|=\min_{p\in\mathcal P_k,p(0)=1}\|p(\mathcal A)b\|\le\|b\|\min_p\max_j|p(\lambda_j)|$ since the eigenvectors are orthonormal.
- For indefinite spectrum in $[\mu_1,\mu_2]\cup[\nu_1,\nu_2]$ with $\mu_1<\mu_2<0<\nu_1<\nu_2$ and equal lengths $\mu_2-\mu_1=\nu_2-\nu_1$: $$\min_p\max_j|p(\lambda_j)|\le2\left(\frac{\sqrt{|\mu_1\nu_2|}-\sqrt{|\mu_2\nu_1|}}{\sqrt{|\mu_1\nu_2|}+\sqrt{|\mu_2\nu_1|}}\right)^{\lceil k/2\rceil}.$$
- Symmetric case ($\mu_1=-1$, $\nu_2=1$, $-\mu_2=\nu_1=\kappa^{-1}$): $2\big(\frac{\kappa-1}{\kappa+1}\big)^{\lceil k/2\rceil}$. So MINRES at step $2k$ with condition $\kappa$ behaves like CG at step $k$ with condition $\kappa^2$ (roughly twice the steps of CG).
## 26. Block-diagonal preconditioning (L12)
- Preconditioner must be SPD (to preserve symmetry): $M=EE^T$ and apply MINRES to $E^T\mathcal AE$ (same spectrum as $M\mathcal A$). Only the action of $M$ is needed in the algorithm. Equivalent to minimizing the residual in the $M$-norm.
- Ideal block-diagonal preconditioner: $$\mathcal D=\begin{pmatrix}A^{-1}&0\\0&-S^{-1}\end{pmatrix}.$$
- Theorem: the preconditioned operator $\mathcal D\mathcal A$ has only three distinct eigenvalues, $1,\ \frac{1+\sqrt5}{2},\ \frac{1-\sqrt5}{2}$ (roots of $\lambda^2-\lambda-1=0$). So MINRES converges in at most 3 iterations.
- Generalized system $\begin{pmatrix}A&B^T\\B&-C\end{pmatrix}$ with $C$ symmetric positive semidefinite, and $S=-(C+BA^{-1}B^T)$. Then $$\sigma(\mathcal D\mathcal A)\subset\Big[-1,\tfrac{1-\sqrt5}{2}\Big]\cup\Big[1,\tfrac{1+\sqrt5}{2}\Big],$$ which is independent of $h$ and $\kappa(A)$ (uniform, optimal).
- Practical version: replace $A^{-1}$ by a spectrally equivalent $\mathcal B_A$ and $-S^{-1}$ by $\mathcal B_S$ and use $\mathrm{blockdiag}(\mathcal B_A,\mathcal B_S)$. The preconditioner stays optimal, only constants change. So the indefinite problem reduces to two SPD preconditioning problems.
- Stokes Schur complement (Thm): with $M_p$ the pressure mass matrix, $$\beta^2\,q^TM_pq\le q^T(-S)q\le q^TM_pq,$$ where $\beta$ is the discrete inf-sup constant. So $-S$ is spectrally equivalent to $M_p$, hence $\mathcal B_S=M_p^{-1}$ (or even $\mathrm{diag}(M_p)^{-1}$ since $\kappa(M_p)=O(1)$). For $A$ (vector Laplacian) use multigrid or domain decomposition.
# Part V: Preconditioning SPD Systems and Domain Decomposition (Lectures 13-14)
## 27. Spectral equivalence and generalized eigenvalues (L13)
- Condition number: $\kappa(A)=\lambda_{\max}/\lambda_{\min}$. For FEM stiffness matrices $\kappa(A)=O(h^{-2})$.
- For SPD $A,M$ the following problems share the same real positive eigenvalues: $Ax=\lambda Mx$, $M^{-1}Ax=\lambda x$, $M^{-1/2}AM^{-1/2}x=\lambda x$, $A^{1/2}M^{-1}A^{1/2}x=\lambda x$. Also $\lambda_{\min}=\inf_x\frac{x^TAx}{x^TMx}$ and $\lambda_{\max}=\sup_x\frac{x^TAx}{x^TMx}$.
- Spectral equivalence: $cx^TMx\le x^TAx\le Cx^TMx$ for all $x$. Then $\sigma(M^{-1}A)\subset[c,C]$ and $\kappa(M^{-1}A)\le C/c$.
- Preconditioner goal: $\mathcal B=M^{-1}$ with (1) $A\sim M$ spectrally equivalent with $C/c$ small, (2) cheap action of $M^{-1}$. Extremes: $M=I$ (cheap, no improvement) and $M=A$ (perfect but as costly as the original problem).
- Uniform (parameter-independent) preconditioner: $c,C$ independent of $h$. Then $\kappa(\mathcal BA)=O(1)$ and PCG converges in a mesh-independent number of iterations.
- PCG: $M^{-1}A$ is self-adjoint and positive definite in the $M$-inner product $(x,y)_M=y^TMx$. PCG is CG in this inner product. Each iteration needs one application of $A$ and one of $\mathcal B=M^{-1}$ (never $M$).
- Change of inner product: if $T$ is SPD in the $A$-inner product $(x,y)_A=y^TAx$ then its eigenvalues are bounded by $\frac{(Tx,x)_A}{(x,x)_A}$.
## 28. Additive Schwarz / parallel subspace correction (L13)
- Decompose $V=\sum_{j=0}^JV_j$ (sum, not necessarily direct). Define $Q_j:V\to V_j$ ($L^2$ or Euclidean projection), $P_j:V\to V_j$ (elliptic projection, $(AP_ju,v_j)=(Au,v_j)$), and $A_j$ = restriction of $A$ to $V_j$.
- Key identity: $A_jP_j=Q_jA$, i.e. $P_j=A_j^{-1}Q_jA$.
- Additive Schwarz preconditioner and preconditioned operator: $$\mathcal B=\sum_{j=0}^JA_j^{-1}Q_j,\qquad \mathcal BA=\sum_jP_j=:P.$$
- Matrix form with inclusion $I_j:V_j\hookrightarrow V$ ($Q_j=I_j^T$): $\mathcal B=\sum_jI_jA_j^{-1}I_j^T$. Each subproblem is of size $\dim V_j$ (embarrassingly parallel).
- If $V_j=\mathrm{span}\{e_i:i\in I_j\}$ then $A_j$ is the principal submatrix on $I_j$. With $N$ one-dimensional spaces $V_j=\mathrm{span}\{e_j\}$ the preconditioner is Jacobi (diagonal scaling).
- Inexact local solves: replace $A_j^{-1}$ by SPD $\mathcal B_j$ giving $\tilde{\mathcal B}=\sum_j\mathcal B_jQ_j$ and $T=\sum_jT_j$ with $T_j=\mathcal B_jQ_jA$. Each $T_j$ is SPD in the $A$-inner product.
- Lemma (fundamental, Xu / Dryja-Widlund): $$(T^{-1}v,v)_A=\inf_{\sum v_j=v}\sum_{j=0}^J(T_j^{-1}v_j,v_j)_A.$$ For exact solves $T_j=P_j$ this equals $\inf_{\sum v_j=v}\sum_j(v_j,v_j)_A$. If $\mathcal B_j\sim A_j^{-1}$ the same holds up to constants.
- Strategy to bound $\kappa(\mathcal BA)$: bound $\lambda_{\max}(T)$ via the finite overlap (strengthened Cauchy-Schwarz) property, and $\lambda_{\min}(T)$ via a stable decomposition: for every $v$ there exist $v_j\in V_j$ with $\sum v_j=v$ and $\sum_ja(v_j,v_j)\le K\,a(v,v)$. Then $\lambda_{\min}\ge1/K$.
## 29. Two-level method with vertex patches (L14)
- Setup: $V_h$ = P1 FE space on mesh $\mathcal T$ obtained by refining a coarse mesh $\mathcal T_0$. $V_0$ = P1 on the coarse mesh ($V_0\subset V_h$). For each fine vertex $j$, $V_j=\mathrm{span}\{\phi_j\}$ (supported on the vertex patch $\Omega_j$, the union of triangles touching vertex $j$).
- Preconditioner: $$\mathcal B=I_0A_0^{-1}I_0^T+\sum_{j=1}^J I_jA_j^{-1}I_j^T,$$ where the vertex part is simply Jacobi, $\mathcal B_D=\mathrm{diag}(A)^{-1}$. $V_0$ is not needed to span $V_h$ but is essential for effectiveness.
- Upper bound on $\lambda_{\max}$ (finite overlap): $a(\phi_i,\phi_j)\ne0$ only if $i,j$ share a triangle. With maximum vertex valence $D$ a crude argument gives $\lambda_{\max}(T)\le2(D+1)$. Sharper: each triangle lies in exactly 4 subdomains (the coarse domain plus the 3 vertex patches), so $a(Pu,u)\le\sum_ja_j(u,u)\le4a(u,u)$ using $a_j(P_ju,P_ju)\le a_j(u,u)$. Thus $\lambda_{\max}(T)\le4$.
- Lower bound on $\lambda_{\min}$ (stable decomposition): take $v_0$ = coarse interpolant of $v$, which satisfies $a(v_0,v_0)\lesssim a(v,v)$ and $\|v-v_0\|^2_{L^2}\lesssim h^2a(v,v)$. Then $v_j=(v-v_0)(x_j)\phi_j$ and a scaling argument on the reference triangle gives $\|\nabla v_j\|^2_{L^2}\lesssim h^{-2}\|v-v_0\|^2_{L^2(\Omega_j)}+\|\nabla(v-v_0)\|^2_{L^2(\Omega_j)}\lesssim\|\nabla v\|^2_{L^2(\Omega_j)}$. Summing with finite overlap: $\sum_ja(v_j,v_j)\le4a(v,v)$ (independent of $h$).
- Conclusion: $\kappa(\mathcal BA)=O(1)$, independent of $h$. So PCG takes asymptotically constant iterations. Without the coarse space ($\mathcal B_D$ alone, Jacobi) $\kappa=O(h^{-2})$ (no improvement in $h$-scaling).
- Complexity: vertex patch solves are trivial (1x1). The coarse solve $A_0^{-1}$ is the main cost, with $A_0$ about $\frac14$ the size of $A$ in 2D. Not scalable asymptotically. Remedy: apply the construction recursively (multilevel additive Schwarz, BPX, multigrid), which gives $O(1)$ condition and $O(N)$ complexity.
## 30. Homework 4 highlights
- Proposition (L13) proof outline: $M^{-1}A$ is similar to $M^{-1/2}AM^{-1/2}$ (SPD), which is also similar to $A^{1/2}M^{-1}A^{1/2}$ and shares eigenvalues with the generalized problem. Eigenvalues are real and positive. The extremal values are the min/max of the Rayleigh quotient $\frac{x^TAx}{x^TMx}$.
- Two-level code: fine grid has $n_x$ (odd) vertices per direction, coarse grid has $(n_x+1)/2$. $J$ is the inclusion (prolongation) matrix of linear interpolation, with boundary rows/columns zeroed. $B_0=JA_0^{-1}J^T$ (coarse solver) and $B_D=\mathrm{diag}(A)^{-1}$ (vertex patch solver), with $B=B_0+B_D$.
- Expected CG iteration counts for $n_x\in\{11,31,51,71\}$: unpreconditioned CG grows like $O(h^{-1})$ (since $\kappa\sim h^{-2}$). $B_D$ (Jacobi) alone gives essentially the same growth because on a uniform grid $\mathrm{diag}(A)$ is constant. $B=B_0+B_D$ gives a roughly constant iteration count (mesh independent).
- Spectra: $A$ has $\lambda_{\max}=O(1)$ and $\lambda_{\min}=O(h^2)$. $B_DA$ has $\lambda_{\max}\le C$ (bounded by finite overlap) and $\lambda_{\min}=O(h^2)$, so the same poor condition number. $BA$ has $\lambda_{\max}\le4$ (finite overlap) and $\lambda_{\min}\ge c>0$ independent of $h$ (stable decomposition), so $\kappa=O(1)$.
# Quick Reference: Key Formulas and Facts
- Heat equation semi-discrete: $M\dot u+Au=f$; FE stable iff $\Delta t\lesssim h^2$; BE and CN unconditionally stable; CN is 2nd order.
- Cea (non-symmetric): $\|u-u_h\|\le\frac\gamma\alpha\inf\|u-v_h\|$. Lax-Milgram needs continuity plus coercivity only.
- Elasticity: $a(u,v)=\int2\mu\,\varepsilon(u):\varepsilon(v)+\lambda\,\nabla\cdot u\,\nabla\cdot v$. Korn gives coercivity. $\ker\varepsilon=RM$. Locking as $\lambda\to\infty$.
- GMRES: $\min\|\beta e_1-H_my\|$; $AQ_m=Q_{m+1}H_m$; memory $O(mN)$ and work $O(m^2N)$. MINRES uses Lanczos with a three-term recurrence.
- Saddle point: $\mathcal A$ has $n$ positive and $m$ negative eigenvalues; $S=-BA^{-1}B^T$; need coercivity on $\ker B$ and inf-sup.
- Fortin operator implies discrete inf-sup. P1-P0 is unstable; P2-P0 and Taylor-Hood $P_k$-$P_{k-1}$ ($k\ge2$) are stable.
- Ideal block preconditioner gives eigenvalues $1,\frac{1\pm\sqrt5}{2}$. Stokes: $-S\sim M_p$.
- Additive Schwarz: $\mathcal B=\sum I_jA_j^{-1}I_j^T$, $(T^{-1}v,v)_A=\inf_{\sum v_j=v}\sum(v_j,v_j)_A$. Two-level with coarse space gives $\kappa=O(1)$ via finite overlap (upper) and stable decomposition (lower).
