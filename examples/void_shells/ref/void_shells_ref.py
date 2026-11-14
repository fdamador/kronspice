# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent Python reference for examples/void_shells.

Same physics as field.go, written separately (numpy stencils; light modes from
scipy.sparse.linalg.eigsh, shift-invert). TestBubbleAndShell checks the Go
solver against its round-bubble numbers (ether 106.48697, light 285.33622 at
eps = 300, V0 = 20, h = 0.4). Usage:

    python3 examples/void_shells/ref/void_shells_ref.py EPS V0 K STEPS     # one photon in mode K
"""
import numpy as np, sys, json
import scipy.sparse as sp, scipy.sparse.linalg as sla
h=0.4; Lr=24.0; Lz=24.0
nr=int(Lr/h); nz=int(2*Lz/h); N=nr*nz
rho=(np.arange(nr)+0.5)*h; z=(np.arange(nz)+0.5)*h-Lz
P,Z=np.meshgrid(rho,z,indexing='ij'); rp=rho+h/2; rm=rho-h/2; dV=2*np.pi*P*h*h; m=0.5
def lap(f):
    fp=np.vstack([f[1:],f[-1:]]); fm=np.vstack([f[:1],f[:-1]])
    g=(rp[:,None]*(fp-f)-rm[:,None]*(f-fm))/(rho[:,None]*h*h)
    return g+(np.hstack([f[:,1:],f[:,-1:]])-2*f+np.hstack([f[:,:1],f[:,:-1]]))/(h*h)
# light operator
I=np.arange(N).reshape(nr,nz); rows=[];cols=[];vals=[]; dg=np.zeros((nr,nz))
for i in range(nr):
    dg[i,:]+=rp[i]/h**2+2*rho[i]/h**2+rho[i]*m*m/rho[i]**2
    if i>0: dg[i,:]+=rm[i]/h**2
for i in range(nr):
    for j in range(nz):
        k=I[i,j]
        if i+1<nr: rows+=[k];cols+=[I[i+1,j]];vals+=[-rp[i]/h**2]
        if i>0: rows+=[k];cols+=[I[i-1,j]];vals+=[-rm[i]/h**2]
        if j+1<nz: rows+=[k];cols+=[I[i,j+1]];vals+=[-rho[i]/h**2]
        if j>0: rows+=[k];cols+=[I[i,j-1]];vals+=[-rho[i]/h**2]
K0=sp.csr_matrix((vals,(rows,cols)),shape=(N,N))+sp.diags(dg.ravel()); Bm=sp.diags(np.repeat(rho,nz))
def modes(dens,V0,k):
    A=K0+sp.diags((V0*dens*P).ravel())
    lam,vec=sla.eigsh(A,k=k,M=Bm,sigma=0.0,which='LM')
    o=np.argsort(lam); lam=lam[o]; vec=vec[:,o]
    out=[]
    for q in range(k):
        c=vec[:,q].reshape(nr,nz); c/=np.sqrt(np.sum(c**2*dV)); out.append(c)
    return lam,out
def ether_E(psi):
    dens=abs(psi)**2
    gx=(np.vstack([psi[1:],psi[-1:]])-psi)/h; gz=(np.hstack([psi[:,1:],psi[:,-1:]])-psi)/h
    return np.sum((0.5*(abs(gx)**2+abs(gz)**2)+0.25*(dens-1)**2)*dV)
def relax(eps,V0,K,steps=20000,every=25,R0=None):
    # start: a ball void of radius R0 with the K-th mode
    if R0 is None: R0=4+2*K
    r=np.hypot(P,Z)
    psi=np.tanh(np.maximum(r-R0,0)/1.0+1e-3).astype(complex)
    lam,ch=modes(abs(psi)**2,V0,K+2); chi=ch[K-1]; lk=lam[K-1]
    dt=0.15*h*h; hist=[]
    for it in range(steps):
        c=eps*V0/(2*np.sqrt(lk))
        G=-0.5*lap(psi)+0.5*(abs(psi)**2-1)*psi+c*chi**2*psi
        psi=psi-dt*G
        if it%every==every-1:
            lam,ch=modes(abs(psi)**2,V0,K+2)
            # follow the mode by overlap
            q=K-1
            chi=ch[q]*np.sign(np.sum(chi*ch[q]*dV)); lk=lam[q]
        if it%2000==0 or it==steps-1:
            Ee=ether_E(psi); EL=eps*np.sqrt(lk)
            hist.append((it,Ee,EL,(Ee+EL)/EL,q if it>0 else K-1))
    return psi,chi,lk,hist
if __name__=="__main__":
    eps=float(sys.argv[1]); V0=float(sys.argv[2]); K=int(sys.argv[3]); st=int(sys.argv[4])
    psi,chi,lk,hist=relax(eps,V0,K,st)
    for hh in hist: print(hh)
    np.savez(f"/tmp/void_shells_ref_{int(eps)}_{K}.npz",psi=psi,chi=chi,lk=lk)
