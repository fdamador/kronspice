# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent Python reference for Part 13 (the full ether simulation).

Same physics as ethersim.go, written separately: numpy stencils, the light mode
from scipy.sparse.linalg.eigsh (shift-invert) instead of Go's CG inverse
iteration, RK4 in real time. Used to check the Go solver (TestEtherSimReference,
TestEtherSimRealTime use its numbers). Run:

    H=0.4 python3 examples/ether_horizon/ref/ether13_ref.py EPS V0 R0 A0 U T DT EVERY TAG

It writes /tmp/full/rt2_TAG.npz and .json; edit the paths to keep them.
"""
import numpy as np, sys, json, time
import scipy.sparse as sp, scipy.sparse.linalg as sla
import os
h=float(os.environ.get('H','0.4')); Lr=20.0; Lz=26.0
nr=int(Lr/h); nz=int(2*Lz/h); N=nr*nz
rho=(np.arange(nr)+0.5)*h; z=(np.arange(nz)+0.5)*h-Lz
P,Z=np.meshgrid(rho,z,indexing='ij'); rp=rho+h/2; rm=rho-h/2; dV=2*np.pi*P*h*h; m=0.5
def lap(f):
    fp=np.vstack([f[1:],f[-1:]]); fm=np.vstack([f[:1],f[:-1]])
    g=(rp[:,None]*(fp-f)-rm[:,None]*(f-fm))/(rho[:,None]*h*h)
    return g+(np.hstack([f[:,1:],f[:,-1:]])-2*f+np.hstack([f[:,:1],f[:,:-1]]))/(h*h)
def dz(f): return (np.hstack([f[:,1:],f[:,-1:]])-np.hstack([f[:,:1],f[:,:-1]]))/(2*h)
# light operator: rho-weighted symmetric stiffness, Dirichlet outside
idx=np.arange(N).reshape(nr,nz)
rows=[];cols=[];vals=[]
def add(a,b,v): rows.append(a);cols.append(b);vals.append(v)
diag0=np.zeros((nr,nz))
for i in range(nr):
    for j in range(nz):
        k=idx[i,j]; d=0.0
        # rho direction
        d+=rp[i]/h**2
        if i+1<nr: add(k,idx[i+1,j],-rp[i]/h**2)
        if i>0: d+=rm[i]/h**2; add(k,idx[i-1,j],-rm[i]/h**2)
        d+=2*rho[i]/h**2
        if j+1<nz: add(k,idx[i,j+1],-rho[i]/h**2)
        if j>0: add(k,idx[i,j-1],-rho[i]/h**2)
        d+=rho[i]*m*m/rho[i]**2
        diag0[i,j]=d
K0=sp.csr_matrix((vals,(rows,cols)),shape=(N,N))+sp.diags(diag0.ravel())
Bm=sp.diags(np.repeat(rho,nz))
def light(dens,V0,x0=None,sig=None):
    A=K0+sp.diags((V0*dens*P).ravel())
    lam,vec=sla.eigsh(A,k=1,M=Bm,sigma=sig if sig else 0.0,which='LM',v0=x0)
    chi=vec[:,0].reshape(nr,nz); chi/=np.sqrt(np.sum(chi**2*dV)); 
    if chi[np.unravel_index(np.argmax(abs(chi)),chi.shape)]<0: chi=-chi
    return lam[0],chi
def ring(R0,a0):
    d1=np.hypot(P-R0,Z); d2=np.hypot(P+R0,Z)
    ph=((P-R0)+1j*Z)/np.maximum(d1,1e-9)*((P+R0)-1j*Z)/np.maximum(d2,1e-9)
    return ph*np.tanh(np.maximum(d1-a0,0)/1.0+1e-3)
def wind(psi,b,c,zc=0):
    i1=min(int(b/h),nr-1); j0,j1=max(int((zc-c+Lz)/h),1),min(int((zc+c+Lz)/h),nz-2)
    path=[psi[0,j] for j in range(j1,j0,-1)]+[psi[i,j0] for i in range(0,i1)]+[psi[i1,j] for j in range(j0,j1)]+[psi[i,j1] for i in range(i1,0,-1)]
    ph=np.angle(np.array(path)); d=np.diff(np.concatenate([ph,ph[:1]])); d=(d+np.pi)%(2*np.pi)-np.pi
    return round(d.sum()/(2*np.pi),2)
def energies(psi,eps,lam,U):
    dens=abs(psi)**2
    gx=(np.vstack([psi[1:],psi[-1:]])-psi)/h; gz=(np.hstack([psi[:,1:],psi[:,-1:]])-psi)/h
    Ee=np.sum((0.5*(abs(gx)**2+abs(gz)**2)+0.25*(dens-1)**2)*dV); Pz=np.sum(np.imag(np.conj(psi)*dz(psi))*dV)
    return Ee, eps*np.sqrt(lam), Pz
def shape(psi,chi):
    dens=abs(psi)**2; void=dens<0.5
    if not void.any(): return dict(void=False)
    w=np.where(void,dV,0); zc=np.sum(Z*w)/w.sum()
    return dict(void=True,zc=zc,inner=P[void].min(),outer=P[void].max(),height=Z[void].max()-Z[void].min()+h,
                axis=bool(void[0].any()),wind=wind(psi,min(P[void].max()+4,Lr-1),(Z[void].max()-Z[void].min())/2+5,zc),
                light_in=float(np.sum(np.where(void,chi**2*dV,0))),vol=float(w.sum()))
def run(eps,V0,R0,a0,U,T,dt,every,tag):
    psi=ring(R0,a0).astype(complex)
    lam,chi=light(abs(psi)**2,V0)
    E0=None; out=[]; frames=[]; t=0.0; n=int(round(T/dt))
    def rhs(ps,ch,lm):
        c=eps*V0/(2*np.sqrt(lm))
        return -1j*(-0.5*lap(ps)+0.5*(abs(ps)**2-1)*ps+c*ch**2*ps) - U*dz(ps)
    for k in range(n+1):
        if k%every==0:
            Ee,EL,Pz=energies(psi,eps,lam,U); E=Ee+EL-U*Pz
            if E0 is None: E0=E
            sh=shape(psi,chi); sh.update(t=t,Ee=Ee,EL=EL,Pz=Pz,drift=(E-E0)/abs(E0),g0=(Ee+EL)/EL,lam=lam)
            out.append(sh); frames.append((t,abs(psi)**2,chi**2))
            print(json.dumps({a:(round(float(b),4) if not isinstance(b,bool) else b) for a,b in sh.items()}),flush=True)
            if not sh['void'] or not np.isfinite(Ee): break
        # RK4 with the light frozen over the step, re-solved after it (Born–Oppenheimer)
        k1=rhs(psi,chi,lam); k2=rhs(psi+0.5*dt*k1,chi,lam); k3=rhs(psi+0.5*dt*k2,chi,lam); k4=rhs(psi+dt*k3,chi,lam)
        psi=psi+dt/6*(k1+2*k2+2*k3+k4)
        lam,chi=light(abs(psi)**2,V0,chi.ravel(),lam*0.95); t+=dt
    np.savez_compressed(f"/tmp/full/rt2_{tag}.npz",t=[f[0] for f in frames],dens=np.array([f[1] for f in frames],dtype=np.float32),chi2=np.array([f[2] for f in frames],dtype=np.float32))
    json.dump(out,open(f"/tmp/full/rt2_{tag}.json","w"))
if __name__=="__main__":
    eps,V0,R0,a0,U,T,dt=map(float,sys.argv[1:8]); run(eps,V0,R0,a0,U,T,dt,int(sys.argv[8]),sys.argv[9])
