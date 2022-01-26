#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 19:47:16 2021

@author: Gabriel Martinez

Program used to calculate badn structure of CdTe under isotropic strain
"""

import ase 
from ase.build import bulk
from ase.visualize import view
from ase.io import write, Trajectory
from gpaw import GPAW, PW, FermiDirac, mpi
from ase.optimize import BFGS

import numpy as np
import matplotlib.pyplot as plt

from ase.dft.kpoints import bandpath
from gpaw.spinorbit import soc_eigenstates
#============================================================================
def optimizacion(estr): # structure optimization
    calc=GPAW( mode=PW(650),
         xc='PBE',
         kpts=(14,14,14),
         occupations=FermiDirac(0.01),
        txt='CdTe_PW_opt.txt')
    estr.calc=calc
    opt = BFGS(estr)
    traj = Trajectory('path_opt.traj', 'w', estr)
    opt.attach(traj)
    opt.run(0.005)

def autoCons(estr,st): #self-consistent calculus DFT
    calc=GPAW( mode=PW(650),
         xc='PBE',
         kpts=(14,14,14),
         occupations=FermiDirac(0.01),
        txt='CdTe_PW_opt_estr_{0:.2f}.txt'.format(st))
    estr.calc=calc
    estr.get_potential_energy()
    calc.write('CdTe_Pw_estr_{0:.2f}.gpw'.format(st))

def calcBandas(est,st,nomgpw): # band calculation CdTe
    calc=GPAW(nomgpw).fixed_density(
        nbands=40,
        symmetry='off',
        kpts={'path': 'WLGXWK', 'npoints': 200},
        convergence={'bands': 8},
        txt= 'CdTe_PW_opt_estr_{0:.2f}_bandas.txt'.format(st))
    calc.write('CdTe_bandas_estr_{0:.2f}.gpw'.format(st))
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lat=6.6; b=1.2
ang=90
fin=25
inicio=-5
intr=0.01
listaEst=[]
strain =np.arange(inicio,fin,1)
print("begining self-consistent calculation")
for i in strain:
    st=intr*i; latD=(1+st)*lat
    print("def: {0:.2f}".format(st))
    CdTe=bulk('CdTe', 'zincblende', a=latD)
    autoCons(CdTe,st)
    calcBandas(CdTe,st, 'CdTe_Pw_estr_{0:.2f}.gpw'.format(st))
    write("CdTe_{0:.2f}.xyz".format(st),CdTe)
    listaEst.append(CdTe)