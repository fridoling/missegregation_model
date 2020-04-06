#!/usr/bin/env python
# coding: utf-8

import numpy as np
import copy
import os


class Cell:
    def __init__(self, n_chroms=16, ploidy=1, strain='None', missegregation=0.001, ap_loss=0.0, ms8=0.0, ap8_gain=0.0, base_fertility=1.0):
        self.strain = strain
        self.dead = False
        self.aneuploid = False
        self.karyotype = np.ones(n_chroms, dtype=int)*ploidy
        self.base_fertility = base_fertility
        self.fertility = base_fertility
        self.missegregation = missegregation
        self.base_missegregation = missegregation 
        self.ap_loss = ap_loss
        self.ap8_gain = ap8_gain
        self.ms8 = ms8
    def set_karyotype(self, karyotype):
        self.karyotype = karyotype
        if any(karyotype!=1):
            self.aneuploid = True
        if karyotype[8]>1:
            self.missegregation = self.ms8
        else:
            self.missegregation = self.base_missegregation
        if any(karyotype==0):
            self.dead = True
            self.fertility = 0
        else:
            self.update_fertility(karyotype)
    def update_fertility(self, karyotype):
        n_additional_chroms = np.sum(karyotype) - len(karyotype)
        if self.karyotype[8]>1:
            ap8_fac = 1+self.ap8_gain
        else:
            ap8_fac = 1-self.ap_loss            
        total_loss = (1-self.ap_loss)**(n_additional_chroms-1)*ap8_fac
        self.fertility = self.base_fertility*total_loss
        

def divide_cell(Cell, **kwargs):
    daughter1 = copy.deepcopy(Cell)
    daughter2 = copy.deepcopy(Cell)
    kary_mother = Cell.karyotype
    rand_vec = np.random.rand(np.sum(kary_mother))
    kary_daughter1 = kary_mother.copy()
    kary_daughter2 = kary_mother.copy()
    missegregated = (rand_vec<Cell.missegregation).astype(int)
    missegregated*= np.random.choice([-1, 1], len(missegregated))
    k = 0
    for i in range(len(kary_mother)):
        for j in range(kary_mother[i]):
            kary_daughter1[i]+=missegregated[k]
            kary_daughter2[i]-=missegregated[k]            
            k+=1
    daughter1.set_karyotype(kary_daughter1, **kwargs)
    daughter2.set_karyotype(kary_daughter2, **kwargs)
    return daughter1, daughter2

class Population:
    def __init__(self, size=100, **kwargs):
        self.Cells = []
        self.size = size
        for i in range(size):
            self.Cells.append(Cell(**kwargs))
    def add_cell(self, Cell):
        self.Cells.append(Cell)
        self.size+=1
    def remove_cell(self):
        if self.size>1:
            cell_ind = np.random.choice(self.size)
            self.Cells.pop(cell_ind)
            self.size-=1            
    def get_aneuploidy(self, chrom):
        if len(self.Cells)==0:
            return(0)
        chrNumber = np.array([Cell.karyotype[chrom] for Cell in self.Cells])
        return np.bincount(chrNumber)/self.size
    def get_dead(self):
        deadCells = 0
        if len(self.Cells)==0:
            return(0)
        for Cell in self.Cells:
            if Cell.dead:
                deadCells += 1 
        return(deadCells/self.size)

def propagate_population(Pop,  max_size=10000, **kwargs):
    Pop_out = Population(size=0, **kwargs)
    for Mother in Pop.Cells:
        if Mother.dead:
            continue
        if Mother.fertility>np.random.random():
            Daughter1, Daughter2 = divide_cell(Mother, **kwargs)
            Pop_out.add_cell(Daughter1)
            Pop_out.add_cell(Daughter2)
    while Pop_out.size>max_size:
        Pop_out.remove_cell()
    return(Pop_out)

def run_simulation(i, j, params):
    m = params['m_vec'][i]
    f4 = params['frac4_vec'][j]
    pop_size = params['pop_size']
    n_gen = params['n_gen']
    fertility = params['fertility']
    pop4_size = np.int64(f4*pop_size)
    pop1_size = pop_size-pop4_size       
    pop = Population(size=0)
    traj = np.zeros((n_gen, 2))        
    for n in range(pop1_size):
        pop.add_cell(Cell(ap_loss = params['ap_loss'], 
                          strain = 'cdc20x1',
                          missegregation=m,
                          base_fertility=fertility*params['fert_factor'],
                          ap8_gain = params['ap8_gain'],
                          ms8 = params['ms8']
                        ))
    for n in range(pop4_size):
        pop.add_cell(Cell(ap_loss = params['ap_loss'],
                          strain = 'cdc20x4',
                          missegregation=m*params['ms_factor'],
                          base_fertility=fertility,
                          ap8_gain = params['ap8_gain'],
                          ms8 = params['ms8']
                        ))
    traj[0,] = [pop1_size, pop4_size]
    for n in np.arange(1,n_gen):
        pop = propagate_population(pop, max_size=pop_size)
        c1_size = len([cell for cell in pop.Cells if cell.strain=='cdc20x1'])
        c4_size = pop.size - c1_size
        traj[n,0] = c1_size
        traj[n,1] = c4_size
    return(traj)