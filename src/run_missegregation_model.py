#!/usr/bin/env python
# coding: utf-8

n_gen = 100
pop_size = 100
m_vec = np.arange(0.01,0.1,0.01)
frac4_vec = np.arange(0.1,1.0,0.1)
ms_factor = 2
fertility = 1
fert_factor = 0.9

trajs = {}
trajs['cdc20x1'] = np.zeros((len(m_vec), len(frac4_vec), n_gen))
trajs['cdc20x4'] = np.zeros((len(m_vec), len(frac4_vec), n_gen))
res = np.zeros((len(m_vec), len(frac4_vec)))

for i,m in zip(range(len(m_vec)), m_vec):
    for j,f4 in zip(range(len(frac4_vec)), frac4_vec):
        pop4_size = np.int64(f4*pop_size)
        pop1_size = pop_size-pop4_size       
        pop_in = Population(size=0)
        pops = [pop_in]
        traj_c1 = np.zeros(n_gen)        
        traj_c4 = np.zeros(n_gen)
        for n in range(pop1_size):
            pop_in.add_cell(Cell(ap_loss = 0.1, 
                              missegregation=m,
                              base_fertility=fertility*fert_factor
                             ))
        for n in range(pop4_size):
            pop_in.add_cell(Cell(ap_loss = 0.1,
                              missegregation=m*ms_factor,
                              base_fertility=fertility))   
        for n in range(n_gen):
            pop_out = propagate_population(pop_in, max_size=pop_size)
            pops.append(pop_out)
            c1_size = len([cell for cell in pop_out.Cells if cell.base_missegregation==m])
            c4_size = pop_out.size - c1_size
            pop_in = pop_out
            traj_c1[n] = c1_size
            traj_c4[n] = c4_size
            if c1_size==0:
                res[i,j] = 1
                break
            elif c4_size==0:
                break
            elif n==(n_gen-1):
                res[i,j] = 0.5
        trajs['cdc20x1'][i,j,] = traj_c1
        trajs['cdc20x4'][i,j,] = traj_c4      