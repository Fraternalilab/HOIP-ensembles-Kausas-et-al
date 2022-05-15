# Setting misc options 
bg_color white  
unset opaque_background 
unset depth_cue 
set antialias, 2 
set antialias_shader, 2 
set hash_max, 300 

# Loading cluster leaders 
load cluster_leader_3.pdb 
load cluster_leader_1.pdb 
load cluster_leader_2.pdb 
load cluster_leader_21.pdb 
load cluster_leader_86.pdb 
load cluster_leader_163.pdb 

# More misc options 
set ray_trace_gain, 0 
set ray_trace_mode, 1 
set ray_trace_color, black
set transparency, 0.6 

# Domain selection 
select RING1, resi 699-751 
select L1, resi 752-795 
select IBR, resi 796-841 
select L2, resi 842-868 
select RING2, resi 869-1071 
select C885, resi 885 
select ZINCS, name ZN  

# Domain colors 
color marine, RING1 
color gray80, L1 
color tv_orange, IBR 
color gray80, L2 
color forest, RING2 
color tv_yellow, C885 

as cartoon 
set cartoon_side_chain_helper, 1 
show spheres, C885 
show spheres, ZINCS 
# Cluster colors 
# Align each cluster leader on IBR 
select alignment, model cluster_leader_3 and resi 796-841 
align cluster_leader_1, alignment 
align cluster_leader_2, alignment 
align cluster_leader_21, alignment 
align cluster_leader_86, alignment 
align cluster_leader_163, alignment 
