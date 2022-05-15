# Setting misc options 
hide all 
bg_color white 
unset opaque_background 
unset depth_cue 

# Loading cluster leaders 
load cluster_leader_0.pdb 
load cluster_leader_15.pdb 
load cluster_leader_12.pdb 
load cluster_leader_17.pdb 
load cluster_leader_7.pdb 
load cluster_leader_26.pdb 
load cluster_leader_11.pdb 
load cluster_leader_5.pdb 
load cluster_leader_24.pdb 
load cluster_leader_19.pdb 
load cluster_leader_3.pdb 
load cluster_leader_20.pdb 
load cluster_leader_10.pdb 
load cluster_leader_30.pdb 
load cluster_leader_14.pdb 
load cluster_leader_18.pdb 
load cluster_leader_4.pdb 
load cluster_leader_23.pdb 
load cluster_leader_29.pdb 
load cluster_leader_25.pdb 
load cluster_leader_8.pdb 
load cluster_leader_1.pdb 
load cluster_leader_22.pdb 
load cluster_leader_16.pdb 
load cluster_leader_31.pdb 
load cluster_leader_21.pdb 
load cluster_leader_28.pdb 
load cluster_leader_27.pdb 
load cluster_leader_6.pdb 
load cluster_leader_32.pdb 
load cluster_leader_2.pdb 
load cluster_leader_13.pdb 
load cluster_leader_9.pdb 

# More misc options 
as cartoon 
set ray_trace_gain, 0 
set ray_trace_mode, 1 

# Domain selection 
select RING1, resi 699-751 
select L1, resi 752-795 
select IBR, resi 796-841 
select L2, resi 842-867 
select RING2, resi 868-1071 
select ZINCS, name ZN 
select C885, resi 885 

# Domain colors 
color skyblue, RING1 
color gray80, L1 
color tv_orange, IBR 
color gray80, L2 
color forest, RING2 
color tv_yellow, C885 
show spheres, C885  
show spheres, ZINCS 

# Cluster colors 
# Align each cluster leader on IBR 
select alignment, model cluster_leader_0 and resi 796-841 
align cluster_leader_0, alignment 
align cluster_leader_15, alignment 
align cluster_leader_12, alignment 
align cluster_leader_17, alignment 
align cluster_leader_7, alignment 
align cluster_leader_26, alignment 
align cluster_leader_11, alignment 
align cluster_leader_5, alignment 
align cluster_leader_24, alignment 
align cluster_leader_19, alignment 
align cluster_leader_3, alignment 
align cluster_leader_20, alignment 
align cluster_leader_10, alignment 
align cluster_leader_30, alignment 
align cluster_leader_14, alignment 
align cluster_leader_18, alignment 
align cluster_leader_4, alignment 
align cluster_leader_23, alignment 
align cluster_leader_29, alignment 
align cluster_leader_25, alignment 
align cluster_leader_8, alignment 
align cluster_leader_1, alignment 
align cluster_leader_22, alignment 
align cluster_leader_16, alignment 
align cluster_leader_31, alignment 
align cluster_leader_21, alignment 
align cluster_leader_28, alignment 
align cluster_leader_27, alignment 
align cluster_leader_6, alignment 
align cluster_leader_32, alignment 
align cluster_leader_2, alignment 
align cluster_leader_13, alignment 
align cluster_leader_9, alignment 

