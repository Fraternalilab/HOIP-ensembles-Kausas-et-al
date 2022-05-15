# Setting misc options 
hide all 
bg_color white 
unset opaque_background 
unset depth_cue 

# Loading cluster leaders 
load MP1.pdb 
load MP2.pdb 
load MP3.pdb 
load MP4.pdb
load MP5.pdb 
load MP6.pdb 
load MP7.pdb 
load MP8.pdb 

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

# Align each cluster leader on IBR 
select alignment, model basis_set_1 and resi 796-841 
align MP1, alignment 
align MP2, alignment 
align MP3, alignment 
align MP4, alignment 
align MP5, alignment 
align MP6, alignment 
align MP7, alignment 
align MP8, alignment
