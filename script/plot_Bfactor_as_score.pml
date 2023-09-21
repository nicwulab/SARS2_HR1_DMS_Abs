set cartoon_transparency, 0
set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
bg_color white

load PDB/6vxx_score.pdb
create S2, 6vxx_score and chain A and resi 883-1034
spectrum b, red white blue, minimum=10, maximum=90
color white, chain A and (not resi 883-1034)
color white, chain B
color white, chain C
set cartoon_transparency, 0.6
set cartoon_transparency, 0, chain A and resi 883-1034

set_view (\
    -0.957329690,   -0.010144172,    0.288813770,\
     0.287655443,   -0.129408702,    0.948948264,\
     0.027748875,    0.991538882,    0.126804829,\
     0.000258466,   -0.000234723, -299.691131592,\
   212.878936768,  218.726196289,  195.179901123,\
   256.356872559,  342.894775391,  -20.000000000 )
hide all
show cartoon, S2
png graph/structure_adj_bind_score.png, dpi=300
