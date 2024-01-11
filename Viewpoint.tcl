
mol new step4.1_equilibration.gro
display resetview
animate style Once
set top [molinfo top]
menu graphics on
mol modselect 0 $top residue 0 to 694
mol modstyle 0 $top NewCartoon 0.300000 10.000000 4.100000 0
#In any publication of scientific results based in part or
#completely on the use of the program STRIDE, please reference:
# Frishman,D & Argos,P. (1995) Knowledge-based secondary structure
# assignment. Proteins: structure, function and genetics, 23, 566-579.

mol material Opaque
mol addrep $top
mol modselect 1 $top residue 695 to 701
mol modstyle 1 $top VDW 1.000000 12.000000
mol color Name
mol representation VDW 1.000000 12.000000
mol selection residue 695 to 701
mol material Opaque
mol addrep $top
mol modselect 2 $top index 3379 or index 6987 or index 8491 or index 8541 or index 2836
mol modcolor 2 $top ColorID 7
mol color ColorID 7
mol representation VDW 1.000000 12.000000
mol selection index 3379 or index 6987 or index 8491 or index 8541 or index 2836
mol material Opaque
mol addrep $top
mol modselect 3 $top index 2849
mol modcolor 3 $top ColorID 4
mol color ColorID 4

source "C:/Program Files/VMD/FMAP_v1-master/funnel_gui/funnel.tcl"
funnel_tk