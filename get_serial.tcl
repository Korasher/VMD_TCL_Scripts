# Robert Coffman. July 12, 2023. Script to print out the serial ids of atoms I want to use for the CV's I currently have designed. 
# they should all be the same between structures. I just want to make sure they actually are. 




set allbinding [atomselect top "residue 462 or residue 469 or residue 470 or residue 501 or residue 558 or residue 560 or residue 562 to 564 or residue 566 or residue 594 or residue 678"]
set heavybinding [atomselect top "(residue 462 or residue 469 or residue 470 or residue 501 or residue 558 or residue 560 or residue 562 to 564 or residue 566 or residue 594 or residue 678) and not hydrogen"]

set ligand [atomselect top "residue 695 to 701"]

set allbindserial [$allbinding get serial]
set heavybindserial [$heavybinding get serial]
set ligandserial [$ligand get serial]

set allbindserial [join $allbindserial ,]
set heavybindserial [join $heavybindserial ,]
set ligandserial [join $ligandserial ,]

puts "all_binding $allbindserial"
puts "heavy_binding $heavybindserial "
puts "ligand $ligandserial"



