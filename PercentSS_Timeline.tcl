mol new ref.pdb waitfor all
mol addfile trajout.xtc waitfor all

set outfile [open ./percent_helix.dat w]
set lookup {H G I}
set frame_num [molinfo top get numframes]
set full [atomselect top "chain B and (resid 127 to 146) and name CA"]
set len [llength [$full get resid]]
$full delete

for {set i 0} {$i < $frame_num} {incr i} {
    animate goto $i
    set sel [atomselect top "chain B and (resid 127 to 146) and name CA"]
    mol ssrecalc top
    set struc_string [$sel get structure]
    set helix 0
    foreach letter $lookup {
        set temp [expr {[llength [split $struc_string $letter]] - 1}]
        incr helix $temp
    }
    set percent [expr {double($helix) / double($len) * 100}]
    puts $outfile "$i\t$percent"
    $sel delete
}
close $outfile