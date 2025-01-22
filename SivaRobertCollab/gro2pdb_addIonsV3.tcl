# Robert Coffman Jan 29, 2024. Divalent Cation insertion into Integrin structures and preparation for use in CHARMM-GUI.
# Updated Feb 9, 2024 to use the Oreint package to align the protein with the xy and z axes.
# Updated April 12, 2024 to use a more complete structure of integrin Alpha2b Beta3 so that the ion at the genu
# is also added.

# Either copy the 8t2v.pdb structure file into the folder that has the .gro file you want to process
# or set a hard path to where the 8t2v.pdb file is in the line below "mol new 8t2v.pdb". 
# e.g. "mol new /mnt/c/user/downloads/8t2v.pdb waitfor all" without the quotes

# delete all loaded molecules in VMD. load required packages and load the structure files you will work with
mol delete all
package require topotools
package require Orient
namespace import Orient::orient
# structure with the ions to be added to head (alpha chain beta-propeller domain and beta chain beta-I domain)
mol new 8t2v.pdb waitfor all
mol new lframe_box.gro waitfor all
set fid [open "gro2pdb_addIonsV3.log" w]
puts $fid "Details of atoms that were within 0.85 angstrom of any of the added ions. Note: I do not know if this \
will cause problems or if minimization can take care of it."
# rename the chains and segnames for later use
set chaina [atomselect top "resid 1 to 1008"]
set chainb [atomselect top "resid 1009 to 1770"]
$chaina set chain A; $chainb set chain B; $chaina set segname PROA; $chainb set segname PROB
# renumber the residues in chain B so I can line them up and add ions.
for {set x 1009} {$x <= 1770} {incr x} {
	set redo [atomselect top "resid $x"]
	$redo set resid [expr {$x - 1008}]
	$redo delete
}
# reorient the protein so that the z-dimension is aligned with the longest length of the protein
set all [atomselect top all]
set I [draw principalaxes $all]
set A [orient $all [lindex $I 2] {0 0 1}]
$all move $A
set I [draw principalaxes $all]
set A [orient $all [lindex $I 1] {0 1 0}]
$all move $A
$chaina delete; $chainb delete; $all delete;
# get a list of loaded molecules. This helps with not needing to restart VMD every time this is run.
set loadedmols [molinfo list]
# nested list of each ions binding site and reference to the ion itself
set residlist { { {(resid 434 or resid 432 or resid 426 or resid 428 or resid 430) and chain A} {resid 1104 and ion and chain A} } \
{ {(resid 371 or resid 365 or resid 367 or resid 369 or resid 373) and chain A} {resid 1103 and ion and chain A} } \
{ {(resid 299 or resid 297 or resid 303 or resid 305 or resid 301) and chain A} {resid 1101 and ion and chain A} } \
{ {(resid 245 or resid 250 or resid 247 or resid 252 or resid 243) and chain A} {resid 1102 and ion and chain A} } \
{ {(resid 219 or resid 217 or resid 158 or resid 220 or resid 215) and chain B} {resid 805 and ion and chain B} } \
{ {(resid 127 or resid 126 or resid 123) and chain B} {resid 806 and ion and chain B} } \
{ {(resid 607 or resid 642 or resid 605 or resid 602) and chain A} {resid 1109 and ion and chain A} } \
{ {(resid 220 or resid 121) and chain B} {resid 804 and ion and chain B} } } 
# loop through each ion binding sites one by one creating new structures each time with ions in the right place. 
# this could potentially be refined by only aligning the oxygen atoms of the residues that are in contact with the ion.
for {set x 0} {$x < [llength $residlist]} {incr x} {
	# align the binding site
	set closed_ca1 [atomselect [lindex $loadedmols 0] "[lindex $residlist $x 0] and not hydrogen"]
	set string1 [atomselect top "[lindex $residlist $x 0] and not hydrogen"]
	set M [measure fit $closed_ca1 $string1]
	set fitclosed [atomselect [lindex $loadedmols 0] all]	
	$fitclosed move $M
	# combine the divalent cation into the structure produced by Siva
	set sellist {}
	set newstring [atomselect top all]
	set divalent [atomselect [lindex $loadedmols 0] "[lindex $residlist $x 1]"]
	$divalent set chain C
	lappend sellist $newstring
	lappend sellist $divalent
	set mol [::TopoTools::selections2mol $sellist]
	animate write pdb combined${x}.pdb $mol
	mol new combined${x}.pdb waitfor all
	# cleanup
	$closed_ca1 delete; $string1 delete; $fitclosed delete; $newstring delete; $divalent delete
	unset mol; unset sellist; unset M
	#Delete the written file after structure is loaded into memory
	file delete combined${x}.pdb
}
# write final structure
animate write pdb lframe_box.pdb
mol delete all
mol new lframe_box.pdb waitfor all
set overlaptest [atomselect top "protein and within 0.85 of ions"]
if { [$overlaptest num] != 0 } {
	puts "Warning there appears to be protein atoms that are very close to the \
	newly added ions. See gro2pdb_addIonsV3.log for details"
	puts $fid "Resid [$overlaptest get resid] Resname [$overlaptest get resname] Chain [$overlaptest get chain] \
	Segname [$overlaptest get segname] SerialNumber(plumed index number) [$overlaptest get serial]"
	set howclose [atomselect top "protein and within 0.5 of ions"]
	if { [$howclose num] != 0 } {
		puts $fid "These atoms are within 0.5 angstroms of an ion"
		puts $fid "Resid [$howclose get resid] Resname [$howclose get resname] Chain [$howclose get chain] \
		Segname [$howclose get segname] SerialNumber(plumed index number) [$howclose get serial]"
	}
	$howclose delete
}
close $fid; $overlaptest delete
