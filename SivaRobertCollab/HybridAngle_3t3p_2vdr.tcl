
# find the chain names, or if there aren't any, the segnames that are in the current molecule
proc getchain_segnames {molid} {
	global errorlog
	set test [atomselect $molid "all"]
	set chains [lsort -unique [$test get chain]]
	set segnames [lsort -unique [$test get segname]]
	#puts $chains
	#puts "$segnames [llength $segnames]"
	if {[string trim $segnames] != "{}"} {
		return [list segnames $segnames]
	} elseif {[string trim $chains] != "{}" } {
		return [list chains $chains]
	} else {
		puts $errorlog "ERROR could not find any segnames or chains"
		return
	}
}

# rename chain/segname 1 to chain A and segname PROA and rename chain/segname 2 to chain B and segname PROB.
proc setchainnames {chains y} {
	global RMSDlog
	if { [lindex $chains 0] eq "segnames"} {
		if {[llength [lindex $chains 1]] > 1 } {
			puts $RMSDlog "found more than 1 segname"
			set test [atomselect $y "segname [lindex [lindex $chains 1] 0]"]
			$test set chain A
			$test set segname PROA
			set test2 [atomselect $y "segname [lindex [lindex $chains 1] 1]"]
			$test2 set chain B
			$test2 set segname PROB
			$test delete
			$test2 delete
		} elseif { [llength [lindex $chains 1]] == 1 } {
			puts $RMSDlog "found only 1 segname"
			set test [atomselect $y "segname [lindex $chains 1]"]
			$test set chain A
			$test set segname PROA
			$test delete
		}
	} elseif { [lindex $chains 0] eq "chains"} {
		if {[llength [lindex $chains 1]] > 1 } {
			puts $RMSDlog "found more than 1 chain"
			set test3 [atomselect $y "chain [lindex [lindex $chains 1] 0]"]
			$test3 set segname PROA
			$test3 set chain A
			set test4 [atomselect $y "chain [lindex [lindex $chains 1] 1]"]
			$test4 set segname PROB
			$test4 set chain B
			$test3 delete
			$test4 delete
		} elseif { [llength [lindex $chains 1]] == 1} {
			puts $RMSDlog "found only 1 chain"
			set test3 [atomselect $y "chain [lindex $chains 1]"]
			$test3 set segname PROA
			$test3 set chain A
			$test3 delete
		}
	}
}

# Function to calculate the vector between two points
proc vector_subtract {pointA pointB} {
    return [list [expr {[lindex $pointA 0] - [lindex $pointB 0]}] \
                 [expr {[lindex $pointA 1] - [lindex $pointB 1]}] \
                 [expr {[lindex $pointA 2] - [lindex $pointB 2]}]]
}

# Function to calculate the dot product of two vectors
proc dot_product {vectorA vectorB} {
    return [expr {[lindex $vectorA 0] * [lindex $vectorB 0]} + \
                 {[lindex $vectorA 1] * [lindex $vectorB 1]} + \
                 {[lindex $vectorA 2] * [lindex $vectorB 2]}]
}

# Function to calculate the magnitude of a vector
proc magnitude {vector} {
    return [expr {sqrt([dot_product $vector $vector])}]
}

# Which molecules do you want to be references? 0 = first molecule loaded, 1 = first and second molecule loaded, etc.
set number 0
# adjust these as necessary
# where is the iteration file?
set root "D:/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/final_string/analyze_orient"
# assuming the subfolders of iter* are named by increasing number
set subfolders "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
# where are the pdb files that you prepared for this script to use? This is also the path where the RMSD output will be written to and
# where directories will be made for the aligned structures to be written to.
set otherroot "D:/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/PDBstructures_2"

cd $root
# beginning of work code
mol delete all

# load the state1 and other closed structure
cd "${otherroot}"
mol new 3t3p.psf waitfor all
mol addfile 3t3p.pdb waitfor all
mol rename top closed_2vdr
mol new 2vdr.psf waitfor all
mol addfile 2vdr.pdb waitfor all
mol rename top open_2vdr
# # load the first string of pdbs from Siva
# foreach x $subfolders {
# 	cd "${root}/target_md_${x}"
# 	mol new minimize-protein.pdb waitfor all
# 	mol rename top "target_md_${x}"
# }
# load the final_string pdbs from Siva
# foreach x $subfolders {
# 	cd "${root}/target_md_${x}"
# 	mol new minimize-notwater-nocofactors_CA_aligned0Plane.pdb waitfor all
# 	mol rename top "target_md_${x}"
# }
# get a list of all the loaded molecules
set mols [molinfo list]

# move to the correct directory and open files for writing
cd $root
set fileid [open HybridAngles_3t3p_2vdr.dat w]
set RMSDlog [open HybridAngle.log w+]
set errorlog [open HybridAngle.err w+]
puts $RMSDlog "working directory is [pwd]"

# rename the first and second chains or segnames to A and B or PROA and PROB respectively
foreach y $mols {
	set chains [getchain_segnames $y]
	setchainnames $chains $y
}

# start the actual comparisons (align and measure angle)
for {set zim 0 } {$zim <= $number} {incr zim} {
	# Make atom selection for the appropriate domains
	set hy1 [atomselect [lindex $mols $zim] "segname PROB and (resid 58 to 109 or resid 353 to 432) and name CA and not hydrogen"]
	set bi [atomselect [lindex $mols $zim] "segname PROB and resid 110 to 353 and name CA and not hydrogen"]

	# measure the angle of the hybrid domain in all strings in relation to the position of the hybrid domain in the closed state. state1 (3zdy_2.pdb)
	for { set ye 1} {$ye < [llength $mols] } {incr ye} {
		# Define the hybrid and beta-I (Beta-A) domain in the string molecule
		set hybrid [atomselect [lindex $mols $ye] "segname PROB and (resid 58 to 109 or resid 353 to 432) and name CA and not hydrogen"]
		set stringbi [atomselect [lindex $mols $ye] "segname PROB and resid 110 to 353 and name CA and not hydrogen"]
		# align the beta-I domains
		set fit [measure fit $stringbi $bi];
		set all [atomselect [lindex $mols $ye] all]
		$all move $fit
		
		# updating all the selections after the alignment for paranoia's sake.
		$hy1 update; $bi update; $hybrid update; $stringbi update
		
		# Define the coordinates of the three points
		set point1 [measure center $hy1 weight mass];
		set point2 [measure center $bi weight mass];
		set point3 [measure center $hybrid weight mass];
		$hybrid delete; $stringbi delete; $all delete
		
		# Calculate vectors
		set vectorA [vector_subtract $point1 $point2]
		set vectorB [vector_subtract $point3 $point2]
	
		# Calculate the angle (in radians)
		set cos_theta [expr {[dot_product $vectorA $vectorB] / ([magnitude $vectorA] * [magnitude $vectorB])}]
		set angle [expr {acos($cos_theta)}]
		
		# Convert the angle to degrees
		set angle_degrees [expr {$angle * 180 / 3.141592653589793}]
		
		# Output the result
		puts $fileid "Hybrid Domain Swing-out is $angle_degrees degrees. Reference Molecule is [molinfo [lindex $mols $zim] get name]. Comparison Molecule is [molinfo [lindex $mols $ye] get name]"
	
	}
	puts $fileid ""
	$bi delete; $hy1 delete; 
}
close $fileid
close $RMSDlog
close $errorlog
