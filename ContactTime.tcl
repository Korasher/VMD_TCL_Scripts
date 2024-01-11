



package require math::statistics;
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
# load all of the pdb cluster files across all conformations.
foreach jkl {1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/} {
	cd ${root}${jkl}/cluster/
	set pdbfiles [glob clusters.pdb*.pdb]
	foreach lek $pdbfiles {
		mol new ../step1_pdbreader.psf
		mol addfile $lek waitfor all
	}
}
set loadedmols [molinfo list]
set counter2 0
array set data_select4 {}
unset data_select4
# build master list of all the residues that are detected across all clusters and conformations.
foreach y $loadedmols {
	set selection4 [atomselect $y "(not hydrogen and within 5.0 of ((residue 698 to 700) and not hydrogen)) and not (residue 695 to 701)" frame 0]
	set numframes [molinfo $y get numframes]
	for {set x 0} {$x < $numframes} {incr x} {
		$selection4 frame $x
		$selection4 update
		if {[$selection4 num] != 0} {
			set residues4 [lsort -unique [$selection4 get residue]]
		} else {
				set residues4 ""
		}
		foreach kq $residues4 {
			append data_select4($kq) " 0"
		}
	}
	$selection4 delete
}
# Master list done
set resnidch1 ""
# Now get the Resname-Resid-Chain formatting done and output to a file.
foreach miv [array names data_select4] {
	set get1 [atomselect top "residue $miv and name CA"]
	set resn1 [$get1 get resname]
	set seg1 [$get1 get segname]
	set resi1 [$get1 get resid] 
	$get1 delete
	append resnidch1 "${resn1}-${resi1}-[string index $seg1 end] "
}
cd $root
set fileid1 [open "BoundContactsAllConformers.txt" w+]
puts $fileid1 $resnidch1
close $fileid1
# load the conformation trajectories' one by one
foreach jkm {1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/} {
	cd ${root}${jkm}/cluster/
	mol delete all
	set pdbfiles [glob clusters.pdb*.pdb]
	foreach lek $pdbfiles {
		mol new ../step1_pdbreader.psf
		mol addfile $lek waitfor all
	}
	set loadedmols [molinfo list]
	set outfile [string replace "$jkm" 9 9 ]
	set fileid2 [open "${outfile}.txt" w+]	
	set totalframes 0
	array set data_select4 {}
	unset data_select4
	# go through the conformations' loaded molecules (trajectories)
	# keep track of how many frames are in all 3 clusters and create the array list of residues and 0's to count the length later.
	foreach y $loadedmols {
		set selection4 [atomselect $y "(not hydrogen and within 5.0 of ((residue 698 to 700) and not hydrogen)) and not (residue 695 to 701)" frame 0]
		set numframes [molinfo $y get numframes]
		set totalframes [expr {$totalframes + $numframes}]
		# 
		for {set x 0} {$x < $numframes} {incr x} {
			$selection4 frame $x
			$selection4 update
			if {[$selection4 num] != 0} {
				set residues4 [lsort -unique [$selection4 get residue]]
			} else {
				set residues4 ""
			}
			foreach kq $residues4 {
				append data_select4($kq) " 0"
			}
		}
		$selection4 delete
	}
	puts $fileid2 "#Total structures in clusters: $totalframes, Selection Criteria:, (not hydrogen and within 5.0 of ((residue 698 to 700) and not hydrogen)) and not (residue 695 to 701)"
	puts $fileid2 "#Residue, Resname-Resid-Chain, Contact Time (% Frames)"
	# now that the list is complete with the number of zeros for the number of frames the residue was in contact
	# output that information to a file.
	foreach pot1 [array names data_select4] {
		set getresidue1 [atomselect $y "residue $pot1 and name CA"]
		set resname1 [$getresidue1 get resname]
		set segname1 [$getresidue1 get segname]
		set resid1 [$getresidue1 get resid] 
		$getresidue1 delete
		set resnameidchain1 "${resname1}-${resid1}-[string index $segname1 end]"
		puts $fileid2 "$pot1, $resnameidchain1, [expr {([llength $data_select4($pot1)] / double($totalframes)) * 100}]"
	}
	close $fileid2
}
cd $root
cd ..