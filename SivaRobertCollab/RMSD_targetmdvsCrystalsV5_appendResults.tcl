# Created July 15, 2024 by Robert Coffman
# this script compares residues that are common between the desired reference structure and all the loaded structures. It always compares the first chain 
# to the first chain in the structure file and so on, regardless of how they are named. I.e. chain I chain A segname PROA etc

# this script will not work on PDB files downloaded straight from rcsb.org due to multiple copies of residues in the files.
# I used charmm-gui to extract the protein fragments I wanted. Other software can probably be used to prepare the pdb files as long as unique chain or segnames are written
# for the different protein chains.

# starting at the end of hte loaded molecules. How many do you want to compare to all the rest? 0 = only last molecule, 1 = last molecule and next to last molecule..etc.
set number 18
# adjust these as necessary
# where is the iteration file?
set root "D:/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/string/iter0_onlyprotein/iter0"
# assuming the subfolders of iter* are named by increasing number
set subfolders "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
# where are the pdb files that you prepared for this script to use? This is also the path where the RMSD output will be written to and
# where directories will be made for the aligned structures to be written to.
set otherroot "D:/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/PDBstructures_2"



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
# compare molecule1 to molecule2 and return a nested list of two selection texts with the resids they have in common between the respective chains.
proc getcommonresids {mol1 mol2} {
	global errorlog
	set chainBmasterlist ""
	set chainAmasterlist ""
	set test [atomselect $mol1 "chain A and alpha"]
	set test2 [atomselect $mol1 "chain B and alpha"]
	set test3 [atomselect $mol2 "chain A and alpha"]
	set test4 [atomselect $mol2 "chain B and alpha"]
	set l1 [$test get resid]
	set l2 [$test2 get resid]
	set l3 [$test3 get resid]
	set l4 [$test4 get resid]
	foreach res $l1 {
		foreach res3 $l3 {
			if {$res == $res3} {
				lappend chainAmasterlist $res
			}
		}
	}
	foreach res2 $l2 {
		foreach res4 $l4 {
			if {$res2 == $res4} {
				lappend chainBmasterlist $res2
			}
		}
	}
	if { [string length $chainAmasterlist] == 0 } {
		set chainAmasterlist ""
	}
	lsort -unique $chainAmasterlist
	if { [string length $chainBmasterlist] > 0 && [string length $chainAmasterlist] > 0} {
		lsort -unique $chainBmasterlist
		set commonresids [list $chainAmasterlist $chainBmasterlist]
		unset chainBmasterlist
		set existsbchain 1
		set existsachain 1
	} elseif { [string length $chainAmasterlist] > 0 && [string length $chainBmasterlist] == 0 } { 
		set commonresids [list $chainAmasterlist]
		set existsachain 1
		set existsbchain 0
	} elseif { [string length $chainAmasterlist] == 0 && [string length $chainBmasterlist] > 0 } {
		set commonresids [list $chainBmasterlist]
		set existsachain 0
		set existsbchain 1
	} else {
		set existsbchain 0
		set existsachain 0
	}
	# compile a selection text so that I can align these resids and get an RMSD.
	if { $existsachain == 1 && $existsbchain == 1} {
		# create the beginning of the selection text for the first chain
		set selectAtext "resid [lindex [lindex $commonresids 0] 0]"
		# append to the selection text for however many resids there are
		for { set ya 1} {$ya < [llength [lindex $commonresids 0]]} {incr ya} {
			append selectAtext " or resid [lindex [lindex $commonresids 0] $ya]"
		}
		set selectBtext "resid [lindex [lindex $commonresids 1] 0]"
		# append to the selection text for however many resids there are
		for { set yz 1} {$yz < [llength [lindex $commonresids 1]]} {incr yz} {
			append selectBtext " or resid [lindex [lindex $commonresids 1] $yz]"
		}	
		# create the final text to use for selection
		set finalselecttext "alpha and ((chain A and ( $selectAtext )) or (chain B and ( $selectBtext )))"
		set finalalltext "((chain A and ( $selectAtext )) or (chain B and ( $selectBtext )))"
	# create the beginning of the selection text for the second chain
	} elseif { $existsachain == 1 && $existsbchain == 0 } {
		# create the beginning of the selection text for the first chain
		set selectAtext "resid [lindex [lindex $commonresids 0] 0]"
		# append to the selection text for however many resids there are
		for { set ya 1} {$ya < [llength [lindex $commonresids 0]]} {incr ya} {
			append selectAtext " or resid [lindex [lindex $commonresids 0] $ya]"
		}
		set finalselecttext "alpha and ((chain A and ( $selectAtext )) )"
		set finalalltext "((chain A and ( $selectAtext )) )"
		#puts "didn't find common resids for chain B of molids $mol1 ${mol2}. Double check that this is accurate"
	} elseif { $existsachain == 0 && $existsbchain == 1} {
		set selectBtext "resid [lindex [lindex $commonresids 1] 0]"
		# append to the selection text for however many resids there are
		for { set yz 1} {$yz < [llength [lindex $commonresids 1]]} {incr yz} {
			append selectBtext " or resid [lindex [lindex $commonresids 1] $yz]"
		}	
		# create the final text to use for selection
		set finalselecttext "alpha and (chain B and ( $selectBtext ))"
		set finalalltext "chain B and ( $selectBtext )"
	} else { 
		set finalselecttext ""
		set finalalltext ""
		puts $errorlog "didn't find common resids for chain A or chain B for molids $mol1 ${mol2}. Double check that this is accurate"
	}
	return [list $finalselecttext $finalalltext]
}
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
# align two molecules and return the RMSD
proc getRMSD {finalselecttext finalalltext ye} {
	global refmol; global mols; global errorlog	
	# select a molecule to compare the other molecules to
	set align1 [atomselect [lindex $mols $refmol] $finalselecttext]
	# make the comparison reporting the RMSD to stdout and saving it to a variable.
	set align2 [atomselect [lindex $mols $ye] $finalselecttext]
	if { [$align1 num] != [$align2 num] } {
		puts $errorlog "the atom selections did not select the same number of atoms. This should not happen! molid [lindex $mols $ye] vs [lindex $mols $refmol]"
		puts $errorlog "Reference atom count [$align1 num] Experiment atom count [$align2 num]"
		puts $errorlog "$finalselecttext"
		return NaN
		}
	set fit [measure fit $align2 $align1]
	set all [atomselect [lindex $mols $ye] all]
	$all move $fit
	$align2 update
	
	set pdbname [split [molinfo [lindex $mols $ye] get name] .]
	set refpdb [split [molinfo [lindex $mols $refmol] get name] .]
	animate write pdb "./[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $pdbname 0]_aligned_CA.pdb"  sel [atomselect [lindex $mols $ye] $finalselecttext] [lindex $mols $ye]
	animate write psf "./[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $pdbname 0]_aligned_CA.psf"  sel [atomselect [lindex $mols $ye] $finalselecttext] [lindex $mols $ye]
	animate write pdb "./[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $refpdb 0]_[lindex $pdbname 0]_CA.pdb" sel [atomselect [lindex $mols $refmol] $finalselecttext] [lindex $mols $refmol]
	
	#animate write pdb "./all_atom_aligned_structures/[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $pdbname 0]_aligned.pdb" sel [atomselect [lindex $mols $ye] $finalalltext] [lindex $mols $ye]
	#animate write psf "./all_atom_aligned_structures/[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $pdbname 0]_aligned.psf"  sel [atomselect [lindex $mols $ye] $finalalltext] [lindex $mols $ye]
	#animate write pdb "./all_atom_aligned_structures/[molinfo [lindex $mols $refmol] get name]_alignedStructures/[lindex $refpdb 0]_[lindex $pdbname 0].pdb" sel [atomselect [lindex $mols $refmol] $finalalltext] [lindex $mols $refmol]
	#puts "[molinfo [lindex $mols $ye] get name] vs $refpdb RMSD is [measure rmsd $align1 $align2]"
	set rmsds " [measure rmsd $align1 $align2]"
	$align2 delete
	$all delete
	$align1 delete
	return $rmsds
}


# beginning of work code
mol delete all
set RMSDlog [open RMSDlog.log w]
set errorlog [open errorlog.err w]
# load the Crystal, CryoEM and NMR experimental structures
cd "${otherroot}"
puts $RMSDlog "working directory is [pwd]"
set exp_pdb [ls *.pdb]
set exp_psf [ls *.psf]
set exp_pdb [lreplace $exp_pdb 0 0]
set exp_psf [lreplace $exp_psf 0 0]
foreach pdbs $exp_pdb psfs $exp_psf {
	mol new $psfs waitfor all
	mol addfile $pdbs waitfor all
	mol rename top $pdbs
}
#set exp_pdb [ls 1tye*.pdb]
#set exp_psf [ls 1tye*.psf]
#set exp_pdb [lreplace $exp_pdb 0 0]
#set exp_psf [lreplace $exp_psf 0 0]
#foreach pdbs $exp_pdb psfs $exp_psf {
#	mol new $psfs waitfor all
#	mol addfile $pdbs waitfor all
#	mol rename top $pdbs
#}
#set exp_pdb [ls 2v*.pdb]
#set exp_psf [ls 2v*.psf]
#set exp_pdb [lreplace $exp_pdb 0 0]
#set exp_psf [lreplace $exp_psf 0 0]
#foreach pdbs $exp_pdb psfs $exp_psf {
#	mol new $psfs waitfor all
#	mol addfile $pdbs waitfor all
#	mol rename top $pdbs
#}
#mol new 1kup_model_1.psf
#mol addfile 1kup_model_1.pdb

# Remnant from comparing to the previous papers predicted structures.
#set otherroot2 "C:/Users/and_r/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/BiophysJ_Fig8_abf_PDB_string-20240717T200842Z-001/BiophysJ_Fig8_abf_PDB_string/bent-int1/structs"
#set othersubfolder2 "0 1 2 3 4 5 6 7 8 9"
#set otherroot3 "C:/Users/and_r/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/BiophysJ_Fig8_abf_PDB_string-20240717T200842Z-001/BiophysJ_Fig8_abf_PDB_string/bent-int2/structs"
#set otherroot4 "C:/Users/and_r/OneDrive - University of Utah/BidoneLab/integrin/siva_robertCollab/BiophysJ_Fig8_abf_PDB_string-20240717T200842Z-001/BiophysJ_Fig8_abf_PDB_string/bent-open/structs"
#foreach it $othersubfolder2 {
#	cd "${otherroot2}"
#	mol new run-tmd-aa-bent-int1-${it}_v2.pdb waitfor all
#	mol rename top run-tmd-aa-bent-int1-${it}
#}
#foreach it $othersubfolder2 {
#	cd "${otherroot3}"
#	mol new run-tmd-aa-bent-int2-${it}_v2.pdb waitfor all
#	mol rename top run-tmd-aa-bent-int2-${it}
#}
#foreach it $othersubfolder2 {
#	cd "${otherroot4}"
#	mol new run-tmd-aa-bent-open-${it}_v2.pdb waitfor all
#	mol rename top run-tmd-aa-bent-open-${it}
#}

# load the target_md pdbs from Siva
foreach x $subfolders {
	cd "${root}/target_md_${x}"
	mol new minimize-protein.pdb waitfor all
	mol rename top "target_md_${x}"
}
# get a list of all the loaded molecules
set mols [molinfo list]
# rename the first and second chains or segnames to A and B or PROA and PROB respectively
foreach y $mols {
	set chains [getchain_segnames $y]
	setchainnames $chains $y
}
cd $root
# start the actual comparisons (align, RMSD, writepdb)
for {set zim 0 } {$zim <= $number} {incr zim} {
	set refmol end-${zim}
	file mkdir "[molinfo [lindex $mols $refmol] get name]_alignedStructures"
	file mkdir "./all_atom_aligned_structures/[molinfo [lindex $mols $refmol] get name]_alignedStructures"
	# compare the resids in the desired reference molecule to the comparison molecule in a loop until all molecules have been compared to the reference.
	for { set ye 0} {$ye < [expr {[llength $mols] - ($number + 1) }] } {incr ye} {
		puts $RMSDlog "comparing molid [lindex $mols $ye] to molid [lindex $mols $refmol ]"
		# get the selection text that uses the resids that are common between the two molecules.
		set finalselecttext [lindex [getcommonresids [lindex $mols $ye] [lindex $mols $refmol ]] 0]
		set finalalltext [lindex [getcommonresids [lindex $mols $ye] [lindex $mols $refmol ]] 1]
		if { [string length $finalselecttext] > 0 } {
			# calculate the RMSD and write the aligned PDBs out
			append rmsds " [molinfo [lindex $mols $refmol] get name] vs [molinfo [lindex $mols $ye] get name] [getRMSD $finalselecttext $finalalltext $ye] \n"
			unset finalselecttext
			unset finalalltext
		} else {
			append rmsds " [molinfo [lindex $mols $refmol] get name] vs [molinfo [lindex $mols $ye] get name] No Common Residues \n"
			unset finalselecttext
			unset finalalltext
		}	
	}
	cd $root 
	set fileid [open RMSD_[molinfo [lindex $mols $refmol] get name]_v3.dat w]
	puts $fileid $rmsds
	close $fileid
	unset rmsds
}
close $RMSDlog
close $errorlog
puts "finished"
# NOTES: 
# The second and third versions of integrins were extracted from the following pdb files on 10/1/2024 and used in the analysis.
# 3fcs has two versions of integrin used the first.
# 3fcu has three versions of integrin used the first
# 3nid, 3nif and 3nig, 3t3m, 3t3p, 4z7n, 4z7o, 4z7q, 5hdb, 7tct, 7td8, 7tho, 7tmz, 7tpd, 7u9f, 7u9v, 
# 7uk9, 7u60, 7ubr, 7ucy, 7udg, 7udh, 7ueo, 7ufh, 7uh8, 7uje, 7ujk, 7uk9, 7uko, 7ukp, 7ukt has two versions of integrin used the first
