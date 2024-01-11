# For identifying contacts betwen RGD ligand and the integrin receptor.

# resname-resid-segname(only the last letter)
# segname PROA,etc.
# Resname VAL
# Resid




package require math::statistics;
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
foreach jkl {1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/} {
#foreach jkl {4.open35A/} {}
	cd ${root}${jkl}/cluster/
	mol delete all
	set pdbfiles XboundxtractCAT2.dcd
	foreach lek $pdbfiles {
		mol new ../step1_pdbreader.psf
		mol addfile ../${lek} waitfor all
		#mol addfile [lindex $pdbfiles 0] waitfor all
	}
	set loadedmols [molinfo list]
	set counter2 0
	foreach y $loadedmols {
		set counter [lindex $pdbfiles $counter2]		
		# I want to select only the heavy atoms that are within 5 of the ligand heavy atoms
		set selection1 [atomselect $y "(not hydrogen and within 5.0 of (residue 698 and not hydrogen)) and not (residue 695 to 701)" frame 0]
		set selection2 [atomselect $y "(not hydrogen and within 5.0 of (residue 699 and not hydrogen)) and not (residue 695 to 701)" frame 0]
		set selection3 [atomselect $y "(not hydrogen and within 5.0 of (residue 700 and not hydrogen)) and not (residue 695 to 701)" frame 0]
		set selection4 [atomselect $y "(not hydrogen and within 5.0 of ((residue 698 to 700) and not hydrogen)) and not (residue 695 to 701)" frame 0]
		set lig698 [atomselect $y "residue 698 and not hydrogen" frame 0]
		set lig699 [atomselect $y "residue 699 and not hydrogen" frame 0]
		set lig700 [atomselect $y "residue 700 and not hydrogen" frame 0]
		set fileid1 [open "BoundCs698${counter}.txt" w+]
		set fileid2 [open "BoundCs699${counter}.txt" w+]
		set fileid3 [open "BoundCs700${counter}.txt" w+]
		set fileid4 [open "log${counter}.log" w+]
		set fileid5 [open "means${counter}.log" w+]
		set fileid6 [open "BoundCTime${counter}.txt" w+]
		#set fileid7 [open "BoundCs699Lifetime${counter}.txt" w+]
		#set fileid8 [open "BoundCs700Lifetime${counter}.txt" w+]
		puts $fileid4 "$pdbfiles"
		#
		set numframes [molinfo $y get numframes]
		puts $fileid1 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 698 and not hydrogen)) and not (residue 695 to 701)
#Residue, MeanDistance, StandardDeviation"
		puts $fileid2 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 699 and not hydrogen)) and not (residue 695 to 701)
#Residue, MeanDistance, StandardDeviation"
		puts $fileid3 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 700 and not hydrogen)) and not (residue 695 to 701)
#Residue, MeanDistance, StandardDeviation"
		puts $fileid6 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 698 and not hydrogen)) and not (residue 695 to 701)"
		#puts $fileid7 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 699 and not hydrogen)) and not (residue 695 to 701)"
		#puts $fileid8 "#Total structures in cluster ${counter}: $numframes, Selection Criteria:, (not hydrogen and within 5.0 of (residue 700 and not hydrogen)) and not (residue 695 to 701)"
		puts $fileid4 "Numframes $numframes"
		array set data_select1 {}
		array set data_select2 {}
		array set data_select3 {}
		unset data_select1
		unset data_select2
		unset data_select3
		# go through all frames and build an array with the keys of the array as all the detected residues in this conformation
		for {set x 0} {$x < $numframes} {incr x} {
			$selection1 frame $x
			$selection2 frame $x
			$selection3 frame $x
			$selection4 frame $x
			$lig698 frame $x
			$lig699 frame $x
			$lig700 frame $x
			$selection1 update
			$selection2 update
			$selection3 update
			$selection4 update
			$lig698 update
			$lig699 update
			$lig700 update
			if {[$selection1 num] != 0} {
				set residues1 [lsort -unique [$selection1 get residue]]
			} else {
				set residues1 ""
			}
			if {[$selection2 num] != 0} {
				set residues2 [lsort -unique [$selection2 get residue]]
			} else {
				set residues2 ""
			}
			if {[$selection3 num] != 0} {
				set residues3 [lsort -unique [$selection3 get residue]]
			} else {
				set residues3 ""
			}
			if {[$selection4 num] != 0} {
				set residues4 [lsort -unique [$selection4 get residue]]
			} else {
				set residues4 ""
			}
			puts $fileid4 "### $x all detected residues are $residues1
$residues2
$residues3
			"
			# go through each residue
			foreach ki $residues1 {
				set r698 ""
				# select the heavy atoms of that residue
				set heavyatoms1 [atomselect $y "residue $ki and not hydrogen" frame $x]
				# measure the distance between all of those heavy atoms and the 3 ligand residues heavy atoms
				# and append all of those distances to a separate list
				puts $fileid4 "starting work on residue $ki"
				foreach ind [$heavyatoms1 get index ] {
					foreach ligind [$lig698 get index] {
						append r698 "[measure bond "$ind $ligind" frame $x] "
						#puts [measure bond "$ind $ligind"]
					}
				}
				$heavyatoms1 delete

				puts $fileid4 "all distances to R698 are $r698"
				# find the minimum distance in each ligand residues' list
				set val1 99999999999
				set lep1 ""
				foreach lep1 $r698 {
					if {$lep1 < $val1} {
						set val1 $lep1
					}
				}
				puts $fileid4 "#minimum is $val1"
				
				puts $fileid4 " "
				append data_select1($ki) " $val1"
				unset heavyatoms1; unset val1;  unset r698; 
			}
			
			unset residues1
			# go through each residue
			foreach kp $residues2 {
				set r699 ""
				# select the heavy atoms of that residue
				set heavyatoms2 [atomselect $y "residue $kp and not hydrogen" frame $x]
				# measure the distance between all of those heavy atoms and the 3 ligand residues heavy atoms
				# and append all of those distances to a list
				puts $fileid4 "starting work on residue $kp"
				foreach inds [$heavyatoms2 get index ] {
					foreach liginds [$lig699 get index] {
						append r699 "[measure bond "$inds $liginds" frame $x] "
					}
				}
				$heavyatoms2 delete
				# find the minimum distance in each ligand residues' list
				set val2 99999999
				set lep2 ""
				foreach lep2 $r699 {
					if {$lep2 < $val2} {
						set val2 $lep2
					}
				}
				puts $fileid4 "all distances to R699 are $r699"
				puts $fileid4 "#minimums are $val2"
				puts $fileid4 " "
				# Filter out too long minimum distances. if the min distance is <= 5.0 append it to an array where the keys are residue numbers
				append data_select2($kp) " $val2"
				unset heavyatoms2; unset val2;  unset r699; 
			}
			unset residues2
			# go through each residue			
			foreach kq $residues3 {
				set r700 ""
				# select the heavy atoms of that residue
				set heavyatoms3 [atomselect $y "residue $kq and not hydrogen" frame $x]
				# measure the distance between all of those heavy atoms and the 3 ligand residues heavy atoms
				# and append all of those distances to a list
				puts $fileid4 "starting work on residue $kq"
				foreach inda [$heavyatoms3 get index ] {
					foreach liginda [$lig700 get index] {
						append r700 "[measure bond "$inda $liginda" frame $x] "
					}
				}
				$heavyatoms3 delete
				# find the minimum distance in each ligand residues' atomlist
				set val3 99999999
				set lep3 ""
				foreach lep3 $r700 {
					if {$lep3 < $val3} {
						set val3 $lep3
					}
				}
				puts $fileid4 "all distances to R700 are $r700"
				puts $fileid4 "#minimums are $val3"
				puts $fileid4 " "
				# 
				append data_select3($kq) " $val3"
				unset heavyatoms3; unset val3;  unset r700; 
			}
			unset residues3
			# creating list of all of the residues that are near all 3 RGD ligand residues
			foreach kq $residues4 {
				append data_select4($kq) " 0"
			}
		}
		set resnidch1 ""
		foreach miv [array names data_select4] {
			set get1 [atomselect $y "residue $miv and name CA"]
			set resn1 [$get1 get resname]
			set seg1 [$get1 get segname]
			set resi1 [$get1 get resid] 
			$get1 delete
			append resnidch1 "${resn1}-${resi1}-[string index $seg1 end] "
		}
		
		
		
		# get the lifetimes
		puts $fileid6 "#Residue, Resname-Resid-Chain, frames"
		puts $fileid6 "$resnidch1"
		foreach pot1 [array names data_select1] {
			set getresidue1 [atomselect $y "residue $pot1 and name CA"]
			set resname1 [$getresidue1 get resname]
			set segname1 [$getresidue1 get segname]
			set resid1 [$getresidue1 get resid] 
			$getresidue1 delete
			set resnameidchain1 "${resname1}-${resid1}-[string index $segname1 end]"
			puts $fileid6 "$pot1, $resnameidchain1, [llength $data_select1($pot1)]"
			
		}
		puts $fileid6 ""
		puts $fileid6 "#Residue, Resname-Resid-Chain, frames"
		puts $fileid6 "$resnidch1"
		foreach pot2 [array names data_select2] {
			set getresidue2 [atomselect $y "residue $pot2 and name CA"]
			set resname2 [$getresidue2 get resname]
			set segname2 [$getresidue2 get segname]
			set resid2 [$getresidue2 get resid] 
			$getresidue2 delete
			set resnameidchain2 "${resname2}-${resid2}-[string index $segname2 end]"
			puts $fileid6 "${pot2}, ${resnameidchain2}, [llength $data_select2($pot2)]"
		}
		puts $fileid6 ""
		puts $fileid6 "#Residue, Resname-Resid-Chain, frames"
		puts $fileid6 "$resnidch1"
		foreach pot3 [array names data_select3] {
			set getresidue3 [atomselect $y "residue $pot3 and name CA"]
			set resname3 [$getresidue3 get resname]
			set segname3 [$getresidue3 get segname]
			set resid3 [$getresidue3 get resid] 
			$getresidue3 delete
			set resnameidchain3 "${resname3}-${resid3}-[string index $segname3 end]"
			puts $fileid6 "$pot3, $resnameidchain3, [llength $data_select3($pot3)]"
		}
		
		puts $fileid6 ""
#loop through all frames again to measure the distance of the detected residues throughout the captured frames.
		array set new_data1 {}
		array set new_data2 {}
		array set new_data3 {}
		unset new_data1
		unset new_data2
		unset new_data3
		for {set x 0} {$x < $numframes } {incr x} {
			foreach mno [array names data_select1] {
				set dist1 ""
				set newselection1 [atomselect $y "residue $mno and not hydrogen" frame $x]
				foreach nop [$newselection1 get index] {
					foreach opq [$lig698 get index] {
						lappend dist1 [measure bond "$nop $opq" frame $x]
					}
				}
				#find minimum for this residue 
				set comp1 9999999
				foreach pqr $dist1 {
					if {$pqr < $comp1} {
						set comp1 $pqr
					}
				}
				append new_data1($mno) " $comp1"
				$newselection1 delete
			}
			foreach qrs [array names data_select2] {
				set dist2 ""
				set newselection2 [atomselect $y "residue $qrs and not hydrogen" frame $x]
				foreach rst [$newselection2 get index] {
					foreach stu [$lig699 get index] {
						lappend dist2 [measure bond "$stu $rst" frame $x]
					}
				}
				#find minimum for this residue 
				set comp2 9999999
				foreach tuv $dist2 {
					if {$tuv < $comp2} {
						set comp2 $tuv
					}
				}
				append new_data2($qrs) " $comp2"
				$newselection2 delete
			}
			foreach uvw [array names data_select3] {
				set dist3 ""
				set newselection3 [atomselect $y "residue $uvw and not hydrogen" frame $x]
				foreach vwx [$newselection3 get index] {
					foreach wxy [$lig700 get index] {
						lappend dist3 [measure bond "$vwx $wxy" frame $x]
					}
				}
				#find minimum for this residue 
				set comp3 9999999
				foreach xyz $dist3 {
					if {$xyz < $comp3} {
						set comp3 $xyz
					}
				}
				append new_data3($uvw) " $comp3"
				$newselection3 delete
			}
		}
		puts "finished all frames
starting statistical calculations"		
		# now that all frames have been gone through 
		# find the mean value for the distance from each residue 
		#I also need to figure out how to output all the values without finding the mean so I can do statistics or do the stats here.
		puts $fileid5 "$x [array get new_data1]"
		puts [array names new_data1]
		puts [array names new_data2]
		puts [array names new_data3]
		set list1 ""; set list2 ""; set list3 ""
		set output1 ""; set output2 ""; set output3 ""
		set stdev1 ""; set stdev2 ""; set stdev3 ""
		foreach jikl [ array names new_data1 ] {
			set list1 $new_data1($jikl)
			puts $fileid4 "$jikl $list1"
			set output1 [expr { [tcl::mathop::+ {*}$list1 0.0] / max(1, [llength $list1])}]
			set stdev1 [::math::statistics::stdev $list1]
			puts $fileid1 "$jikl $output1 $stdev1"
		}
		foreach jkli [ array names new_data2 ] {
			set list2 $new_data2($jkli)
			puts $fileid4 "$jkli $list2"
			set output2 [expr { [tcl::mathop::+ {*}$list2 0.0] / max(1, [llength $list2])}]
			set stdev2 [::math::statistics::stdev $list2]
			puts $fileid2 "$jkli $output2 $stdev2"
		}
		foreach jlki [ array names new_data3 ] {
			set list3 $new_data3($jlki)
			puts $fileid4 "$jlki $list3"
			set output3 [expr { [tcl::mathop::+ {*}$list3 0.0] / max(1, [llength $list3])}]
			set stdev3 [::math::statistics::stdev $list3]
			puts $fileid3 "$jlki $output3 $stdev3"
		}	
	incr counter2
	$selection1 delete; $selection2 delete; $selection3 delete; $selection4 delete; $lig698 delete; $lig699 delete; $lig700 delete
	unset list1; unset list2; unset list3; unset output1; unset output2; unset output3; unset lig698; unset lig699; unset lig700; unset selection1
	close $fileid1; close $fileid2; close $fileid3; close $fileid4; close $fileid5; close $fileid6; #close $fileid7; close $fileid8
	}
 
}
#puts [info local]
#foreach var [info local] {
#unset $var
#}
cd C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/
return

