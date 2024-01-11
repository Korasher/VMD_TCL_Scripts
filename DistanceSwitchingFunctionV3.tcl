# The transformed distance can be compared with a reference value in order to calculate the squared distance between two contact maps.

proc CMdist {refdist currentdist} {
	set cmdist [expr {pow(($refdist - $currentdist),2)}]
	return $cmdist
}
proc square {dist} {
	set square [expr {pow(($dist),2)}]
	return $square
}
proc switching {framenum atomselection1 atomselection2} {
	$atomselection1 frame $framenum
	$atomselection2 frame $framenum
	set atoms2 [veclength [vecsub [measure center $atomselection1] [measure center $atomselection2]]]
	set switch1 [expr { (1-pow((($atoms2 - 0)/15),6))/(1-pow((($atoms2 - 0)/15),12)) }]
	return $switch1
	
}
proc refdist {frame atomselection3 atomselection4} {
	$atomselection3 frame $frame
	$atomselection4 frame $frame
	set atoms1 [veclength [vecsub [measure center $atomselection3] [measure center $atomselection4]]]
	return $atoms1
}
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
set folders "1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/"
set subfolders "1 2 3 4 5 6 7 8"
set rswitch1 0.982227
set rswitch2 0.997942
set rswitch3 0.999417
set rswitch4 0.688831
set rswitch5 0.783917
set rswitch6 0.999967
set rsquare1 [square $rswitch1]
set rsquare2 [square $rswitch2]
set rsquare3 [square $rswitch3]
set rsquare4 [square $rswitch4]
set rsquare5 [square $rswitch5]
set rsquare6 [square $rswitch6]
set refsum [expr {pow(($rswitch1 + $rswitch2 + $rswitch3 + $rswitch4 + $rswitch5 + $rswitch6),2)}]
set squarerefsum [expr { $rsquare1 + $rsquare2 + $rsquare3 + $rsquare4 + $rsquare5 + $rsquare6}]
foreach x $folders {
	cd ${root}${x}
	set fid [open switchv123.csv w+]
	puts $fid $squarerefsum
	mol delete all
	mol new ref.pdb waitfor all
	set atoms11 [atomselect top "serial 3380"]
	set atoms12 [atomselect top "serial 10614"]
	set atoms21 $atoms11
	set atoms22 [atomselect top "serial 10617"]
	set atoms31 [atomselect top "serial 6988"]
	set atoms32 [atomselect top "serial 10638"]
	set atoms41 [atomselect top "serial 8492"]
	set atoms42 $atoms32
	set atoms51 [atomselect top "serial 8542"]
	set atoms52 [atomselect top "serial 10621"]
	set atoms61 [atomselect top "serial 2837"]
	set atoms62 $atoms22
	foreach y $subfolders {
		if { [catch {cd ${root}${x}${y}} ] } {
			puts "couldn't find folder ${x}${y}"
			continue
		}
		mol addfile trajout.xtc waitfor all
		set numframes [molinfo top get numframes]
		puts $numframes
		for {set z 1} {$z <= $numframes} {incr z} {
			set switch1 [switching $z $atoms11 $atoms12]
			set switch2 [switching $z $atoms21 $atoms22]
			set switch3 [switching $z $atoms31 $atoms32]
			set switch4 [switching $z $atoms41 $atoms42]
			set switch5 [switching $z $atoms51 $atoms52]
			set switch6 [switching $z $atoms61 $atoms62]
		
			set cursum [expr {pow(($switch1 + $switch2 + $switch3 + $switch4 + $switch5 + $switch6),2)}]
			set squaresum [expr { [square $switch1] + [square $switch2] + [square $switch3] + [square $switch4] + [square $switch5] + [square $switch6]}]
			#set sqdiff1 [expr { [square $switch1] - $rsquare1}]
			#set sqdiff2 [expr { [square $switch2] - $rsquare2}]
			#set sqdiff3 [expr { [square $switch3] - $rsquare3}]
			#set sqdiff4 [expr { [square $switch4] - $rsquare4}]
			#set sqdiff5 [expr { [square $switch5] - $rsquare5}]
			#set sqdiff6 [expr { [square $switch6] - $rsquare6}]
			
			# 1 sum all(square individual transform) - sum all(square individual ref transform)
			set output1 [expr {$squarerefsum - $squaresum}]
			# 2 ((sum(all current transform))^2) - ((sum(ref transform))^2)
			#set output2 [expr {$cursum - $refsum}]
			# 3 sumall(Square Individual - Square Ref)
			#set output3 [expr {$sqdiff1 + $sqdiff2 + $sqdiff3 + $sqdiff4 + $sqdiff5 + $sqdiff6}]
			puts $fid "${output1}"
		}
		animate delete beg 1 end -1
	}
	$atoms11 delete
	$atoms12 delete
	$atoms22 delete
	$atoms31 delete
	$atoms32 delete
	$atoms41 delete
	$atoms51 delete
	$atoms52 delete
	$atoms61 delete
	close $fid
	mol delete all
}
cd $root
return