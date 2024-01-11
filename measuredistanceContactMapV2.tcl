# Modified 12/12/2023 to be used for the reference states after minimization
# #### 12-12-2023 Ignore: To be used after you have extracted frames using matlab Metadynamics_extractframes.m and vmdExtractFrames.tcl. In that order.
set folders "1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/"
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
mol delete all
foreach x $folders {
	cd ${root}$x
	mol new step3_input.psf waitfor all
	mol addfile step4.0_minimization.gro waitfor all
	mol rename top $x
	#3380,10614
	#3380,10617
	#6988,10638
	#8492,10638
	#8542,10621
	#2837,10617
	set atoms1L [atomselect top "serial 3380 or serial 10614"];	set atoms1P [$atoms1L get index]
	set atoms2L [atomselect top "serial 3380 or serial 10617"];	set atoms2P [$atoms2L get index]
	set atoms3L [atomselect top "serial 6988 or serial 10638"];	set atoms3P [$atoms3L get index]
	set atoms4L [atomselect top "serial 8492 or serial 10638"];	set atoms4P [$atoms4L get index]
	set atoms5L [atomselect top "serial 8542 or serial 10621"];	set atoms5P [$atoms5L get index]
	set atoms6L [atomselect top "serial 2837 or serial 10617"];	set atoms6P [$atoms6L get index]
	set numframes [molinfo top get numframes]
	set fileid [open "RefstateDistances.dat" w+]; puts $fileid "Frame Atoms1-3380,10614 Atoms2-3380,10617 Atoms3-6988,10638 Atoms4-8492,10638 Atoms5-8542,10621 Atoms6-2837,10617"
	for {set x 0} {$x < $numframes} {incr x} {
		set atoms1distance [measure bond $atoms1P frame $x]; set atoms2distance [measure bond $atoms2P frame $x]; set atoms3distance [measure bond $atoms3P frame $x]
		set atoms4distance [measure bond $atoms4P frame $x]; set atoms5distance [measure bond $atoms5P frame $x]; set atoms6distance [measure bond $atoms6P frame $x]
		puts $fileid "$x $atoms1distance $atoms2distance $atoms3distance $atoms4distance $atoms5distance $atoms6distance"
		}
	close $fileid
	$atoms1L delete; $atoms2L delete; $atoms3L delete; $atoms4L delete; $atoms5L delete; $atoms6L delete
	unset atoms1L; unset atoms1P; unset atoms2L; unset atoms2P; unset atoms3L; unset atoms3P; unset atoms4L; unset atoms4P; unset atoms5L; unset atoms5P; unset atoms6L
	unset atoms6P; unset atoms1distance; unset atoms2distance; unset atoms3distance; unset atoms4distance; unset atoms5distance; unset atoms6distance
	#mol delete all
	}
cd $root