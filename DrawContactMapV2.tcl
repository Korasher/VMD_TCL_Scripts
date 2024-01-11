# Modified 12/12/2023 to be used for the reference states after minimization
# #### 12-12-2023 Ignore: To be used after you have extracted frames using matlab Metadynamics_extractframes.m and vmdExtractFrames.tcl. In that order.
set folders "1.bent35A 2.int135A 3.int235A 4.open35A"
set root1 "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
set root2 "G:/My Drive/PaperFiguresV2/Figure3/"
mol delete all
cd ${root1}
play ${root2}CV2ContactMapBondsV5HiResV2.vmd
set loadedmols [molinfo list]
for {set o 1} {$o <=3} {incr o} {
	mol delete [lindex $loadedmols $o]
}

display resetview
graphics top delete all

#set dcontactsl "3382 3382 6991 8495 8545 2837"
#set dcontactsp "10618 10621 10642 10642 10625 10621"
#set acontacsl "3380 3380 6988 8492 8542 2837"
#set acontactsp "10614 10617 10638 10638 10621 10617"

set dcontactsl "8545"
set dcontactsp "10625"
set acontacsl "8542"
set acontactsp "10621"
scale by 2.5

foreach w $dcontactsl x $dcontactsp y $acontacsl z $acontactsp {
	# Desired
	display resetview
	scale by 2.5
	graphics top delete all
	set atoms1L [atomselect top "serial $w"]; set atoms2L [atomselect top "serial $x"];
	set x1 [$atoms1L get {x y z}]; set y1 [$atoms2L get {x y z}]
	graphics top color blue
	graphics top line [lindex $x1 0] [lindex $y1 0] width 4 style dashed
	# Actual
	set atoms3L [atomselect top "serial $y"]; set atoms4L [atomselect top "serial $z"];
	set x2 [$atoms3L get {x y z}]; set y2 [$atoms4L get {x y z}]
	graphics top color red
	graphics top line [lindex $x2 0] [lindex $y2 0] width 4 style dashed
	$atoms1L delete; $atoms2L delete; $atoms3L delete ; $atoms4L delete
	render snapshot ${w}_${x}_1.bmp
	rotate x by 30
	render snapshot ${w}_${x}_2.bmp
	rotate x by -30
	rotate y by 30
	render snapshot ${w}_${x}_3.bmp
	rotate y by -30
}
# ACTUAL CONTACTS
	#3380,10614 
	#3380,10617 
	#6988,10638
	#8492,10638
	#8542,10621
	#2837,10617
# DESIRED CONTACTS
	#3382,10618	(Pretty far off but might facilitate the desired contact)
	#3382,10621	(Pretty close)
	#6991,10642	(Pretty far off)
	#8495,10642	(Pretty far off)
	#8545,10625	 (Can't tell from pictures)
	#2837,10621	(Pretty darn Close)
	
