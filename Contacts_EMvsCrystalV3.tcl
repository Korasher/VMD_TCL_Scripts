# for measuring the distance between atoms in the open or closed conformation and puts all gathered data into one file. These are the indices determined by ref.pdb
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
set folders "1.bent35A 2.int135A 3.int235A 4.open35A"
mol delete all
foreach x $folders {
	cd ${root}${x}
	mol new step1_pdbreader.psf
	mol addfile allcat.xtc waitfor all
	mol rename top $x
	set fileid [open "Contacts.csv" w+]
	puts $fileid "Crystal Structure Contact Analysis
	Name,Frame,a224to408_1,a224to408_2,a224to408_3,a224to408_4,a189to408_1,a189to408_2,b122to410_1,b122to410_2,b123to410_1,b123to410_2,b214to410_1,b214to410_2,b214to410_3,b214to410_4,b215to410_1,b215to410_2,b215to410_3,b215to410_4,b215to410_5,b215to410_6"
	puts -nonewline $fileid "[molinfo top get name],"
	set numframes [molinfo top get numframes]
	set 214 [atomselect top "index 8472 or index 8478 or index 8477 or index 8492"]
	set 410_1 [atomselect top "index 10640"]
	set 410_2 [atomselect top "index 10641"]
	for {set x 1} {$x <= $numframes} {incr x} {
		$214 frame $x
		$410_1 frame $x
		$410_2 frame $x
		set 224to408_1 [measure bond {3380 10617} frame $x];   # alpha Asp-224 OD1 to NH1 of Arg-408
		set 224to408_2 [measure bond {3381 10620} frame $x];   # alpha ASP-224 OD2 to NH2 of Arg-408
		set 224to408_3 [measure bond {3381 10617} frame $x];	  # alpha ASP-224 OD2 to NH1 of Arg-408
		set 224to408_4 [measure bond {3380 10620} frame $x];   # alpha ASP-224 OD1 to NH2 of Arg-408
		set 189to408_1 [measure bond {2836 10617} frame $x];	  # alpha Tyr-189 O to NH1 of Arg-408
		set 189to408_2 [measure bond {2836 10620} frame $x];	  # alpha Tyr-189 O to NH2 of Arg-408
		set 122to410_1 [measure bond {6990 10641} frame $x];   # beta  Tyr-122 N   to OD2 Asp-410
		set 122to410_2 [measure bond {6990 10640} frame $x];   # beta  Tyr-122 N   to OD1 Asp-410
		set 123to410_1 [measure bond {7011 10640} frame $x];   # beta  Ser-123 N   to OD1 of ASP-410
		set 123to410_2 [measure bond {7011 10641} frame $x];   # beta  Ser-123 N   to OD2 of ASP-410		
		set 214to410_1 [measure bond {8470 10641} frame $x]; 	  								    # beta Arg-214 N   					   to OD2 Asp-410
		set 214to410_2 [measure bond {8470 10640} frame $x]; 	  								    # beta Arg-214 N   					   to OD1 Asp-410
		set 214to410_3 [veclength [vecsub [measure center $214 weight mass] [measure center $410_1 weight mass]]] ; # beta Arg-214 COM of CA, CB, CG and C to OD1 Asp-410
		set 214to410_4 [veclength [vecsub [measure center $214 weight mass] [measure center $410_2 weight mass]]] ; # beta Arg-214 COM of CA, CB, CG and C to OD2 Asp-410	
		set 215to410_1 [measure bond {8494 10641} frame $x];  # beta ASN-215 N   to OD2 ASP-410 (Backbone. Going to assume this is the incorrect contact due to comment in Springer 2008 that "Asn-215 should probably be flipped")
		set 215to410_2 [measure bond {8494 10640} frame $x] ; # beta ASN-215 N   to OD1 ASP-410 	
		set 215to410_3 [measure bond {8503 10641} frame $x] ; # beta ASN-215 ND2 to OD2 ASP-410 
		set 215to410_4 [measure bond {8503 10640} frame $x] ; # beta ASN-215 ND2 to OD1 ASP-410 
		set 215to410_5 [measure bond {8498 10636} frame $x] ;  # beta Asn-215 CB  to CB  Asp-410
		set 215to410_6 [measure bond {8498 10639} frame $x] ;  # beta Asn-215 CB  to CG  Asp-410	
		puts $fileid "${x},${224to408_1},${224to408_2},${224to408_3},${224to408_4},${189to408_1},${189to408_2},${122to410_1},${122to410_2},${123to410_1},${123to410_2},${214to410_1},${214to410_2},${214to410_3},${214to410_4},${215to410_1},${215to410_2},${215to410_3},${215to410_4},${215to410_5},${215to410_6}"
	}
	$214 delete
	$410_1 delete
	$410_2 delete
	close $fileid
	mol delete all
}
# NOTES:			
# A Valen after the D (Asp) residue (if present) consistently coordinates bonding with ADMIDAS Ca2+ through a water molecule in Beta3.			
# alphaIIb ASP-232 sidechain holds in place two water molecules which H-bonds to Arginine			
# The conformation of the ligand backbone is stabilized as it crosses the interface between the alphaIIb and Beta3 subunits by two previously unremarked water molecules with strong density that are held in place by hydrogen bonds to the side chain of alphaIIb Asp-232 and  backbone of beta3 Ala-218 and that hydrogen bond to the carbonyl oxygen of the Arg of RGD ( Fig. 3 b ).			

# residue 698 to 700
