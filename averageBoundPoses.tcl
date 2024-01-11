#
# for measuring the distance between atoms in the open or closed conformation and puts all gathered data into one file. These are the indices determined by ref.pdb
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
set folders "1.bent35A 2.int135A 3.int235A 4.open35A"
#set pdbs "boundxtractCAT.pdb unboundxtractCAT.pdb"
mol delete all
foreach y $folders {
	cd ${root}${y}
	mol new step1_pdbreader.psf
	for {set x 1} {$x <= 16} {incr x} {
		if { [file exists pose${x}.dcd ] == 0 } {
			puts "couldn't find file pose${x}.dcd"
			continue
		}
		set existfile [catch {file size pose${x}.dcd} size]
		if {$existfile != 0 || $size < 3} {
			puts "pose${x}.dcd file is empty in folder $y"
			continue
		} else {
			mol addfile pose${x}.dcd waitfor all
			set selection1 [atomselect top "residue 0 to 701"]
			set avgpos [measure avpos $selection1 first 1 last -1]
			$selection1 set {x y z} $avgpos
			$selection1 writepdb boundpose${x}.pdb
			$selection1 delete
		}
	}	
}
# NOTES:			
# A Valen after the D (Asp) residue (if present) consistently coordinates bonding with ADMIDAS Ca2+ through a water molecule in Beta3.			
# alphaIIb ASP-232 sidechain holds in place two water molecules which H-bonds to Arginine			
# The conformation of the ligand backbone is stabilized as it crosses the interface between the alphaIIb and Beta3 subunits by two previously unremarked water molecules with strong density that are held in place by hydrogen bonds to the side chain of alphaIIb Asp-232 and  backbone of beta3 Ala-218 and that hydrogen bond to the carbonyl oxygen of the Arg of RGD ( Fig. 3 b ).			

# residue 698 to 700
