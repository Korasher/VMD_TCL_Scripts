# For measuring salt bridges throughout your trajectory


set root "/uufs/chpc.utah.edu/common/home/bidone-group2/robert/funnel_metad/23Sep1/"
set folders "1.bent35A 2.int135A 3.int235A 4.open35A"
set pdbs "boundxtractCAT2.pdb unboundxtractCAT2.pdb allcat.xtc"
mol delete all
package require saltbr
foreach x $folders {
	cd ${root}${x}
	mol delete all
	mol new step1_pdbreader.psf
	file mkdir saltbridges
	cd saltbridges
	foreach y $pdbs {
		animate delete all
		mol addfile ../${y} waitfor all
		set selection1 [atomselect top "residue 0 to 694 or residue 698 to 700"]
		saltbr -sel $selection1 -log saltbridge.log
		$selection1 delete
	}
}
# NOTES:			
# A Valen after the D (Asp) residue (if present) consistently coordinates bonding with ADMIDAS Ca2+ through a water molecule in Beta3.			
# alphaIIb ASP-232 sidechain holds in place two water molecules which H-bonds to Arginine			
# The conformation of the ligand backbone is stabilized as it crosses the interface between the alphaIIb and Beta3 subunits by two previously unremarked water molecules with strong density that are held in place by hydrogen bonds to the side chain of alphaIIb Asp-232 and  backbone of beta3 Ala-218 and that hydrogen bond to the carbonyl oxygen of the Arg of RGD ( Fig. 3 b ).			

# residue 698 to 700
