#
# for measuring the distance between atoms in the open or closed conformation and puts all gathered data into one file. These are the indices determined by ref.pdb
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
set folders "1.bent35A 2.int135A 3.int235A 4.open35A"
#set pdbs "boundxtractCAT.pdb unboundxtractCAT.pdb"
mol delete all
package require timeline
foreach y $folders {
	cd ${root}${y}
	mol new step1_pdbreader.psf
	for {set x 1} {$x <= 16} {incr x} {
		if { [file exists boundpose${x}.pdb ] == 0 } {
			puts "couldn't find file boundpose${x}.pdb"
			continue
		}
		set existfile [catch {file size boundpose${x}.pdb} size]
		if {$existfile != 0 || $size < 3} {
			puts "boundpose${x}.pdb file is empty in folder $y"
			continue
		} else {
			animate delete all
			set fileid [open pose${x}secondarystructure.txt w]
			mol addfile boundpose${x}.pdb waitfor all
			set selection1 [atomselect top "residue 0 to 694 and alpha"]
			foreach resid [$selection1 get resid] segname [$selection1 get segname] {
				set selection [atomselect top "resid $resid and alpha and segname $segname"]
				set structure [$selection get structure]
				puts $fileid "$resid P $segname 0 $structure"
				#puts "$resid P $segname 0 $structure"
				$selection delete
			}
			$selection1 delete
		}
		close $fileid
	}
}
return
# NOTES:			
# A Valen after the D (Asp) residue (if present) consistently coordinates bonding with ADMIDAS Ca2+ through a water molecule in Beta3.			
# alphaIIb ASP-232 sidechain holds in place two water molecules which H-bonds to Arginine			
# The conformation of the ligand backbone is stabilized as it crosses the interface between the alphaIIb and Beta3 subunits by two previously unremarked water molecules with strong density that are held in place by hydrogen bonds to the side chain of alphaIIb Asp-232 and  backbone of beta3 Ala-218 and that hydrogen bond to the carbonyl oxygen of the Arg of RGD ( Fig. 3 b ).			

# residue 698 to 700

mol new titin.psf type psf first 0 last -1 step 1 filebonds 1 \
        autobonds 1 waitfor all
mol addfile titin.dcd type dcd first 0 last -1 step 1 filebonds 1 \
        autobonds 1 waitfor all
set myMolid 0

set outFilename titin-per-res-rmsd.tml

proc myRmsdBatchCalc {filename molid} {
Since we are making a per-residue data file, usesFreeSelection is set to 0. (See the follwinng subsection for more on "Free Selection".)
  # 0 for per-residue calcs
  # 1 for per-selection calcs
  set usesFreeSelection 0
Set all the following values, which will satisfy all user-set header data Here, we set set the chain and the number of residues we examine manually. This molecule has only 1 chain; for molecules with more chains we would have to generate data for all chains.
  # include units in title
  set dataTitle "res. RMSD (A)"
  set firstFrame 0
  set lastFrame 96
  # the sample molecule is a 1-chain, 1-segment molecule, so 
  # we will do a simple loop over residue numbers
  set theChain "T"
  set theSeg "TIT"
  set firstRes 1
  set lastRes 89

  set numFrames [expr $lastFrame - $firstFrame + 1]
Set the number of items (referred to as selection groups)
  set numSelectionGroups [expr $lastRes-$firstRes + 1]
Check the filename and open the file for writing
if {$filename == ""  } {
    die "usage: myRMSDBatchCalc FILENAME MOLID\n FILENAME\
           cannot be an empty string."
   }

set outDataFile [open $filename w]
Now write out the .tml data file header.
::timeline::writeDataFileHeader $outDataFile $molid $dataTitle \ 
     $numFrames $numSelectionGroups $usesFreeSelection
Now perform the actual analysis calculations for the whole trajectory, looping over residues and frames, and output the data. For this example, we calculate the RMSD of each individual residue for each frame of the simulation vs. its initial (frame 0) configuration
  set chain $theChain
  set seg $theSeg
  for {set r $firstRes} {$r <= $lastRes} {incr r} {
    set sela [atomselect $molid "resid $r"]
    set selb [atomselect $molid "resid $r"]
    $sela frame 0
    for {set f $firstFrame} {$f<=$lastFrame} {incr f} {
      $selb frame $f
      display update
      set val [measure rmsd $sela $selb]
      set resid $r
      puts $outDataFile "$resid $chain $seg $f $val"
    }
  }
Now that all the data has been written, close the data file and end the procedure.
  close $outDataFile
  return
}
Now we call the procedure we just entered.
myRmsdBatchCalc $outFilename $myMolid
Generate the .tml file by running the .tcl script in VMD, in either an interactive session using the source command or at a command prompt:
vmd -dispdev text -eofexit < titin-per-res-rmsd.tcl > titin-per-res-rmsd.log .