# Metadynamics_extractframes.m is to be used in matlab first.
# If you are going to use the produced .pdb frames in gromacs they a little editing to make them compatible.
# these commands will do the trick in a linux shell (bash). Change the file names for your use. e.g. unboundxtractCAT.pdb
# cat unboundxtractCAT.pdb | grep -v CRYST1 > temp.pdb #(grep -v means to invert the match)
# mv -f temp.pdb unboundxtractCAT.pdb
# sed -i 's/END/ENDMDL/g' unboundxtractCAT.pdb

set inputfile "bound"
set folders "1.bent35A/ 2.int135A/ 3.int235A/ 4.open35A/"
#set folders "1.bent35A/"
set root "C:/Users/and_r/OneDriveUtah/BidoneLab/integrin/funnel_metad/23Sep1/"
proc extract {inputfile folders subfolders} {
global root;
mol delete all
	foreach x $folders {
		cd ${root}${x}
		mol new step1_pdbreader.psf
		foreach y $subfolders {
			if { [catch {cd ${root}${x}${y}} ] } {
				puts "couldn't find folder ${x}${y}"
				continue
			}
			animate delete all
			set existfile [catch {file size ${inputfile}.txt} size]
			if {$existfile != 0 || $size < 3} {
				puts "${inputfile} file is empty in subfolder $y of $x"
				continue
			} else {
				set fp [open "${inputfile}.txt" r]
				set file_data [read $fp]
				close $fp
				set endframe [lindex $file_data end]
				puts $endframe
				
				mol addfile trajout.xtc first 0 last $endframe waitfor all
				set all [atomselect top all]
	
				file mkdir selected_structures/$inputfile/
				foreach z $file_data {
					$all frame $z
					$all writepdb selected_structures/${inputfile}/${z}.pdb
				}
				animate delete all				
				foreach f [glob selected_structures/${inputfile}/*.pdb] {
					mol addfile $f waitfor all
				}
				puts "finished $y"
				animate write dcd selected_structures/${inputfile}/${inputfile}extractedframes2.dcd waitfor all
			}
		}
	
		animate delete all
		foreach a $subfolders {
			if { [catch {cd ${root}${x}${a}} ] } {
				#puts "couldn't find folder ${x}${a}"
				continue
			}
			set existfile2 [catch {file size ./selected_structures/${inputfile}/${inputfile}extractedframes.dcd} size]
			if {$existfile2 != 0 || $size < 3} {
				set pwd [pwd]
				#puts "$y $pwd ${inputfile}extractedframes.dcd is empty or doesn't exist"
				continue
			} else {
				mol addfile ./selected_structures/${inputfile}/${inputfile}extractedframes2.dcd waitfor all
			}
		}
		animate write dcd ${root}${x}${inputfile}xtractCAT2.dcd waitfor all
		animate write pdb ${root}${x}${inputfile}xtractCAT2.pdb waitfor all
		puts "finished writing dcd and pdb of extracted frames"
	}
}

extract $inputfile $folders "1 2 3 4 5 6 7 8" 
mol delete all
foreach k $folders {
	cd ${root}${k}
	mol new step1_pdbreader.psf
	mol addfile ${inputfile}xtractCAT.dcd waitfor all
	mol rename top $k
}
return
