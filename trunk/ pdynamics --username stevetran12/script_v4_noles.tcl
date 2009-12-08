# pDynamics
#	build water box, fixed atom files, ionized psf, pdb files				
# 	you should run this script in terminal or tkconsole with	 	
#	command like:				 	
#	vmd -dispdev text -e script_v3.tcl -args <pdbfile> <outputname> <topology file> > output.log
#   							

puts "Number of arguments $argc"
if { $argc != 4 } {
   puts "USAGE:   vmd -dispdev text < vmd_select.tcl -args <pdbfile> <outputname> <topology file>"
   exit
}

	#original pdb
	set pdbfile [lindex $argv 0]

	#output name
	set outputname [lindex $argv 1]

	#topology file
	set topfile [lindex $argv 2]

	#mutate switch
	set mutate off

	#what to mutate (resid)
	set mutid 64

	#what to mutate it to
	set mutres trp

	#how many oxygens to make
	set oxylim 15

#split segment chains

cat $pdbfile | grep ATOM > $outputname\_pro.pdb
cat $pdbfile | grep HETATM | grep HEM > $outputname\_hem.pdb
cat $pdbfile | grep HETATM | grep OXY > $outputname\_oxy.pdb

#load package plugins to load psf builder, solvation plugin (waterbox), and ionize water box
package require psfgen
package require solvate
package require autoionize

# clear structures and reset psf builder
mol delete all
resetpsf

#correcct atom names in pdb and code
pdbalias atom TYR:151 O OT1
pdbalias atom TYR:151 OXT OT2
pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD
pdbalias residue HEM HEME
pdbalias residue OXY O2

# mol load pdb input
mol load pdb $outputname\_pro.pdb
mol load pdb $outputname\_hem.pdb
mol load pdb $outputname\_oxy.pdb

#load topology file
topology $topfile

#create protein segment
segment PRO {
  auto angles dihedrals 
  pdb $outputname\_pro.pdb
	if { $mutate } {
	mutate $mutid $mutres
	}
}

#create heme segments
segment HEM {
  auto  angles dihedrals 
  pdb $outputname\_hem.pdb
}

#create oxygens coordiante pdbs
for {set x 1} {$x<=$oxylim} {incr x} {
set otwo "O"
append otwo $x
 segment $otwo {
  auto  angles dihedrals 
  pdb $outputname\_oxy.pdb
 }
#}

# Load the coordinates for each segment.
coordpdb $outputname\_pro.pdb PRO
coordpdb $outputname\_hem.pdb HEM
for {set x 1} {$x<=$oxylim} {incr x} {
set otwo "O"
append otwo $x
coordpdb $outputname\_oxy.pdb $otwo
}



#Patches
  patch PHEM PRO:93 HEM:154
  regenerate angles
  patch FHEM HEM:154

# Write out the psf file
writepsf $outputname\_noles.psf

# Guess the positions of missing atoms.  As long as all the heavy
# atoms are present, psfgen usually does a very good job of this.
guesscoord
writepdb $outputname\_noles.pdb

mol delete all
resetpsf
mol load pdb $outputname\_noles.pdb psf $outputname\_noles.psf
set pro_minmax [measure minmax [atomselect top "protein or segname HEM"]]
puts "$pro_minmax"
#load pdbs again and select and split segments for solvation

foreach { pro_min pro_max } $pro_minmax {break}
foreach { pro_minx pro_miny pro_minz } $pro_min {break}
foreach { pro_maxx pro_maxy pro_maxz } $pro_max {break}
puts "protein min xyz: $pro_minx $pro_miny $pro_minz"
puts "protein max xyz: $pro_maxx $pro_maxy $pro_maxz"

set coordx [expr $pro_minx - [expr rand() * 10 ]]

puts "random x coord $coordx"

for {set x 2} {$x<=$oxylim} { set x $x} {

#splitting for each segment specificed an wrapping from vmd scripting tutorail from cincinatti
	#xmin
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec [list  [expr $pro_minx - 3 - [expr rand() * 7 ]] [expr { ($pro_maxy - $pro_miny) * rand() + $pro_miny } ] [expr { ($pro_maxz - $pro_minz) * rand() + $pro_minz } ] ]
		$sel moveby [vecsub  $newvec $coord]
		puts "[vecsub  $newvec $coord]]"
	set x [expr $x+1]

	#xmax
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec [list [expr $pro_maxx + 3 + [expr rand() * 7 ]] [expr { ($pro_maxy - $pro_miny) * rand() + $pro_miny } ] [expr { ($pro_maxz - $pro_minz) * rand() + $pro_minz } ]]
		$sel moveby  [vecsub  $newvec $coord]
		puts "newcoord: [vecsub  $newvec $coord]"
	set x [expr $x+1]

	#ymin
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec  [list [expr { ($pro_maxx - $pro_minx) * rand() + $pro_minx } ] [expr $pro_miny - 3 - [expr rand() * 7 ]] [expr { ($pro_maxz - $pro_minz) * rand() + $pro_minz } ]]
		$sel moveby  [vecsub  $newvec $coord]
		puts "newcoord: [vecsub  $newvec $coord]"
	set x [expr $x+1]

	#ymax
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec [list [expr { ($pro_maxx - $pro_minx) * rand() + $pro_minx } ]  [expr $pro_maxy  + 3 + [expr rand() * 7 ]] [expr { ($pro_maxz - $pro_minz) * rand() + $pro_minz } ]]
		$sel moveby  [vecsub  $newvec $coord]
		puts "newcoord: [vecsub  $newvec $coord]"
	set x [expr $x+1]

	#zmin
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec [list [expr { ($pro_maxx - $pro_minx) * rand() + $pro_minx } ] [expr { ($pro_maxy - $pro_miny) * rand() + $pro_miny } ]   [expr $pro_minz - 3 - [expr rand() * 7 ]] ]
		$sel moveby  [vecsub  $newvec $coord]
		puts "newcoord: [vecsub  $newvec $coord]"
	set x [expr $x+1]

	#zmax
	set otwo "O"
	append otwo $x
	set sel [atomselect top "segname $otwo"]
	foreach coord [[atomselect top "segname $otwo and type O2"] get {x y z}] {break}
	set newvec [list [expr { ($pro_maxx - $pro_minx) * rand() + $pro_minx } ] [expr { ($pro_maxy - $pro_miny) * rand() + $pro_miny } ]   [expr $pro_maxz + 3 + [expr rand() * 7 ]] ]
		$sel moveby  [vecsub  $newvec $coord]
		puts "newcoord: [vecsub  $newvec $coord]"
	set x [expr $x+1]

}
#guess coordinates using psfgen and write out pdb
guesscoord
[atomselect top "all"] writepdb $outputname\_noles.pdb

#reset and Add water
mol delete all

resetpsf



# Read in structure files
readpsf $outputname\_noles.psf 

# Read in coordinates
coordpdb $outputname\_noles.pdb


resetpsf

# solvation package (solvate) for water box
solvate  $outputname\_noles.psf  $outputname\_noles.pdb -o $outputname\_md_wb -t 1
autoionize -psf  $outputname\_md_wb.psf -pdb $outputname\_md_wb.pdb -is 0.1 -o $outputname\_md_ion

#clear molecules and load structures for ionization in water box
resetpsf
mol delete all
mol load pdb $outputname\_md_ion.pdb psf  $outputname\_md_ion.psf

# Read in structure files
readpsf $outputname\_md_ion.psf
# Read in coordinates
coordpdb $outputname\_md_ion.pdb
if { $opatch } {
set file [open $outputname\_noles_box_size.dat "w"]
  set center [measure center [atomselect top "all"]]
  set sides  [vecinvert [eval vecsub [measure minmax [atomselect top "water and noh"]]]]
  #puts $file "BOX SIZE:"
  puts $file "cellBasisVector1 [lindex $sides 0]\t0\t\t0"
  puts $file "cellBasisVector2 0\t\t[lindex $sides 1]\t0"
  puts $file "cellBasisVector3 0\t\t0\t\t[lindex $sides 2]"
  puts $file "cellOrigin       [lindex $center 0]\t[lindex $center 1]\t[lindex $center 2]"
  close $file 
}

#mutate switch and what to mutate to 
regenerate angles

	if { $opatch } {
		puts "patch on"
		patch PLO2 O1:154 HEM:154
	}


#write out structures inside water box
writepsf $outputname\_noles_box.psf
[atomselect top all] writepdb $outputname\_noles_box.pdb
mol delete all

#creating fixed atom files on segments 
#fixed protein structures, backbone, and alha structures
mol load pdb $outputname\_noles_box.pdb
set basename $outputname\_noles_box
[atomselect top all] set beta 0
  [atomselect top "protein or segname HEM O1"] set beta 1
  [atomselect top all] writepdb $basename-fixed-protein-heme-o2.pdb
  
  [atomselect top all] set beta 0
  [atomselect top "backbone or segname HEM O1"] set beta 1
  [atomselect top all] writepdb $basename-fixed-backbone-heme-o2.pdb
  
[atomselect top all] set beta 0
  [atomselect top "name CA"] set beta 1
  [atomselect top all] writepdb $basename-restrain-CA.pdb

exit
