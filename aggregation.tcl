# Aggregation VMD plugin
#
# Author: Jessica A. Nash
#
# $Id:
#
# Performs analysis of aggregation in a molecular dynamics simulation using
# a simple distance cut off 
# 
# GUI and assisting functions based heavily on VMD hbonds plugin.
# usage: 

package require Tk

namespace eval ::Aggregation:: {
  variable defaultDist 3.5
  variable defaultFrames "all"
  variable defaultWrite 0
  variable defaultDetails "cluster-details.dat"
  variable defaultPDB "agg"
  variable defaultROG 0
  variable defaultBoundaries 0
  set chainCount 0
  
}

proc aggregation_tk {} {
  return [Aggregation::aggregation_gui]
}

# Adapted from hbonds gui which was adapted from pmepot gui
proc ::Aggregation::fill_mol_menu {name} {

  variable usableMolLoaded
  variable currentMol
  variable nullMolString
  $name delete 0 end

  set molList ""
  foreach mm [array names ::vmd_initialize_structure] {
    if { $::vmd_initialize_structure($mm) != 0} {
      lappend molList $mm
      $name add radiobutton -variable Aggregation::currentMol \
        -value $mm -label "$mm [molinfo $mm get name]"
    }
  }

  #set if any non-Graphics molecule is loaded
  if {[lsearch -exact $molList $currentMol] == -1} {
    if {[lsearch -exact $molList [molinfo top]] != -1} {
      set currentMol [molinfo top]
      set usableMolLoaded 1
    } else {
      set currentMol $nullMolString
      set usableMolLoaded  0
    }
  }

}

proc ::Aggregation::getoutdir {} {
  variable guiOutdir

  set newdir [tk_chooseDirectory \
    -title "Choose output directory" \
    -initialdir $guiOutdir -mustexist true]

  if {[string length $newdir] > 0} {
    set guiOutdir $newdir 
  } 
}

proc ::Aggregation::write_state {args} {
  variable w
  variable guiWrite
  variable rogWrite

  # Disable the prefix file field
  if {$guiWrite == 0} {
    if {[winfo exists $w.out.all]} {
      $w.out.all.fbdata configure -state disabled
      $w.out.all.datname configure -state disabled
    }
  } else {
    if {[winfo exists $w.out.all]} {
      $w.out.all.fbdata configure -state normal
      $w.out.all.datname configure -state normal
    }
  }
}

proc ::Aggregation::write_state2 {args} {
  variable w
  variable rogWrite

  # Disable the prefix file field
  if {$rogWrite == 0} {
    if {[winfo exists $w.out.all]} {
      $w.out.all.rgdata configure -state disabled
      $w.out.all.rgdatname configure -state disabled
    }
  } else {
    if {[winfo exists $w.out.all]} {
      $w.out.all.rgdata configure -state normal
      $w.out.all.rgdatname configure -state normal
    }
  }
}

###################################################################
###                            GUI 
###################################################################

proc Aggregation::aggregation_gui {} {
  variable w
  variable defaultDist
  variable defaultWrite
  variable defaultFrames
  variable defaultUpdateSel
  variable defaultDetails
  variable defaultPDB
  variable defaultROG
  variable defaultBoundaries
  
  variable nullMolString "none"
  variable currentMol
  variable molMenuButtonText

 trace add variable Aggregation::currentMol write [namespace code {
    variable currentMol
    variable molMenuButtonText
    if { ! [catch { molinfo $currentMol get name } name ] } {
      set molMenuButtonText "$currentMol: $name"
    } else {
      set molMenuButtonText $currentMol
    }
  # } ]
  set currentMol $nullMolString
  variable usableMolLoaded 0
  
  variable atomselectText1 "all not nucleic"
  
    # Add traces to the checkboxes, so various widgets can be disabled
  # appropriately
  if {[llength [trace info variable [namespace current]::guiWrite]] == 0} {
    trace add variable [namespace current]::guiWrite write ::Aggregation::write_state
  }

  if {[llength [trace info variable [namespace current]::rogWrite]] == 0} {
    trace add variable [namespace current]::rogWrite write ::Aggregation::write_state2
  }

  
  # If already initialized, just turn on
  if { [winfo exists .agg] } {
    wm deiconify $w
    return
  }
  set w [toplevel ".agg"]
  wm title $w "Aggregation Analysis"
  wm resizable $w 0 0

  variable guiDist $defaultDist
  variable guiWrite $defaultWrite
  variable guiFrames $defaultFrames
  variable guiClusterDetailsFile $defaultDetails
  variable guiPDBFile $defaultPDB
  variable guiROG "cluster-rog.dat"
  variable rogWrite $defaultROG
  variable boundaries $defaultBoundaries

  ############## frame for input options #################
  labelframe $w.in -bd 2 -relief ridge -text "Input options" -padx 1m -pady 1m
  
  set f [frame $w.in.all]
  set row 0
  
  grid [label $f.mollable -text "Molecule: "] \
    -row $row -column 0 -sticky e
  grid [menubutton $f.mol -textvar Aggregation::molMenuButtonText \
    -menu $f.mol.menu -relief raised] \
    -row $row -column 1 -columnspan 3 -sticky ew
  menu $f.mol.menu -tearoff no
  incr row
  
  fill_mol_menu $f.mol.menu
  trace add variable ::vmd_initialize_structure write [namespace code "
    fill_mol_menu $f.mol.menu
  # " ]

  grid [label $f.sellabel1 -text "Selection : "] \
    -row $row -column 0 -sticky e
  grid [entry $f.sel1 -width 50 \
    -textvariable Aggregation::atomselectText1] \
    -row $row -column 1 -columnspan 3 -sticky ew
  incr row

  grid [label $f.frameslabel -text "Frames: "] \
    -row $row -column 0 -sticky e
  grid [entry $f.frames -width 10 \
    -textvariable Aggregation::guiFrames] \
    -row $row -column 1 -sticky ew
  grid [label $f.framescomment -text "(all, b:e, or b:s:e)"] \
    -row $row -column 2 -columnspan 2 -sticky w
  incr row

  pack $f -side top -padx 0 -pady 0 -expand 1 -fill none

  set f [frame $w.in.cutoffs]
  set row 0

  grid [label $f.ondistlabel -text "Search cut-off distance (A): "] \
    -row $row -column 0 -sticky e
  grid [entry $f.ondist -width 5 \
    -textvariable Aggregation::guiDist] \
    -row $row -column 1 -columnspan 3 -sticky ew
  incr row


  pack $f -side top -padx 0 -pady 5 -expand 1 -fill x

  pack $w.in -side top -pady 5 -padx 3 -fill x -anchor w

  ############## frame for output options #################
  labelframe $w.out -bd 2 -relief ridge -text "Output options" -padx 1m -pady 1m

  set f [frame $w.out.all]
  set row 0

  grid [label $f.label -text "Output directory: "] \
    -row $row -column 0 -columnspan 1 -sticky e
  grid [entry $f.entry -textvariable Aggregation::guiOutdir \
    -width 35 -relief sunken -justify left -state readonly] \
    -row $row -column 1 -columnspan 1 -sticky e
  grid [button $f.button -text "Choose" -command "::Aggregation::getoutdir"] \
    -row $row -column 2 -columnspan 1 -sticky e
  incr row
  grid [label $f.loglabel -text "File Name: "] \
    -row $row -column 0 -sticky e
  grid [entry $f.logname -width 30 \
    -textvariable Aggregation::guiClusterDetailsFile] \
    -row $row -column 1 -columnspan 2 -sticky ew
  incr row
  
  grid [checkbutton $f.check0 -text \
    "Check across periodic boundaries? (currently assumes origin is box corner)" \
    -variable Aggregation::boundaries] \
    -row $row -column 0 -columnspan 3 -sticky w
  incr row

  grid [checkbutton $f.check2 -text \
    "Write radius of gyration details" \
    -variable Aggregation::rogWrite] \
    -row $row -column 0 -columnspan 3 -sticky w
    incr row
  grid [label $f.rgdata -text "ROG File Name " -state disabled] \
    -row $row -column 0 -sticky e
  grid [entry $f.rgdatname -width 30 \
    -textvariable [namespace current]::guiROG  -state disabled] \
    -row $row -column 1 -columnspan 2 -sticky ew
    incr row
  

  grid [checkbutton $f.check3 -text \
    "Write aggregates as pdb" \
    -variable Aggregation::guiWrite] \
    -row $row -column 0 -columnspan 3 -sticky w
    incr row
  grid [label $f.fbdata -text "PDB Prefix " -state disabled] \
    -row $row -column 0 -sticky e
  grid [entry $f.datname -width 30 \
    -textvariable [namespace current]::guiPDBFile  -state disabled] \
    -row $row -column 1 -columnspan 2 -sticky ew
    incr row

  pack $f -side left -padx 0 -pady 5 -expand 1 -fill x
  pack $w.out -side top -pady 5 -padx 3 -fill x -anchor w

  ############## frame for status #################
  set f [frame $w.control]
  button $f.button -text "Perform Aggregation Analysis" -width 20 \
  -command {::Aggregation::aggregation -sel $::Aggregation::atomselectText1 -frames $::Aggregation::guiFrames -dist $::Aggregation::guiDist -outdir $::Aggregation::guiOutdir -filename $::Aggregation::guiClusterDetailsFile -bound $::Aggregation::boundaries -rogflag $::Aggregation::rogWrite -rogname $::Aggregation::guiROG -pdbflag $::Aggregation::guiWrite -pdbprefix $::Aggregation::guiPDBFile  } 

  pack $f $f.button
  
  return $w
}

proc ::Aggregation::aggregation { args } {
  variable currentMol
  variable endFrame
  variable startFrame
  variable incrNum
  
  set nf [molinfo $currentMol get numframes]
  
  foreach {name val} $args {
    switch -- $name {
      -sel { set arg(sel) $val }
      -frames { set arg(frames) $val }
      -dist { set arg(dist) $val }
      -outdir { set arg(outdir) $val }
      -filename { set arg(filename) $val }
      -bound { set arg(bound) $val }
      -rogflag { set arg(rogflag) $val }
      -rogname { set arg(rogname) $val }
      -pdbflag { set arg(pdbflag) $val }
      -pdbprefix { set arg(pdbprefix) $val }
      default { error "unknown argument: $name $val" }
    }
  }
  
  ## IMPORTANT - Assumes box with 90 degree angles
  set boxDim [pbc get]
  set boxDim [lindex $boxDim 0]
  set xDim [lindex $boxDim 0]
  set yDim [lindex $boxDim 1]
  set zDim [lindex $boxDim 2]

  set xLim [expr $xDim/2]
  set yLim [expr $yDim/2]
  set zLim [expr $zDim/2]
  
  
   ## Set variables from user input
  set sel [atomselect top $arg(sel)]
  set fragments [$sel get fragment]
  
  ## Open files for writing data
  set outfile [open $arg(filename) "w"]
  if {$arg(rogflag)==1} {
    set rogFile [open $arg(rogname) "w"]
  }
  
  if {$arg(frames)=="all"} {
    set incrNum 1
    set startFrame 0
    set endFrame $nf
  } else {
    set fl [split $arg(frames) :]
    if {[llength $fl]==3} {
      set startFrame [lindex $fl 0]
      set incrNum [lindex $fl 1]
      set endFrame [lindex $fl 2]
    }
    
    if {[llength $fl]==2} {
      set startFrame [lindex $fl 0]
      set endFrame [lindex $fl 1]
      set incrNum 1
    }

  }
  
  puts $outfile "Frame\tNumber of Aggregates\tAverage Aggregate Size"
  
  ##Run aggregation analysis
  
  for {set x $startFrame} {$x < $endFrame} {incr x $incrNum} {
	set aggList []
	set aggListAll []
	set aggListAll2 []
	set rogAll []
	set lengthAgg []
		
	set aggAvg 0
	set aggCount 0
	set fragments_unique [lsort -integer -increasing -unique $fragments]
	set aggNumPrevious 0
	
		#Look for clusters - Initial Pass
		set Ulength [llength $fragments_unique]

		for {set y 0} {$y<$Ulength} {incr y} {
		   #Look for the fragment # of interest (y) in the list of fragments
		   #already assigned to an aggregate
		   set foundN [lsearch $aggListAll [lindex $fragments_unique $y]]
		   
		   #If the fragment of interest is not in the list, perform 
		   #aggregation analysis
		   if {$foundN==-1} {
			   #Set some variables to start. aggNum gives the number of 
			   #fragments per aggregate
			   #fragment P gives the initial fragment of interest. This 
			   #variable will be appended to as analysis continues
			   set aggNum 1
			   set fragmentP [lindex $fragments_unique $y]
			   
			   #Place holder. initial fragment
			   set original $fragmentP
			   
			   #Analysis flags
			   set doneMoving 0
			   set moveFlag 0
			   set chainCount 0
			   set dim 0
			   set moveTwice 0
			   
			   while {$doneMoving==0} {
				   puts $dim
				   #Look for fragments within specified distance of specified fragments
				   set selString [concat $arg(sel) "and within " $arg(dist) "of ( " $arg(sel) " and fragment " $fragmentP ")"]
				   set aggSel [atomselect top $selString frame $x]
				   
				   set aggRes [$aggSel get fragment]
				   set resSel [atomselect top "fragment $aggRes"]
				   
				   set aggResid [$resSel get resid]
				   
				   #Get unique fragments
				   set aggRes_unique [lsort -unique $aggRes]
				   set aggResid_unique [lsort -unique $aggResid]
				   set aggNumPrevious $aggNum
				   set aggNum [llength $aggRes_unique]
				   set fragmentP $aggRes_unique
				   set aggList [concat $aggList $fragmentP]
				   
				   if {$aggNumPrevious==$aggNum  && $arg(bound)==1} {
					   # If search appears to be finished, need to check
					   # periodic boundaries.
					   set moveSel [atomselect top "fragment $fragmentP" frame $x]
					   array set coords {}
					   set coords(0) [$aggSel get {x}]
					   set coords(1) [$aggSel get {y}]
					   set coords(2) [$aggSel get {z}]
					   
					   set coordsDim $coords($dim)
					   set coordSorted [lsort -increasing $coordsDim]
					   set coorFirst [lindex $coordSorted 0]
					   set coorLast [lindex $coordSorted end]
					   set diffFirst [expr abs($coorFirst-0)]
					   set diffLast [expr abs($coorLast-[lindex $boxDim $dim])]
					   set moveBy [list 0 0 0]
					   set moveNum [lindex $boxDim $dim]
					   
					   if {$diffFirst< [expr 2*$arg(dist)] || $coorFirst<0} {
						   set moveBy [lreplace $moveBy $dim $dim $moveNum]
						   $moveSel moveby $moveBy
						}
						   
					   if {$diffLast< [expr 2*$arg(dist)] || $coorLast>[lindex $boxDim $dim]} {
						   set moveBy [lreplace $moveBy $dim $dim -$moveNum]
						   $moveSel moveby $moveBy
						}
						
					   $moveSel delete
					   
					   if {$dim==2} {
						   set dim -1
						   if {$moveTwice==1} {
							   set doneMoving 1
							}
						   set moveTwice 1
						}
						
						set dim [expr $dim +1 ]
					}

					
					#Condition for non-periodic boundaries
					if {$aggNumPrevious==$aggNum && $arg(bound)==0} {
						set doneMoving 1
					}
					
					$aggSel delete
			
				}


			  if ($x==[expr $nf-$incrNum]) {
				  mol selection $selString
				  set colorNum $aggCount
				  if $colorNum>32 {
				   set colorNum [expr $colorNum-32]
				  }
				  
				  #Color aggregates
				  mol color ColorID $colorNum
				  mol representation Licorice 0.3 10.0 10.0
				  mol addrep top
				}
			
			  
			  set selChainn [atomselect top $selString frame $x]
			  set selSASA [atomselect top $selString frame $x]
		   
			
			  set rog2 [measure rgyr $selChainn]
			  set nLP [llength $fragmentP]
			
				if {[llength $fragmentP]>0} {
				   set rogAll [concat $rogAll $rog2]
				   set lengthAgg [concat $lengthAgg [llength $fragmentP]]
				   if {$arg(pdbflag)==1} { 
					   set pdbString $arg(pdbprefix)
					   append pdbString "_n"
					   append pdbString $nLP
					   append pdbString "_"
					   append pdbString $x
					   append pdbString "_"
					   append pdbString [lindex $fragmentP 1]
					   append pdbString ".pdb"
					   set pdbString [string trim $pdbString " "]
					   $selChainn writepdb $pdbString
					}

				}
			
		  
			  set chainCount [expr $chainCount+1]
			  $selChainn delete



			  set aggListAll [concat $aggListAll $fragmentP]
			  lappend aggListAll2 $aggRes_unique
			  set aggCount [expr $aggCount+1]
			  set aggAvg [expr $aggAvg+$aggNum]

			  if {$aggNumPrevious==1} {
				  set aggNumPrevious 0
			  }
			   
			   
			   #Clears list for next aggregate.
			   set aggList []

			}
		}
	   
		set aggAvgP $aggAvg
		set aggAvg [expr double($aggAvg)/double($aggCount)]
		if {$arg(rogflag)==1} {
		puts $rogFile "$x $lengthAgg $rogAll"
		}
	
		puts $outfile "$x\t$aggListAll2"



	}

close $outfile

if {$arg(rogflag)==1} {
  close $rogFile
}

  
  
}
  

