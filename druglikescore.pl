#!/usr/bin/perl -w 

# compute the druglikeliness score of the SMILES given as input
#
# Reference:
# M.C. Hutter "Separating Drugs from Nondrugs: A Statistical Approach
# Using Atom Pair Distributions" J.Chem.Inf.Model. 47 (2007) 186-194.
#
# usage: druglikescore.pl  infile.smi
#
# output: compoundname <white space>  druglikeliness score
#
# Revision history
# 28/07/2004 now recognizes spaces or tabs between SMILES and compound name

$INFILE = $ARGV[0] ;     

$dative = 1 ;
$neutra = 1 ;  # neutralize acidic oxygens
$aromat = 1 ;  # force aromatic bond information in .hin file
#$hyperc = 1 ;  # HYPERCHEM compatible atom type assignment
               # -NO2 groups:  N1(=O1)(=O1) instead of NO(=O1)(=O1)
               # -N3 groups:   NA=N1=NA     instead of NA=NZ=N1 
# MM force field-like atom types used! ($hyperc = 0)

open(INPU,$INFILE) or die "error opening $INFILE : $!\n"; 

$hinname = " " ;
$natoms = 0           ;  # total number of atoms
$nheavy = 0           ;  # number of non-hydrogen atoms
$nhydr = 0            ;  # number of hydrogen atoms

while(defined($i = <INPU>)) {
# @cols = split(/\s+/,$i);
  @cols = split(/[\s\t]+/,$i); # split at spaces or tabs

  $smi = " "  ;     # compound SMILES
  $cid = "" ;       # compound CID   will be file name of the .hin file

  if ($cols[0] ne "") {
    $smi = $cols[0] ;
  }

  if (defined($cols[1])) {
    if ($cols[1] ne "") {
      $cid = $cols[1] ;
    }
  }

  if (defined($cols[2])) {
    if ($cols[2] ne "") {      # splitted name, e.g. acetic ester
      $cid = $cid.$cols[2] ;
    }
  }

# identifier or name of compound
# if there is no second coloum in the .smi file, its name is used instead
  $lll = length($cid) ;
  if ($lll < 1) {
    $cid = substr($INFILE,0,length($INFILE)-4) ;
  }
  $hinname = $cid ;

# print "$hinname " ;

  if ($smi ne "") { 

# substitute [NH] to N
    $smi =~ s/\[NH\]/N/g ;
#   print "$smi\n" ;

#   now convert SMILES to .hin
#   interprete SMILES: get each character build up connection matrix
    $natoms = 0      ;  # total number of atoms
    $nheavy = 0      ;  # number of non-hydrogen atoms
    $nhydr = 0       ;  # number of hydrogen atoms
    $single = "s"    ;
    $double = "d"    ;
    $triple = "t"    ;
    $aromatic = "a"  ;
    $pos = 0         ;  # current position / atom  in SMILES string
    $bat = 0         ;  # last atom position (also from where branching occurs)
    $cm1 = " "       ; 

# init arrays
    for ($i = 0; $i <= length($smi) ; $i++) {
       $neighbors[$i] = "" ;
       $el[$i] = "" ;    # atom position in SMILES string
       $elem[$i] = "" ;  # atoms numbered consequtively
       $emap[$i] = -1 ;  # map conseq. atoms position to pos. in SMILES
       $chrg[$i] = 0 ;
       $nsub[$i] = 0 ;   # number of substituents of each atom
       $bondlist[$i] = "" ;
    }
# number of ring systems that are open at the same time
# is limited here to 9 !
    for ($i = 0; $i <= 9  ; $i++) { 
       $ring[$i] = -1 ;
    }

    for ($j=1 ; $j <= length($smi) ; $j++) {
      $c1 = substr($smi,$j-1,1) ;   # one character (e.g. C, N, O)
      if ($j < length($smi)) {
        $c2 = substr($smi,$j,1) ;
      }
      else { 
        $c2 = "" ;                  # next character (e.g. Cl, Br, Si)
      }


# atom or branching or bond or ring opening/closing ?
# element
# still to do: treat stereoinformation @ and @@
      if ($c1 eq "@") {           # skip @ and @@
        $pos++ ;
        $cm2 = $cm1 ; # pre-preceeding character
        $cm1 = $c1 ; # preceeding character
        next ;
      }
      if ($c1 eq "H" && $cm1 eq "@") {
        $pos++ ;
        $cm2 = $cm1 ; # pre-preceeding character
        $cm1 = $c1 ; # preceeding character
        next ;
      }

# still to do: treat cis/trans information / and \


      if ($c1 =~ /[a-z]/) {next}  # skip lower case letters  
      if ($c1 eq "[" && $pos == 0) {next}  # if SMILES starts with [
      if ($c1 =~ /[A-Z]/) {       # element (always starts with upper case letter)
        $nheavy++;
        if ($c2 =~ /[a-z]/) {
          $el[$pos] = $c1.$c2 ;
        }
        else {
          $el[$pos] = $c1 ;
        }
#       $neighbors[$pos] = "$pos" ;  # obsolte !!!
        if ($pos > 0) {  # bonded to previous atom
          $neighbors[$bat] = $neighbors[$bat]." ".$pos ; 
          $neighbors[$pos] = $neighbors[$pos]." ".$bat ; 
          $nsub[$pos]++ ;
          $nsub[$bat]++ ;
          if ($c2 =~ /\+/) { $chrg[$pos] =  1 }
          if ($c2 =~ /\-/) { $chrg[$pos] = -1 }
          if ($cm1 eq "[") {$cm1 = $cm2} # cases like C=[N+] 
          if ($cm1 eq "=") {
            $bondlist[$bat] = $bondlist[$bat]." ".$double ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$double ; 
          }
          elsif ($cm1 eq "#") {
            $bondlist[$bat] = $bondlist[$bat]." ".$triple ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$triple ; 
          }
          elsif ($cm1 eq ":") {
            $bondlist[$bat] = $bondlist[$bat]." ".$aromatic ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$aromatic ; 
          }
          else {
            $bondlist[$bat] = $bondlist[$bat]." ".$single ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$single ; 
          }
        }
        $bat = $pos ;
#       if ($c2 =~ /\+/) { $chrg[$nheavy] =  1 }
#       if ($c2 =~ /\-/) { $chrg[$nheavy] = -1 }
        if ($c2 =~ /[a-z]/) { $pos++ }
      }

# open branch
# the branching atom is put on the stack
      if ($c1 =~ /[(]/) {
        push(@blist,$bat) ;
##      print "$c1 $bat\n" ;
      }

# close branch
# the latest branching atom is retrieved from the stack
      if ($c1 =~ /[)]/) {
        $bat = pop(@blist) ;
##      print "$c1 $bat\n" ;
      }

# open/close ring
      if ($c1 =~ /[1-9]/) {
        if ($ring[$c1] == -1) {  # ring opening
          $ring[$c1] = $bat ;
          $nsub[$bat]++ ;
#         print " open ring system $c1 at position $bat\n" ;
        }
        else {                   # ring closure, connect atoms
          $neighbors[$bat] = $neighbors[$bat]." ".$ring[$c1] ; 
          $neighbors[$ring[$c1]] = $neighbors[$ring[$c1]]." ".$bat ; 
          $nsub[$bat]++ ;
#         print "close ring system $c1 between $bat and $ring[$c1]\n" ;
##        print "bat $bat  ring $ring[$c1]\n" ;
          if ($cm1 eq "[") {$cm1 = $cm2} # cases like C=[N+] 
          if ($cm1 eq "=") {
            $bondlist[$bat] = $bondlist[$bat]." ".$double ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$double ; 
          }
          elsif ($cm1 eq "#") {
            $bondlist[$bat] = $bondlist[$bat]." ".$triple ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$triple ; 
          }
          elsif ($cm1 eq ":") {
            $bondlist[$bat] = $bondlist[$bat]." ".$aromatic ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$aromatic ; 
          }
          else {
            $bondlist[$bat] = $bondlist[$bat]." ".$single ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$single ; 
          }

          $ring[$c1] = -1 ; # reset to allow re-use of this ring number
        }
      }

# explicit bonds
#     if ($c1 =~ /[=#:]/) {
##       print "explicit bond $c1\n" ;
#     }

      $cm2 = $cm1 ; # pre-preceeding character
      $cm1 = $c1 ; # preceeding character
      $pos++ ;

   }

# now derive atoms sequentially 
# setup connectivity matrix $neib[][], add hydrogens
# $elem[] $typ[] $charge[] $typ[] $nsubst[] $neib[][] $bond[][]

for ($i = 0; $i < $pos ; $i++) {
  if ($el[$i] =~ /[A-Z]/) {
    $natoms++ ; 
    if ($el[$i] ne "H") {
      $nheavy++ ;
    }
    else {
      $nhydr++ ;
    }
    $emap[$i] = $natoms ;    # atom position in SMILES to conseq. number
    $elem[$natoms] = $el[$i] ;
    $typ[$natoms] = "**" ;
    $charge[$natoms] = $chrg[$i] ;
    $nsubst[$natoms] = $nsub[$i] ;
#   split bond and neighbor lists
    @coln = split(/\s+/,$neighbors[$i]) ;
    @colb = split(/\s+/,$bondlist[$i]) ;
    for ($j = 1; $j <= $nsubst[$natoms] ; $j++) {
##     print "$j $coln[$j] $colb[$j]\n" ;
      $neib[$natoms][$j] = $coln[$j] ;
      $bond[$natoms][$j] = $colb[$j] ;
    }
  }
}

# now correct atom positions in neighbor list 
for ($i=1 ; $i <= $natoms ; $i++) {
  for ($j = 1; $j <= $nsubst[$i] ; $j++) {
    $neib[$i][$j] = $emap[$neib[$i][$j]] ;
  }
}


# Determine valencies and add hydrogens
$newato = $natoms ;
for ($i=1 ; $i <= $natoms ; $i++) {
  $neb = 0 ;  # number of bonds 
  $val = 0 ;  # number of valencies
  if ($nsubst[$i] > 0) {
    for ($k=1 ; $k <= $nsubst[$i] ; $k++) {
      if ($bond[$i][$k] eq "s" ) {
        $neb = $neb + 1 ;
      }
      elsif ($bond[$i][$k] eq "d") {
        $neb = $neb + 2 ;
      }
      elsif ($bond[$i][$k] eq "t") {
        $neb = $neb + 3 ;
      }
      else {   # aromatic
        $neb = $neb + 1 ;
      }
    }
  }
  if ($elem[$i] eq "B" ) {
    if ($charge[$i] == -1) {
      $val = 4 ;
      $typ[$i] = "B4" ;
    }
    else {
      $val = 3 ;
      $typ[$i] = "B3" ;
    }
  }
  if ($elem[$i] eq "C" ) {
    if ($charge[$i] == -1) {
      $val = 3 ;
    }
    else {
      $val = 4 ;
    }
  }
  if ($elem[$i] eq "Si") {
    $val = 4 ;
    $typ[$i] = "SI" ;
  }
  if ($elem[$i] eq "O" ) {
    if ($charge[$i] == -1) {
      $val = 1 ;
#     $val = 0 ;
    }
    elsif ($charge[$i] == 1) {
      $val = 3 ;
    }
    else {
      $val = 2 ;
    }
    if ($neb == 1 && $dative > 0) {  # no hydrogens added to [N+][O-] 
      $nn = $neib[$i][1] ;
      if ($elem[$nn] eq "N" && $nsubst[$nn] == 3) {
        $val = 1 ;
        $typ[$i] = "O1" ;
#       print "atom $i $elem[$i] $typ[$i] $charge[$i] $nsubst[$i]\n" ;
      }
    }
    if ($charge[$i] == -1 && $neutra > 0 && $typ[$i] eq "**") {
      $nn = $elem[$neib[$i][1]] ;
      if ($nn eq "C" || $nn eq "S" || $nn eq "P" ) {
        $charge[$i] = 0 ;   # neutralize acidic oxygens
        $val = 2 ;
      }
    }
  }
  if ($elem[$i] eq "F" ) {
    $val = 1 ;
    $typ[$i] = "F" ;
  }
  if ($elem[$i] eq "Cl") {
    $val = 1 ;
    $typ[$i] = "CL" ;
  }
  if ($elem[$i] eq "Br") {
    $val = 1 ;
    $typ[$i] = "BR" ;
  }
  if ($elem[$i] eq "I" ) {
    $val = 1 ;
    $typ[$i] = "I" ;
  }

#   Nitrogen, Sulfur and Phoshorus are treated individually
  if ($elem[$i] eq "N" ) {
    if ($neb > 3) {
      $val = $neb ;
      $charge[$i] = 1 ;
    }
#   if ($neb == 2) {      # wouldn't add hydogens to -N-
#     $charge[$i] = -1 ;
#     $val = 2 ;
#   }
    if ($charge[$i] == 0) {
      $val = 3 ;
    }
    if ($dative > 0) {   # convert dative to directional bonds
# dative nitro groups to double bonds   [N+]([O-])=O  to N(=O)=O
      $oneib   = 0 ; # number of oxygens as substituents
      $onei[1] = 0 ; 
      $onei[2] = 0 ; 
      $onei[3] = 0 ; 
      $nneib   = 0 ; # number of nitrogens as substituents
      $nnei[1] = 0 ;
      $nnei[2] = 0 ;
      $nnei[3] = 0 ;
      for ($k=1 ; $k <= $nsubst[$i] ; $k++) {
        $nn = $neib[$i][$k] ;
        if ($elem[$nn] eq "O") {
          $oneib++ ;
          $onei[$oneib] = $nn ;
        }
        if ($elem[$nn] eq "N") {
          $nneib++ ;
          $nnei[$oneib] = $nn ;
        }
      }
      if ($nsubst[$i] == 3) {
        if ($oneib == 2 || $oneib == 3) {  # NO2 
          $charge[$i] = 1 ;
#         if ($hyperc == 1) {
#           $typ[$i] = "N1" ;   # Hyperchem assigns N1 to nitro nitrogens
#         }
#         else {
            $typ[$i] = "NO" ;   # MM3 assigns NO
#         }
          for ($k=1 ; $k <= $oneib ; $k++) {
            $typ[$onei[$k]] = "O1" ;
#   print "atom $onei[$k] $charge[$onei[$k]] $nsubst[$onei[$k]]\n" ;
            if ($nsubst[$onei[$k]] == 1) {
#             if ($bond[$onei[$k]][1] eq "s") {
#               $bond[$onei[$k]][1] = "d" ; #  force N=O
#               $charge[$onei[$k]] = 0 ;
#             }
            }
          }
          $val = $neb ;
        }
        if ($oneib == 1) {  # N+O-
          if ($neb == 3) {  
            $typ[$i] = "N2" ;
          }
          else {
            $typ[$i] = "N1" ;
          }
          $charge[$i] = 1 ;
          $typ[$onei[1]] = "O1" ;
          $charge[$onei[1]] = -1 ;
          $val = $neb ;
        }
      }
    }
  }

  if ($elem[$i] eq "S") {
    if ($nsubst[$i] < 2) {
      if ($neb == 1) {
        $val = 1 ;
      }
      else {
        $val = 0 ;
      }
    }
    elsif ($nsubst[$i] == 2) {
      if ($neb == 3) {
        $val = 1 ;
      }
      else {
        $val = 0 ;
      }
    }
    elsif ($nsubst[$i] == 3) {
      if ($neb == 5) {
        $val = 1 ;
      }
      else {
        $val = 0 ;
      }
    }
    elsif ($nsubst[$i] == 4) {
      $val = 0 ;
    }
    else {
      $val = 0 ;
    }
  }

  if ($elem[$i] eq "P") {
    if ($neb < 4) {   
      if ($nsubst[$i] < 3) {   #
        $val = 3 ;   # PX3
        $typ[$i] = "P" ;
      }
      else {
        $val = 5 ;   # PX5
        $typ[$i] = "P5" ;
      }
    }
    elsif ($neb == 4 && $nsubst[$i] == 4) {  # [P+]X4 
      $val = 4 ;
      $charge[$i] = 1 ; 
      $typ[$i] = "P5" ;
    }
    else {  # PX5
      $val = 5 ;
      $typ[$i] = "P5" ;
    }
  }


  if ($neb > 0) {
    $val = $val - $neb ;
  }
  else {
    $val = 0 ;
  }

#   $val is the number of hydrogens to be added
#  print "atom $i $elem[$i] nsub $nsubst[$i] bonds $neb  addhyd $val\n" ;

# add hydrogens
  if ($val > 0) {
    for ($j=1 ; $j <= $val ; $j++) {
      $newato++ ;               # number of new hydrogen atom
      $elem[$newato] = "H" ;
      $typ[$newato] = "**" ;
      $charge[$newato] = 0 ;
      $nsubst[$newato] = 1 ;
      $neib[$newato][1] = $i ;
      $bond[$newato][1] = "s" ;
      $nsubst[$i]++ ;           # update heavy atom
      $neib[$i][$nsubst[$i]] = $newato ;
      $bond[$i][$nsubst[$i]] = "s" ;
    }
  }

}
$natoms = $newato ; # update number of atoms

# generate topological bond matrix
for ($i=1 ; $i <= $natoms ; $i++) {
  for ($j=1 ; $j <= $natoms ; $j++) {
    $nbo[$i][$j] = 0 ;
  }
}
for ($i=1 ; $i <= $natoms ; $i++) {
  $nn = $nsubst[$i] ;
  for ($j=1 ; $j <= $nn ; $j++) {
    $kk = $neib[$i][$j] ;
    $nbo[$i][$kk] = $bond[$i][$j] ; # entries: s,d,t,a
  }
}

# ring perception is done according to Figuera's approach
# see J.Chem.Inf.Comput.Sci. 36 (1996) pp.986-991
# init arrays
  $nr3 = 0 ;  # number of 3-membered rings
  $nr4 = 0 ;  # number of 4-membered rings
  $nr5 = 0 ;  # number of 5-membered rings
  $nr6 = 0 ;  # number of 6-membered rings
  $nra = 0 ;  # number of aromatic rings
for ($i=1 ; $i <= $natoms ; $i++) {
  $ring3[$i] = 0 ;  # atom is part of a 3-membered ring
  $ring4[$i] = 0 ;  # atom is part of a 4-membered ring
  $ring5[$i] = 0 ;  # atom is part of a 5-membered ring
  $ring6[$i] = 0 ;  # atom is part of a 6-membered ring
  # Hier noch moeglich groesse Ringe (7,8) fuer 3D-Generierung
  $aring[$i] = 0 ;  # atom is part of a Hueckel aromatic ring
  $nring3[$i][0] = 0 ; # number and atoms of 3-membered rings
  $nring4[$i][0] = 0 ; # number and atoms 
  $nring5[$i][0] = 0 ; # number and atoms 
  $nring6[$i][0] = 0 ; # number and atoms 
  $atl[$i]   = 0 ;  # atom is not terminating
  $atpath[$i][0]  = 0 ;  # path of each atom 
  $atplen[$i] = 0     ;  # length of each atom path
}

# 1st generate list of non-terminating atoms (speed up ring search)
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($elem[$i] eq "H") {next}  
  if ($elem[$i] eq "F") {next}  
  if ($elem[$i] eq "Cl") {next}  
  if ($elem[$i] eq "Br") {next}  
  if ($elem[$i] eq "I") {next}  
  if ($nsubst[$i] <= 1) {next}
  $hneib = 0 ;
  for ($j=1 ; $j <=$nsubst[$i] ; $j++) {
    $kk = $neib[$i][$j] ;
    if ($elem[$kk] eq "H") {$hneib++} 
  }
  # -CH3, -NH2, -OH
  if ($elem[$i] eq "C" && $hneib >= 3) {next}
  if ($elem[$i] eq "N" && $hneib >= 2) {next}
  if ($elem[$i] eq "O" && $hneib >= 1) {next}
  $atl[$i] = 1 ; # not terminating atom or group
} 

# put all non-terminating atoms into the queue
#
##print "\nnon-terminating atoms are " ;
$atqul = 0 ; # length of atom queue array
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($atl[$i] == 1) {
    $atqul++ ;
    push(@ntqu,$i) ;
##    print "$i " ;
  }
}
##print "\n" ;

for ($rn=1 ; $rn <= $atqul ; $rn++) { # loop over all non-terminating atoms
  @atqu = () ;   # non-terminating atoms (array of type queue)
  @source = () ; # preceeding atom (array of type queue)
  $t1 = shift(@ntqu) ;
  push(@atqu,$t1) ;
  for ($t2=1 ; $t2 <= $natoms ; $t2++) {
    $atpath[$t2][0]  = 0 ;  # path of each atom 
    $atplen[$t2] = 0     ;  # length of each atom path
  }
  $frsrc = 0 ;
  $nxtr  = 0 ;
  for ($i=1 ; $i <= $atqul ; $i++) {  # breadth-first-search loop
    $front = shift(@atqu) ;   # step 1; get front node and its source
##   print "current front node is atom $front\n" ;
    if ($i == 1) { # init path of first atom
      $atpath[$front][0] = $front ;
      $atplen[$front] = 1 ;
    }
    else {$frsrc = shift(@source)} 
##   print "current source atom is $frsrc\n" ;
    for ($j=1 ; $j <= $nsubst[$front] ; $j++) { # step 2
      $m = $neib[$front][$j] ;
      if ($atl[$m] == 0) {next} # skip terminating atoms
      if ($m == $frsrc) {next}  # atom connected to front node and =! source
      if ($atpath[$m][0] == 0) {     
##       print "compute path for atom $m is " ; # path empty, compute path
        for ($k=0 ; $k < $atplen[$front] ; $k++) { 
          $atpath[$m][$k] = $atpath[$front][$k] ;
        } 
        $atpath[$m][$atplen[$front]] = $m ; # add atom m to path
        $atplen[$m] = $atplen[$front] + 1;
        for ($k=0 ; $k < $atplen[$m] ; $k++) {
##         print " $atpath[$m][$k]" ;
        }
##       print "\n" ;
        push(@atqu,$m) ;       # add atom $m to atom queue
        push(@source,$front) ; # add front node of atom $m to queue
##       print "adding atom $m to queue\n" ;
##       print "adding atom $front to source queue\n" ;
      }
      else { # compute intersection of path[front node] and path[atom m]
##       print "path to atom $m exits\n" ;
        $isectl = 0 ; # count of intersecting atoms
        for ($k1=0 ; $k1 < $atplen[$m] ; $k1++) {
          for ($k2=0 ; $k2 < $atplen[$front] ; $k2++) {
            if ($atpath[$m][$k1] == $atpath[$front][$k2]) { 
              $isectl++;
            }
          }
        }
        if ($isectl == 1) { # one atom in common => ring closure
#          print "ring closure between atoms $front and $m\n" ;
          # compute new ring set containing atoms of path[m]+path[front node]
          @newlt = () ; # init with empty array
          for ($k1=0 ; $k1 < $atplen[$m] ; $k1++) {
            push(@newlt, $atpath[$m][$k1]) ;
          }
          for ($k2=0 ; $k2 < $atplen[$front] ; $k2++) {
            push(@newlt, $atpath[$front][$k2]) ;
          }
          @newl1 = sort {$a <=> $b} @newlt ; # sort common atom path
##         print "common path is @newl1\n" ;
#         remove double atoms
          $rsize = 0 ;
          $atold = 0 ;
##          print "ring contains " ;
          for ($k1=0 ; $k1 < $atplen[$m]+$atplen[$front] ; $k1++) {
            $tt = shift(@newl1) ;
            if ($tt != $atold) {
              $newring[$rsize] = $tt ;
              $rsize++ ;
            }
            $atold = $tt ;
          }
          for ($k=0 ; $k < $rsize ; $k++) {
##            print " $newring[$k]"
          }
#         print "\n" ;       
          # assign atoms to corresponding ring systems
          # and check if this ring was already found 
          if ($rsize < 7) {
            if ($rsize == 3) {
              for ($k1=0 ; $k1 < 3 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring3[$tt] = 1 ;
              }
              # loop over all existing 3-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr3 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 3 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 3 ; $l2++) {
                    if ($tt == $nring3[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 3) {last} # ring already present
              }
              if ($atold < 3) { # ring not present yet
                $nr3++ ;
##                print "3-membered ring no $nr3 contains" ;
                for ($l1=0 ; $l1 < 3 ; $l1++) {
                  $nring3[$nr3][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
            if ($rsize == 4) {
              for ($k1=0 ; $k1 < 4 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring4[$tt] = 1 ;
              }
              # loop over all existing 4-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr4 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 4 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 4 ; $l2++) {
                    if ($tt == $nring4[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 4) {last} # ring already present
              }
              if ($atold < 4) { # ring not present yet
                $nr4++ ;
##                print "4-membered ring no $nr4 contains" ;
                for ($l1=0 ; $l1 < 4 ; $l1++) {
                  $nring4[$nr4][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
            if ($rsize == 5) {
              for ($k1=0 ; $k1 < 5 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring5[$tt] = 1 ;
              }
              # loop over all existing 5-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr5 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 5 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 5 ; $l2++) {
                    if ($tt == $nring5[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 5) {last} # ring already present
              }
              if ($atold < 5) { # ring not present yet
                $nr5++ ;
##              print "5-membered ring no $nr5 contains" ;    #remove pr
                for ($l1=0 ; $l1 < 5 ; $l1++) {
                  $nring5[$nr5][$l1] = $newring[$l1] ;
##                print " $newring[$l1]" ;                    #remove pr
                }
##              print "\n" ;                                  #remove pr
              }
            }
            if ($rsize == 6) {
              for ($k1=0 ; $k1 < 6 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring6[$tt] = 1 ;
              }
              # loop over all existing 6-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr6 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 6 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 6 ; $l2++) {
                    if ($tt == $nring6[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 6) {last} # ring already present
              }
              if ($atold < 6) { # ring not present yet
                $nr6++ ;
##                print "6-membered ring no $nr6 contains" ;
                for ($l1=0 ; $l1 < 6 ; $l1++) {
                  $nring6[$nr6][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
          }
          $nxtr = 1 ;
          last ; # start with next atom as root
        }
        if ($nxtr == 1) {last} 
      }
      if ($nxtr == 1) {last} 
    }
    if ($nxtr == 1) {last} # start from next non-terminating atom
  }
# print "###\n" ;
} # end of loop over all non-terminating atoms

# assign aromatic rings  
# $aring[$i] = 1 ;  # atom is part of a Hueckel aromatic ring
#
# Hueckel aromatic rings: only 5- and 6-membered rings 
# check for  
# 6-membered ring: 3 double bonds of C and N, but not exocyclic
# 5-membered ring: 2 double bonds of C and N and
#                    O with two neighbors or
#                    S with two neighbors
#
# This approach does not work for case such as
#         
#   //\ //\   Here, only the right ring is recognizes as
#   |  |  ||  being aromatic, since the left one possess
#   \\/ \\/   two exocyclic double bonds.
#     

for ($i=1 ; $i <= $nr6 ; $i++) { # check all 6-membered rings
  $dbc = 0 ; # double bonds of carbons
  $dbn = 0 ; # double bonds of nitrogens
  for ($j=0 ; $j < 6 ; $j++) { 
    $tt = $nring6[$i][$j] ;
    if ($elem[$tt] eq "O") {last} # these elements do not occur 
    if ($elem[$tt] eq "S") {last} # in 6-membered aromatic rings
    if ($elem[$tt] eq "P") {last}
    if ($elem[$tt] eq "Si") {last}
    if ($elem[$tt] eq "C") {
      if ($nsubst[$tt] != 3) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 3 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 6 ; $k2++) {
              if ($nn == $nring6[$i][$k2]) {$dbc++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "N") {
      if ($nsubst[$tt] != 2) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 2 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 6 ; $k2++) {
              if ($nn == $nring6[$i][$k2]) {$dbn++}
            }
          }
        }
      }
    }
  }
# print "6-membered ring no $i contains $dbc C= and $dbn N= bonds\n" ;
  if ($dbc+$dbn == 6) { # aromatic ring  (double bonds counted twice)
    $nra++ ;
    for ($l1=0 ; $l1 < 6 ; $l1++) {
      $tt = $nring6[$i][$l1] ;
      $aring[$tt] = 1 ;  
    }
    if ($aromat == 1) {
      for ($l1=0 ; $l1 < 6 ; $l1++) { # update bonding information
        $tt = $nring6[$i][$l1] ;
        for ($k1=1 ; $k1 <= $nsubst[$tt] ; $k1++) {
          $nn = $neib[$tt][$k1] ;
          for ($k2=0 ; $k2 < 6 ; $k2++) {
            if ($nn == $nring6[$i][$k2]) {
              $nbo[$tt][$nring6[$i][$k2]] = "a" ;
              $bond[$tt][$k1] = "a" ;
            }
          }
        }
      }
    }
  }
}

for ($i=1 ; $i <= $nr5 ; $i++) { # check all 5-membered rings
  $dbc = 0 ; # double bonds of carbons
  $dbn = 0 ; # double bonds of nitrogens
  $dbo = 0 ; # count of oxygens and sulfurs with two substituents
  for ($j=0 ; $j < 5 ; $j++) { 
    $tt = $nring5[$i][$j] ;
    if ($elem[$tt] eq "P")  {last} # these elements do not occur 
    if ($elem[$tt] eq "Si") {last} # in 6-membered aromatic rings
    if ($elem[$tt] eq "C") {
      if ($nsubst[$tt] != 3) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 3 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 5 ; $k2++) {
              if ($nn == $nring5[$i][$k2]) {$dbc++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "N") {
      if ($nsubst[$tt] != 2) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 2 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 5 ; $k2++) {
              if ($nn == $nring5[$i][$k2]) {$dbn++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "O" || $elem[$tt] eq "S") {
      if ($nsubst[$tt] == 2) {$dbo++}  
    }
  }
# print "5-membered ring no $i contains $dbc C= and $dbn N= bonds\n" ;
# print "                     and $dbo oxygens and sulfurs\n" ;
  if ($dbc+$dbn == 4 && $dbo == 1) { # aromatic ring  
    $nra++ ;
    for ($l1=0 ; $l1 < 5 ; $l1++) {
      $tt = $nring5[$i][$l1] ;
      $aring[$tt] = 1 ;  
    }
    if ($aromat == 1) {
      for ($l1=0 ; $l1 < 5 ; $l1++) { # update bonding information
        $tt = $nring5[$i][$l1] ;
        for ($k1=1 ; $k1 <= $nsubst[$tt] ; $k1++) {
          $nn = $neib[$tt][$k1] ;
          for ($k2=0 ; $k2 < 5 ; $k2++) {
            if ($nn == $nring5[$i][$k2]) {
              $nbo[$tt][$nring5[$i][$k2]] = "a" ;
              $bond[$tt][$k1] = "a" ;
            }
          }
        }
      }
    }
  }
}

# check ring assignment   
#print "\natom  3ring  4ring  5ring  6ring  aromat\n" ;
#for ($k=1 ; $k <= $natoms ; $k++) {
#  print "$k\t$ring3[$k]\t$ring4[$k]\t$ring5[$k]\t$ring6[$k]\t$aring[$k]\n" ;
#}
#print "\n" ;


for ($i=1 ; $i <= $natoms ; $i++) {
  if ($typ[$i] ne "**") {next}      # type already assigned
  $nn = $nsubst[$i] ;               # number of neighbors
  # The sequence is ordered in decreasing frequency of elements 
  # found in molecules to speed up things somewhat

  if ($elem[$i] eq "C" && $nn > 0) {                        # carbon
    if ($nn == 4) {                                         # sp3
      $typ[$i] = "C4" ;
      if ($ring4[$i] == 1) {$typ[$i] = "CB"} # cyclobutane
      if ($ring3[$i] == 1) {$typ[$i] = "CP"} # cyclopropane
      # atom type for 3-membered ring superseeds larger rings
      next ;
    }
    if ($nn == 3) {                                         # sp2
      $typ[$i] = "C3" ;
      if ($ring4[$i] == 1) {$typ[$i] = "CC"} # cyclobutene
      if ($ring3[$i] == 1) {$typ[$i] = "CZ"} # cyclopropene
      $oneib = 0 ;         # >C=O ?
      $nneib = 0 ;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "O" && $bond[$i][$j] eq "d") {
          $oneib++ ;
        }
        if ($elem[$kk] eq "N") {$nneib++}
      }
      if ($oneib > 0) {
        $typ[$i] = "CO" ;
        if ($ring4[$i] == 1) {               # cyclobutanone
          if ($nneib < 1) {$typ[$i] = "CD"}  # except lactim
        }
        if ($ring3[$i] == 1) {               # cyclopropanone
          if ($nneib < 1) {$typ[$i] = "CE"}  # except lactim
        }
        next ;
      }
      if ($aring[$i] == 1) {  # aromatic 
        $typ[$i] = "CA" ;
        next ;
      }
    }
    if ($nn == 2) {                                         # sp
      $typ[$i] = "C2" ;
      next ;
    }
    if ($nn == 1) {   # sp isonitrile [N+]#[C-]
      $typ[$i] = "C2" ;
      next ;
    }
    if ($charge[$i] == 1) {  # carbonium
      $typ[$i] = "C+" ;
      next ;
    }
    if ($typ[$i] eq "**") {  # default
      $typ[$i] = "C?" ;
      next ;
    } 
  }

  if ($elem[$i] eq "O" && $nn > 0) {                        # oxygen
    if ($nn == 2) { 
      $typ[$i] = "O2" ;
      if ($ring3[$i] == 1) {
        $typ[$i] = "OE" ;      # expoxide
        next ;
      }
      if ($ring5[$i] == 1) {
        $oneib = 0;
        for ($j=1 ; $j <= $nn ; $j++) {
          $kk = $neib[$i][$j] ;
#         if ($elem[$kk] eq "C" && $nsubst[$kk] == 2) {  corrected 10.01.14
          if ($elem[$kk] eq "C" && $nsubst[$kk] == 3) {
            $oneib++ ;
          }
        }
#       check if aromatic ring is furane (OF) or contains N or S (O2)
        if ($oneib == 2 && $aring[$i] == 1) {  # Car-O-Car
          for ($j=1 ; $j <= $nr5 ; $j++) {
            $m5 = 0 ;
            for ($k=0 ; $k < 5 ; $k++) {
              $kk = $nring5[$j][$k] ;
              if ($kk == $i) {$m5++}
              if ($kk == $neib[$i][1]) {$m5++}
              if ($kk == $neib[$i][2]) {$m5++}
            }
            if ($m5 == 3) {  # O and its neighbors are in ring $j
              $m5 = $j ;
              last ;
            }
          }
          $nns = 0 ;
          for ($j=0 ; $j < 5 ; $j++) {
            $kk = $nring5[$m5][$j] ;
            if ($elem[$kk] eq "N" || $elem[$kk] eq "S") {$nns++}
          }
          if ($nns == 0) {
            $typ[$i] = "OF" ;  # furane
          }
          else {
            $typ[$i] = "O2" ;  # other heteroaromatic 5-ring
          }
          next ;
        }
      }
    }
    if ($nn == 1) { 
      $kk = $neib[$i][1] ;
      if ($elem[$kk] eq "C" || $elem[$kk] eq "S" || $elem[$kk] eq "P") {
        if ($bond[$i][1] eq "d") {              # acidic oxygens
          $typ[$i] = "O1" ; # carbonyl oxygen
        }
        else {
          $typ[$i] = "OC" ; # carboxylate oxygen
        }
        next ;
      }
      if ($elem[$kk] eq "N") {
        if ($bond[$i][1] eq "s") {
          $typ[$i] = "ON" ;   # amine oxide oxygen
        }
        else {
          $typ[$i] = "O1" ;   # N=O, N(=O)(=O) 
        }
        next ;
      }
    }
    if ($typ[$i] eq "**") {  # default
      $typ[$i] = "O?" ;
      next ;
    } 
  }

  if ($elem[$i] eq "S" && $nn > 0) {                        # sulfur
    if ($nn == 1) {
      $typ[$i] = "S2" ;
      next ;
    }
    if ($nn == 2) {
      $typ[$i] = "S2" ;
      $oneib = 0;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "C" && $nsubst[$kk] == 3) {
          $oneib++ ;
        }
      }
      if ($oneib == 2 && $aring[$i] == 1) {  # Car-S-Car
        $typ[$i] = "SA" ;  # thiophen
        next ;
      }
    }
    if ($nn == 3) {
      $oneib = 0;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "O" && $nsubst[$kk] == 1) {
          $oneib++ ;
        }
      }
      if ($oneib > 0) {  
        $typ[$i] = "SO" ;  # >S=O
        next ;
      }
    }
    if ($nn == 4) {
      $typ[$i] = "S4" ;  # >S(=O)(=O)
      next ;
    }
    if ($typ[$i] eq "**") {  # default
      $typ[$i] = "S?" ;
      next ;
    } 
  }

  if ($elem[$i] eq "N" && $nn > 0) {                        # nitrogen
    if ($nn == 4) {
      $typ[$i] = "N4" ;  # ammonium
      next ;
    }
    if ($nn == 1) {
      $typ[$i] = "N1" ;  # sp1
#     if ($hyperc == 1) {
#       if ($bond[$i][1] eq "d" && $elem[$neib[$i][1]] eq "N") {
#         $typ[$i] = "NA" ; # -N=N=NA assignment of Hyperchem
#       }
#     }
      next ;
    }
    if ($nn == 2 || $nn ==3) {
      if ($nn == 3) {$typ[$i] = "N3"}  # default 
      if ($nn == 2) {$typ[$i] = "N2"}  # default
      $hneib = 0 ;
      $oneib = 0 ;
      $nneib = 0 ;
      $cneib = 0 ;
      $sneib = 0 ;
      for ($j=1 ; $j <= 3 ; $j++) { 
        $hnei[$j] = 0 ; 
        $onei[$j] = 0 ; 
        $nnei[$j] = 0 ; 
        $cnei[$j] = 0 ; 
        $snei[$j] = 0 ;
      }
      for ($j=1 ; $j <= $nn ; $j++) { # determine kind of substituents
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "H") {
          $hneib++ ;
          $hnei[$hneib] = $kk ;
        }
        if ($elem[$kk] eq "O") {
          $oneib++ ;
          $onei[$oneib] = $kk ;
        }
        if ($elem[$kk] eq "N") {
          $nneib++ ;
          $nnei[$nneib] = $kk ;
        }
        if ($elem[$kk] eq "C") {
          $cneib++ ;
          $cnei[$cneib] = $kk ;
        }
        if ($elem[$kk] eq "S") {
          $sneib++ ;
          $snei[$sneib] = $kk ;
        }
      }
      if ($nn == 3) { 
        if ($hneib > 0) {
          if ($cneib == 1) {
            if ($nsubst[$cnei[1]] > 3) {
              $typ[$i] = "N3" ; # nonplanar amine or amide
              next ;
            }
            else {
              $typ[$i] = "N2" ; # planar amine or amide
              next ;
            }
          }
          if ($cneib > 1) {
            if ($cneib == 2 && $hneib == 1) {
              if ($nsubst[$cnei[1]] == 3 || $nsubst[$cnei[2]] == 3) {
                $typ[$i] = "N2" ; # planar amine 
#               check for pyrrole ring
                if ($ring5[$i] == 1) {
                  for ($j=1 ; $j <= $nr5 ; $j++) {
                    $m5 = 0 ;
                    for ($k=0 ; $k < 5 ; $k++) {
                      $kk = $nring5[$j][$k] ;
                      if ($kk == $i) {$m5++}
                      if ($kk == $cnei[1]) {$m5++}
                      if ($kk == $cnei[2]) {$m5++}
                    }
                    if ($m5 == 3) {  # N and its C neighbors are in ring $j
                      $m5 = $j ;
                    last ;
                    }
                  }
                  $nos = 0 ;
                  for ($j=0 ; $j < 5 ; $j++) {
                    $kk = $nring5[$m5][$j] ;
                    if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
                  }
                  if ($nos == 0) {
                    $typ[$i] = "NP" ;   # pyrrole
                    next ;
                  }
                }
              }
              next ;
            }
            for ($k=1 ; $k <= $cneib ; $k++) { # determine kind of substituents
              $tt = $cnei[$k] ;
              if ($nsubst[$tt] == 3) {
                for ($k1=1 ; $k1 <= 3 ; $k1++) {
                  $oo = $neib[$tt][$k1] ;
                  if ($elem[$oo] eq "O" && $bond[$tt][$k1] eq "d") {
                    $typ[$i] = "N2" ; # peptide
                    next ;
                  }
                }
              }
            }
          }
          if ($sneib > 0) {
            if ($nsubst[$snei[1]] > 2) {
              $typ[$i] = "N2" ;
              next ;
            }
            if ($sneib == 2) {
              if ($nsubst[$snei[1]] == 2 && $nsubst[$snei[2]] == 2) {
                $typ[$i] = "N3" ;
                next ;
              }
              else {
                $typ[$i] = "N2" ;
                next ;
              }
            }
          }
          if ($nneib > 0) {
            $typ[$i] = "N2" ; 
            next ;
          }

        }
        else {
          if ($cneib > 0) {
            $typ[$i] = "N2" ; 
            next ;
          }
        }
      }
      if ($nn == 2) { 
        if ($oneib == 2) {     # NO2
          $typ[$i] = "NO" ;    # nitro 
          next ;
        }   
        if ($bond[$i][1] eq "t" && $bond[$i][2] eq "s") {
          $typ[$i] = "N1" ;    #   #[N-]  e.g. isonitrile
          next ;
        }
        if ($bond[$i][1] eq "s" && $bond[$i][2] eq "t") {
          $typ[$i] = "N1" ;    #   [N-]#[C+]  isonitrile
          next ;
        }
        if ($nneib == 2) {     
          if ($bond[$i][1] eq "d" && $bond[$i][2] eq "d") {
            # middle atom in -N=N=N 
#           if ($hyperc == 1) {
              $typ[$i] = "N1" ;  # Hyperchem assigns N1
#           }
#           else {
#             $typ[$i] = "NZ" ;  # MM3 assigns NZ
#           }
            next ; 
          }
          else {
            $typ[$i] = "NA" ;    # middle atom in -N=N-N= 
            next ;
          }
        }
        if ($cneib == 2) {     
          if ($bond[$i][1] eq "d" && $bond[$i][2] eq "s") {
            if ($nsubst[$cnei[2]] == 4) {
              $typ[$i] = "NC" ;    # imine C=N-C
              next ;
            }
            elsif ($nsubst[$cnei[2]] == 3) {
              $typ[$i] = "N2" ;    # C=N-Car 
              next ;
            }
            else {
              $typ[$i] = "NA" ;    
              next ;
            }
          }
          if ($bond[$i][1] eq "s" && $bond[$i][2] eq "d") {
            if ($nsubst[$cnei[1]] == 4) {
              $typ[$i] = "NC" ;    # imine C-N=C
              next ;
            }
            elsif ($nsubst[$cnei[1]] == 3) {
              $typ[$i] = "N2" ;    # Car-N=C    
              next ;
            }
            else {
              $typ[$i] = "NA" ;    
              next ;
            }
          }

#         check for pyrrole if explicit double bonds are present
          if ($ring5[$i] == 1) {
#   print "nsub N $nsubst[$i] C1 $nsubst[$cnei[1]] C2 $nsubst[$cnei[2]] \n" ;
            if ($nsubst[$cnei[1]] == 3 && $nsubst[$cnei[2]] == 3 &&                             $nsubst[$i] == 3) {
              for ($j=1 ; $j <= $nr5 ; $j++) {
                $m5 = 0 ;
                for ($k=0 ; $k < 5 ; $k++) {
                  $kk = $nring5[$j][$k] ;
                  if ($kk == $i) {$m5++}
                  if ($kk == $cnei[1]) {$m5++}
                  if ($kk == $cnei[2]) {$m5++}
                }
                if ($m5 == 3) {  # N and its C neighbors are in ring $j
                  $m5 = $j ;
                last ;
                }
              }
              $nos = 0 ;
              for ($j=0 ; $j < 5 ; $j++) {
                $kk = $nring5[$m5][$j] ;
                if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
              }
              if ($nos == 0) {
                $typ[$i] = "NP" ;   # pyrrole
                next ;
              }
            }
          }



          if ($aring[$cnei[1]] == 1 && $aring[$cnei[2]] == 1) {
            if ($ring5[$i] == 1) {
#             check if pyrolle (NP) or other heteroaromatic ring (NA)
              for ($j=1 ; $j <= $nr5 ; $j++) {
                $m5 = 0 ;
                for ($k=0 ; $k < 5 ; $k++) {
                  $kk = $nring5[$j][$k] ;
                  if ($kk == $i) {$m5++}
                  if ($kk == $cnei[1]) {$m5++}
                  if ($kk == $cnei[2]) {$m5++}
                }
                if ($m5 == 3) {  # N and its C neighbors are in ring $j
                  $m5 = $j ;
                  last ;
                }
              }
              $nos = 0 ;
              for ($j=0 ; $j < 5 ; $j++) {
                $kk = $nring5[$m5][$j] ;
                if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
              }
              if ($nos == 0) {
                $typ[$i] = "NP" ;   # pyrrole
              }
              else {
                $typ[$i] = "NA" ;   # other aromatic nitrogen
              }
              next ;
            }
          }
        }
        if ($nneib == 1 && $oneib == 1) {
          if ($nbo[$i][$nnei[1]] eq "d" && $nbo[$i][$onei[1]] eq "s") {
            $typ[$i] = "NB" ;  # middle atom in azoxy -N=N-O-
            next ;
          }
        }
        if ($cneib == 1 && $oneib == 1) {
          if ($nbo[$i][$cnei[1]] eq "d" && $nbo[$i][$onei[1]] eq "s") {
            $typ[$i] = "NC" ;  # oxime -C=N-OH 
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "a" && $nbo[$i][$onei[1]] eq "a") {
            $typ[$i] = "NC" ;  # aromatic oxazole -C=N-O-
            next ;
          }
        }
        if ($cneib == 1 && $nneib == 1) {
          if ($nbo[$i][$cnei[1]] eq "d" && $nbo[$i][$nnei[1]] eq "s") {
            $typ[$i] = "NA" ;
#           $typ[$i] = "N2" ;
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "s" && $nbo[$i][$nnei[1]] eq "d") {
            $typ[$i] = "NA" ;
#           $typ[$i] = "N2" ;
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "a" && $nbo[$i][$nnei[1]] eq "a") {
            $typ[$i] = "NA" ;
            next ;
          }
        }
      }
    }
    if ($typ[$i] eq "**") {  # default
      $typ[$i] = "N?" ;
      next ;
    } 
  }

}

# Now assign hydrogen atoms types
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($typ[$i] ne "**") {next}      # type already assigned
  $nn = $nsubst[$i] ;               # number of neighbors
  if ($elem[$i] eq "H" && $nn > 0) {                        # hydrogen
    $kk = $neib[$i][1] ;    # neighbor        
    $neibel = $elem[$kk] ;  # element of neighbor
    if ($neibel ne "N" && $neibel ne "O") {
      $typ[$i] = "H" ;    # H attached to C,Si,B,P,S,halogens
      next ;
    }
    if ($neibel eq "N") {           # connected to nitrogen     
      $typ[$i] = "HN" ;             # default
      if ($nsubst[$kk] > 3) {
        $typ[$i] = "HB" ;           # ammonium
        next ;
      }
      if ($typ[$kk] eq "N3" || $typ[$kk] eq "NC") {
        $typ[$i] = "HN" ;           # nonplanar amine 
        next ;
      }
      if ($typ[$kk] eq "N2" || $typ[$kk] eq "NP") {
        $typ[$i] = "HV" ;           # planar amine 
        next ;
      }
      next ;
    }
    if ($neibel eq "O") {           # connected to oxygen
                            # atom types HO,HX,HV
      $typ[$i] = "HO" ;             # default
      $l1 = $neib[$kk][1] ; # next neighbor to O
      $l2 = $neib[$kk][2] ;
      if ($l1 == $i) {
        $ll = $l2 ;
      }
      else {
        $ll = $l1 ;
      }
      if ($elem[$ll] ne "C") {next}
      if ($nsubst[$ll] == 4) {
        $typ[$i] = "HO" ;  # alcohol oxygen C-C-O-H
        next ;
      }
      if ($nsubst[$ll] != 3) {next}
      for ($j=1 ; $j <= $nsubst[$ll] ; $j++) {
        $mm = $neib[$ll][$j] ;
        if ($elem[$mm] eq "C" && $nbo[$ll][$mm] eq "d") {
          $typ[$i] = "HV" ; # enol oxygen  C=C-O-H 
          next ;
        }
        if ($elem[$mm] eq "O" && $nbo[$ll][$mm] eq "d") {
          $typ[$i] = "HX" ; # carboxylic oxygen  C(=O)OH
          next ;
        }
      }
      next ;
    }
  }
}

### testprint for hinfile
#  print "\n" ;
#  for ($j=1 ; $j <= $natoms ; $j++) {
#    print "atom $j - $elem[$j] $typ[$j] - $charge[$j] 0 0 0 $nsubst[$j] " ;
#    if ($nsubst[$j] > 0) {
#      for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
#        print "$neib[$j][$k] $bond[$j][$k] ";
#      }
#    }
#    print "\n" ;
#  }   
#  print "\n" ;

}

}

close (INPU) or die "error closing $INFILE : $!\n";

# $natoms     # number of atoms
# $elem[i]    # element
# $typ[i]     # atom type
# $nsubst[i]  # numer of neighbors including hydrogens
# $neib[i][j] # list of neighbors


# List of available atom types
# this block incl. $maxt must be identical to that in genmat.pl !
# here: those from the MM+ force field plus defaults for each element

$t[0]   = "H?" ; # H any hydrogen atom
$t[1]   = "C4" ; # Csp3
$t[2]   = "C3" ; # Csp2 alkene
$t[3]   = "CO" ; # Csp2 carbonyl
$t[4]   = "C2" ; # Csp  alkyne and C=C=O
$t[5]   = "CP" ; # Csp3 cyclopropane
$t[6]   = "C+" ; # Carboniumion
$t[7]   = "CZ" ; # Csp2 cyclopropene
$t[8]   = "CA" ; # Csp2 aromatic
$t[9]   = "CB" ; # Csp3 cyclobutane
$t[10]  = "CC" ; # Csp2 cyclobutene
$t[11]  = "CD" ; # Csp2 cyclobutanone C=O
$t[12]  = "CE" ; # Csp2 cyclopropanone C==
$t[13]  = "C?" ; # C any other carbon

$t[14]  = "N3" ; # Nsp3 
$t[15]  = "N2" ; # Nsp2 amide
$t[16]  = "N1" ; # Nssp
$t[17]  = "NA" ; # N -N= azo, pyridine
$t[18]  = "NH" ; # Nsp3 ammonium 
$t[19]  = "NP" ; # Nsp2 pyrrole
$t[20]  = "NB" ; # N -N=N-O azoxy
$t[21]  = "NZ" ; # N azide central N
$t[22]  = "NO" ; # N NO2 nitro
$t[23]  = "NC" ; # N =N- imine, oxime
$t[24]  = "N?" ; # N any other nitrogen

$t[25]  = "O2" ; # O C-O-H, C-O-C
$t[26]  = "O1" ; # O =O carbonyl
$t[27]  = "OF" ; # Osp2 furan 
$t[28]  = "OC" ; # O carboxylate 
$t[29]  = "OE" ; # O epoxy
$t[30]  = "ON" ; # O amine oxide 
$t[31]  = "O?" ; # O any other oxygen

$t[32]  = "B3" ; # B trigonal  
$t[33]  = "B4" ; # B tetragonal

$t[34]  = "SI" ; # Si silane  

$t[35]  = "P" ;  # P >P- phosphine
$t[36]  = "P5" ; # P (V)
$t[37]  = "P?" ; # P any other phosphorus

$t[38]  = "S2" ; # S -S- sulfide
$t[39]  = "S+" ; # S+ >S+ Sulfonium
$t[40]  = "SO" ; # S >S=O sulfoxide
$t[41]  = "S4" ; # S >SO2 sulfone
$t[42]  = "SA" ; # Ssp2 thiophene
$t[43]  = "S?" ; # S any other sulfur

$t[44]  = "F"  ; # F 
$t[45]  = "CL" ; # Cl
$t[46]  = "BR" ; # Br
$t[47]  = "I"  ; # I

$t[48]  = "NN" ; # no atom 

$maxt = 48 ;

# end of atom types block

#
# set up matrices
#

for ($k=1 ; $k <= 5 ; $k++) {
  for ($i=0 ; $i <= $maxt; $i++) {
    for ($nl=1 ; $nl <= $maxt ; $nl++) {
      if ($k == 1) {$s[$nl][$i] = 0.0}
      if ($k == 2) {$s2[$nl][$i] = 0.0}
      if ($k == 3) {$s3[$nl][$i] = 0.0}
      if ($k == 4) {$s4[$nl][$i] = 0.0}
      if ($k == 5) {$s5[$nl][$i] = 0.0}
    }
  }
}

$p0[1]  =  0.0699344 ; # C4
$p0[2]  = -0.0064087 ;
$p0[3]  =  0.0029980 ;
$p0[4]  = -0.0039266 ;
$p0[5]  =  0.0015467 ;
$p0[6]  =  0.0 ; 
$p0[7]  =  0.0000334 ;
$p0[8]  = -0.0648419 ;
$p0[9]  = -0.0005799 ;
$p0[10] = -0.0001926 ;
$p0[11] = -0.0000565 ;
$p0[12] =  0.0 ; 
$p0[13] =  0.0 ;
$p0[14] =  0.0097141 ;
$p0[15] =  0.0140481 ;
$p0[16] = -0.0028588 ;
$p0[17] = -0.0043664 ;
$p0[18] = -0.0001130 ;
$p0[19] = -0.0000904 ;
$p0[20] =  0.0000502 ;
$p0[21] =  0.0 ;
$p0[22] = -0.0000565 ;
$p0[23] =  0.0007383 ;
$p0[24] =  0.0 ;
$p0[25] = -0.0106598 ;
$p0[26] = -0.0006456 ;
$p0[27] =  0.0004365 ;
$p0[28] =  0.0001690 ;
$p0[29] =  0.0000750 ;
$p0[30] =  0.0 ;
$p0[31] = -0.0001591 ;
$p0[32] = -0.0006047 ;
$p0[33] =  0.0 ;
$p0[34] = -0.0025820 ;
$p0[35] = -0.0002260 ;
$p0[36] = -0.0005758 ;
$p0[37] =  0.0 ;
$p0[38] =  0.0018685 ;
$p0[39] =  0.0000167 ;
$p0[40] = -0.0000715 ;
$p0[41] = -0.0010295 ;
$p0[42] =  0.0003611 ;
$p0[43] =  0.0 ;
$p0[44] =  0.0034341 ;
$p0[45] = -0.0002826 ;
$p0[46] = -0.0039239 ;
$p0[47] = -0.0005672 ;
$p0[48] =  0.0 ;


$s[1][1] = -0.279 ;
$s[1][2] = -0.163 ;
$s[1][3] = -0.076 ;
$s[1][4] = -0.103 ;
$s[1][5] = -0.746 ;
$s[1][8] =  0.162 ;
$s[1][9] = -1.373 ;
$s[1][10] = -2.673 ;
$s[1][14] =  0.194 ;
$s[1][15] =  0.889 ;
$s[1][16] =  5.871 ;
$s[1][17] =  0.285 ;
$s[1][18] = -9.293 ;
$s[1][19] =  0.518 ;
$s[1][25] =  0.130 ;
$s[1][29] =  0.243 ;
$s[1][32] = -2.571 ;
$s[1][34] = -2.948 ;
$s[1][36] = -0.229 ;
$s[1][38] =  0.410 ;
$s[1][39] =  7.736 ;
$s[1][40] = -0.028 ;
$s[1][41] =  0.667 ;
$s[1][44] =  0.113 ;
$s[1][45] = -0.981 ;
$s[1][46] = -9.094 ;
$s[1][47] = -8.318 ;
$s[2][1] = -0.163 ;
$s[2][2] = -0.092 ;
$s[2][3] =  0.071 ;
$s[2][4] = -0.390 ;
$s[2][5] =  0.456 ;
$s[2][8] =  0.787 ;
$s[2][9] = -8.999 ;
$s[2][14] = -0.981 ;
$s[2][15] =  0.182 ;
$s[2][16] =  1.617 ;
$s[2][17] = -0.352 ;
$s[2][19] = -0.142 ;
$s[2][23] =  0.034 ;
$s[2][25] =  0.069 ;
$s[2][27] =  0.906 ;
$s[2][31] =  6.624 ;
$s[2][34] = -7.969 ;
$s[2][36] =  6.624 ;
$s[2][38] = -0.237 ;
$s[2][40] =  8.729 ;
$s[2][41] =  0.812 ;
$s[2][42] =  0.498 ;
$s[2][44] =  6.493 ;
$s[2][45] = -0.548 ;
$s[2][46] =  1.046 ;
$s[2][47] = -1.962 ;
$s[3][1] = -0.076 ;
$s[3][2] =  0.071 ;
$s[3][3] = -0.449 ;
$s[3][4] = -8.963 ;
$s[3][5] = -0.681 ;
$s[3][7] =  7.326 ;
$s[3][8] = -0.625 ;
$s[3][9] = -7.879 ;
$s[3][15] =  0.465 ;
$s[3][17] = -0.828 ;
$s[3][25] = -0.188 ;
$s[3][26] = -0.038 ;
$s[3][28] = -1.784 ;
$s[3][36] =  6.624 ;
$s[3][38] =  1.435 ;
$s[3][45] = -9.157 ;
$s[3][46] = -7.109 ;
$s[4][1] = -0.103 ;
$s[4][2] = -0.390 ;
$s[4][3] = -8.963 ;
$s[4][4] = -1.228 ;
$s[4][5] =  6.624 ;
$s[4][8] =  0.720 ;
$s[4][15] = -1.616 ;
$s[4][16] = -1.646 ;
$s[4][26] = -8.097 ;
$s[4][34] = -8.999 ;
$s[4][38] = -2.537 ;
$s[4][47] =  6.624 ;
$s[5][1] = -0.746 ;
$s[5][2] =  0.456 ;
$s[5][3] = -0.681 ;
$s[5][4] =  6.624 ;
$s[5][5] =  0.113 ;
$s[5][8] =  7.043 ;
$s[5][14] =  8.253 ;
$s[5][15] =  8.660 ;
$s[5][25] =  5.754 ;
$s[5][44] =  6.624 ;
$s[5][45] =  7.326 ;
$s[7][3] =  7.326 ;
$s[7][7] =  7.326 ;
$s[7][8] =  7.326 ;
$s[8][1] =  0.162 ;
$s[8][2] =  0.787 ;
$s[8][3] = -0.625 ;
$s[8][4] =  0.720 ;
$s[8][5] =  7.043 ;
$s[8][7] =  7.326 ;
$s[8][8] =  0.138 ;
$s[8][9] =  6.624 ;
$s[8][15] = -0.237 ;
$s[8][16] = -0.101 ;
$s[8][17] =  0.173 ;
$s[8][19] = -0.152 ;
$s[8][23] =  5.444 ;
$s[8][25] = -0.141 ;
$s[8][28] =  9.566 ;
$s[8][32] = -7.777 ;
$s[8][34] = -7.439 ;
$s[8][35] = -9.999 ;
$s[8][36] = -2.731 ;
$s[8][38] =  0.637 ;
$s[8][40] =  0.198 ;
$s[8][41] =  0.025 ;
$s[8][44] = -0.032 ;
$s[8][45] =  0.982 ;
$s[8][46] =  0.575 ;
$s[8][47] = -0.583 ;
$s[9][1] = -1.373 ;
$s[9][2] = -8.999 ;
$s[9][3] = -7.879 ;
$s[9][8] =  6.624 ;
$s[9][9] = -1.145 ;
$s[9][11] = -8.586 ;
$s[9][14] = -7.879 ;
$s[9][25] = -8.586 ;
$s[10][1] = -2.673 ;
$s[10][10] = -1.967 ;
$s[11][9] = -8.586 ;
$s[11][26] = -7.879 ;
$s[14][1] =  0.194 ;
$s[14][2] = -0.981 ;
$s[14][5] =  8.253 ;
$s[14][9] = -7.879 ;
$s[14][14] = -2.537 ;
$s[14][15] = -1.940 ;
$s[14][32] = -9.706 ;
$s[14][34] = -8.999 ;
$s[14][36] =  7.736 ;
$s[14][45] =  8.236 ;
$s[14][46] =  8.594 ;
$s[14][47] =  6.624 ;
$s[15][1] =  0.889 ;
$s[15][2] =  0.182 ;
$s[15][3] =  0.465 ;
$s[15][4] = -1.616 ;
$s[15][5] =  8.660 ;
$s[15][8] = -0.237 ;
$s[15][14] = -1.940 ;
$s[15][15] = -0.277 ;
$s[15][17] =  1.367 ;
$s[15][19] =  6.624 ;
$s[15][25] =  0.093 ;
$s[15][26] = -0.826 ;
$s[15][34] = -9.232 ;
$s[15][36] = -0.405 ;
$s[15][38] = -0.821 ;
$s[15][41] =  1.673 ;
$s[15][45] = -1.514 ;
$s[15][46] = -7.948 ;
$s[15][47] =  6.624 ;
$s[16][1] =  5.871 ;
$s[16][2] =  1.617 ;
$s[16][4] = -1.646 ;
$s[16][8] = -0.101 ;
$s[16][16] = -0.855 ;
$s[16][17] = -1.266 ;
$s[16][20] =  6.624 ;
$s[16][25] =  2.652 ;
$s[16][26] =  0.494 ;
$s[17][1] =  0.285 ;
$s[17][2] = -0.352 ;
$s[17][3] = -0.828 ;
$s[17][8] =  0.173 ;
$s[17][15] =  1.367 ;
$s[17][16] = -1.266 ;
$s[17][17] =  0.093 ;
$s[17][20] =  7.326 ;
$s[17][25] =  4.763 ;
$s[17][26] =  7.107 ;
$s[17][34] = -7.879 ;
$s[17][38] =  6.624 ;
$s[17][41] =  7.326 ;
$s[17][45] =  0.659 ;
$s[17][46] =  6.624 ;
$s[17][47] =  6.624 ;
$s[18][1] = -9.293 ;
$s[19][1] =  0.518 ;
$s[19][2] = -0.142 ;
$s[19][8] = -0.152 ;
$s[19][15] =  6.624 ;
$s[19][32] = -7.879 ;
$s[20][16] =  6.624 ;
$s[20][17] =  7.326 ;
$s[20][25] =  7.736 ;
$s[22][25] = -7.879 ;
$s[22][26] = -7.879 ;
$s[23][2] =  0.034 ;
$s[23][8] =  5.444 ;
$s[23][25] =  0.299 ;
$s[23][28] =  6.624 ;
$s[25][1] =  0.130 ;
$s[25][2] =  0.069 ;
$s[25][3] = -0.188 ;
$s[25][5] =  5.754 ;
$s[25][8] = -0.141 ;
$s[25][9] = -8.586 ;
$s[25][15] =  0.093 ;
$s[25][16] =  2.652 ;
$s[25][17] =  4.763 ;
$s[25][20] =  7.736 ;
$s[25][22] = -7.879 ;
$s[25][23] =  0.299 ;
$s[25][25] =  0.168 ;
$s[25][32] = -2.998 ;
$s[25][34] = -3.833 ;
$s[25][35] = -9.293 ;
$s[25][36] = -0.510 ;
$s[25][40] = -9.293 ;
$s[25][41] = -1.608 ;
$s[25][45] = -5.777 ;
$s[26][3] = -0.038 ;
$s[26][4] = -8.097 ;
$s[26][11] = -7.879 ;
$s[26][15] = -0.826 ;
$s[26][16] =  0.494 ;
$s[26][17] =  7.107 ;
$s[26][22] = -7.879 ;
$s[26][36] = -0.707 ;
$s[26][40] = -0.275 ;
$s[26][41] = -0.043 ;
$s[27][2] =  0.906 ;
$s[28][3] = -1.784 ;
$s[28][8] =  9.566 ;
$s[28][23] =  6.624 ;
$s[28][28] =  7.326 ;
$s[28][34] = -9.520 ;
$s[29][1] =  0.243 ;
$s[31][2] =  6.624 ;
$s[31][36] =  7.326 ;
$s[31][41] = -2.669 ;
$s[32][1] = -2.571 ;
$s[32][8] = -7.777 ;
$s[32][14] = -9.706 ;
$s[32][19] = -7.879 ;
$s[32][25] = -2.998 ;
$s[34][1] = -2.948 ;
$s[34][2] = -7.969 ;
$s[34][4] = -8.999 ;
$s[34][8] = -7.439 ;
$s[34][14] = -8.999 ;
$s[34][15] = -9.232 ;
$s[34][17] = -7.879 ;
$s[34][25] = -3.833 ;
$s[34][28] = -9.520 ;
$s[34][34] = -9.999 ;
$s[34][38] = -7.879 ;
$s[34][45] =-10.413 ;
$s[35][8] = -9.999 ;
$s[35][25] = -9.293 ;
$s[36][1] = -0.229 ;
$s[36][2] =  6.624 ;
$s[36][3] =  6.624 ;
$s[36][8] = -2.731 ;
$s[36][14] =  7.736 ;
$s[36][15] = -0.405 ;
$s[36][25] = -0.510 ;
$s[36][26] = -0.707 ;
$s[36][31] =  7.326 ;
$s[36][38] = -2.537 ;
$s[36][45] = -9.706 ;
$s[36][46] = -7.879 ;
$s[38][1] =  0.410 ;
$s[38][2] = -0.237 ;
$s[38][3] =  1.435 ;
$s[38][4] = -2.537 ;
$s[38][8] =  0.637 ;
$s[38][15] = -0.821 ;
$s[38][17] =  6.624 ;
$s[38][34] = -7.879 ;
$s[38][36] = -2.537 ;
$s[38][38] =  8.438 ;
$s[38][40] = -7.879 ;
$s[38][45] = -8.586 ;
$s[39][1] =  7.736 ;
$s[39][45] =  6.624 ;
$s[40][1] = -0.028 ;
$s[40][2] =  8.729 ;
$s[40][8] =  0.198 ;
$s[40][25] = -9.293 ;
$s[40][26] = -0.275 ;
$s[40][38] = -7.879 ;
$s[41][1] =  0.667 ;
$s[41][2] =  0.812 ;
$s[41][8] =  0.025 ;
$s[41][15] =  1.673 ;
$s[41][17] =  7.326 ;
$s[41][25] = -1.608 ;
$s[41][26] = -0.043 ;
$s[41][31] = -2.669 ;
$s[41][44] = -9.520 ;
$s[41][45] = -9.706 ;
$s[42][2] =  0.498 ;
$s[44][1] =  0.113 ;
$s[44][2] =  6.493 ;
$s[44][5] =  6.624 ;
$s[44][8] = -0.032 ;
$s[44][41] = -9.520 ;
$s[45][1] = -0.981 ;
$s[45][2] = -0.548 ;
$s[45][3] = -9.157 ;
$s[45][5] =  7.326 ;
$s[45][8] =  0.982 ;
$s[45][14] =  8.236 ;
$s[45][15] = -1.514 ;
$s[45][17] =  0.659 ;
$s[45][25] = -5.777 ;
$s[45][34] =-10.413 ;
$s[45][36] = -9.706 ;
$s[45][38] = -8.586 ;
$s[45][39] =  6.624 ;
$s[45][41] = -9.706 ;
$s[45][47] =  6.624 ;
$s[46][1] = -9.094 ;
$s[46][2] =  1.046 ;
$s[46][3] = -7.109 ;
$s[46][8] =  0.575 ;
$s[46][14] =  8.594 ;
$s[46][15] = -7.948 ;
$s[46][17] =  6.624 ;
$s[46][36] = -7.879 ;
$s[47][1] = -8.318 ;
$s[47][2] = -1.962 ;
$s[47][4] =  6.624 ;
$s[47][8] = -0.583 ;
$s[47][14] =  6.624 ;
$s[47][15] =  6.624 ;
$s[47][17] =  6.624 ;
$s[47][45] =  6.624 ;

$s2[1][1] = -0.314 ;
$s2[1][2] =  0.024 ;
$s2[1][3] =  0.303 ;
$s2[1][4] =  0.019 ;
$s2[1][5] = -0.701 ;
$s2[1][8] =  0.224 ;
$s2[1][9] = -1.632 ;
$s2[1][10] = -2.350 ;
$s2[1][11] = -7.508 ;
$s2[1][14] =  0.099 ;
$s2[1][15] =  0.867 ;
$s2[1][16] = -0.384 ;
$s2[1][17] =  0.061 ;
$s2[1][18] = -9.383 ;
$s2[1][19] =  0.253 ;
$s2[1][22] = -7.508 ;
$s2[1][23] =  1.049 ;
$s2[1][25] =  0.282 ;
$s2[1][26] = -0.002 ;
$s2[1][27] = -0.774 ;
$s2[1][28] = -9.684 ;
$s2[1][29] =  0.223 ;
$s2[1][31] = -2.032 ;
$s2[1][32] = -3.735 ;
$s2[1][34] = -3.023 ;
$s2[1][35] = -8.959 ;
$s2[1][36] =  0.239 ;
$s2[1][38] =  0.463 ;
$s2[1][39] =  6.213 ;
$s2[1][40] = -1.748 ;
$s2[1][41] =  0.831 ;
$s2[1][42] =  0.083 ;
$s2[1][44] =  1.223 ;
$s2[1][45] = -0.921 ;
$s2[1][46] =  0.215 ;
$s2[1][47] =  0.100 ;
$s2[2][1] =  0.024 ;
$s2[2][2] = -0.182 ;
$s2[2][3] =  0.128 ;
$s2[2][4] = -0.453 ;
$s2[2][5] =  0.399 ;
$s2[2][8] =  0.620 ;
$s2[2][9] = -9.807 ;
$s2[2][11] = -7.508 ;
$s2[2][14] = -0.356 ;
$s2[2][15] =  0.166 ;
$s2[2][16] =  0.783 ;
$s2[2][17] =  0.576 ;
$s2[2][19] = -0.360 ;
$s2[2][20] =  8.059 ;
$s2[2][23] =  0.051 ;
$s2[2][25] =  0.148 ;
$s2[2][26] =  0.187 ;
$s2[2][27] =  1.009 ;
$s2[2][28] =  9.528 ;
$s2[2][29] = -0.174 ;
$s2[2][32] = -8.233 ;
$s2[2][34] = -8.559 ;
$s2[2][36] = -1.295 ;
$s2[2][38] =  0.657 ;
$s2[2][40] = -0.174 ;
$s2[2][41] =  1.188 ;
$s2[2][42] =  0.577 ;
$s2[2][44] = -0.034 ;
$s2[2][45] = -0.744 ;
$s2[2][46] =  0.784 ;
$s2[2][47] = -1.731 ;
$s2[3][1] =  0.303 ;
$s2[3][2] =  0.128 ;
$s2[3][3] = -0.060 ;
$s2[3][4] = -1.628 ;
$s2[3][5] = -0.574 ;
$s2[3][7] =  6.927 ;
$s2[3][8] = -0.320 ;
$s2[3][9] = -9.383 ;
$s2[3][14] = -0.635 ;
$s2[3][15] =  1.232 ;
$s2[3][16] = -0.155 ;
$s2[3][17] = -0.203 ;
$s2[3][19] =  0.539 ;
$s2[3][23] =  1.635 ;
$s2[3][25] =  0.677 ;
$s2[3][26] = -0.357 ;
$s2[3][27] =  0.551 ;
$s2[3][28] =  6.927 ;
$s2[3][29] =  7.871 ;
$s2[3][34] = -9.609 ;
$s2[3][36] =  6.927 ;
$s2[3][38] =  1.034 ;
$s2[3][40] =  6.927 ;
$s2[3][41] =  0.979 ;
$s2[3][42] = -1.088 ;
$s2[3][44] = -1.901 ;
$s2[3][45] = -0.801 ;
$s2[3][46] = -9.742 ;
$s2[3][47] = -8.233 ;
$s2[4][1] =  0.019 ;
$s2[4][2] = -0.453 ;
$s2[4][3] = -1.628 ;
$s2[4][4] =-10.269 ;
$s2[4][5] =  7.345 ;
$s2[4][8] =  0.773 ;
$s2[4][14] = -0.975 ;
$s2[4][15] = -0.101 ;
$s2[4][17] = -3.170 ;
$s2[4][23] =  6.213 ;
$s2[4][25] = -0.173 ;
$s2[4][26] = -8.317 ;
$s2[4][34] = -8.233 ;
$s2[4][36] = -7.508 ;
$s2[4][38] =  6.213 ;
$s2[4][41] = -7.508 ;
$s2[4][42] = -7.508 ;
$s2[4][46] = -7.508 ;
$s2[4][47] =  6.213 ;
$s2[5][1] = -0.701 ;
$s2[5][2] =  0.399 ;
$s2[5][3] = -0.574 ;
$s2[5][4] =  7.345 ;
$s2[5][5] = -1.002 ;
$s2[5][8] =  2.042 ;
$s2[5][14] =  2.061 ;
$s2[5][15] =  9.464 ;
$s2[5][23] =  6.213 ;
$s2[5][25] =  0.153 ;
$s2[5][26] = -0.791 ;
$s2[5][38] =  6.213 ;
$s2[5][44] =  6.927 ;
$s2[5][45] =  7.641 ;
$s2[7][3] =  6.927 ;
$s2[7][7] =  6.927 ;
$s2[7][8] =  8.059 ;
$s2[7][26] =  6.927 ;
$s2[8][1] =  0.224 ;
$s2[8][2] =  0.620 ;
$s2[8][3] = -0.320 ;
$s2[8][4] =  0.773 ;
$s2[8][5] =  2.042 ;
$s2[8][7] =  8.059 ;
$s2[8][8] =  0.030 ;
$s2[8][9] =  7.641 ;
$s2[8][14] = -0.464 ;
$s2[8][15] = -0.059 ;
$s2[8][16] = -0.127 ;
$s2[8][17] = -0.049 ;
$s2[8][19] =  0.143 ;
$s2[8][23] =  0.032 ;
$s2[8][25] = -0.239 ;
$s2[8][26] = -0.497 ;
$s2[8][27] =  7.112 ;
$s2[8][28] =  6.648 ;
$s2[8][32] = -8.554 ;
$s2[8][34] = -7.866 ;
$s2[8][35] =-10.409 ;
$s2[8][36] = -1.872 ;
$s2[8][38] =  0.680 ;
$s2[8][40] =  0.482 ;
$s2[8][41] =  0.066 ;
$s2[8][44] =  0.088 ;
$s2[8][45] =  0.770 ;
$s2[8][46] =  0.503 ;
$s2[8][47] = -0.403 ;
$s2[9][1] = -1.632 ;
$s2[9][2] = -9.807 ;
$s2[9][3] = -9.383 ;
$s2[9][8] =  7.641 ;
$s2[9][9] = -1.212 ;
$s2[9][11] = -8.233 ;
$s2[9][14] = -0.592 ;
$s2[9][25] = -8.959 ;
$s2[9][26] = -8.658 ;
$s2[9][29] = -7.508 ;
$s2[10][1] = -2.350 ;
$s2[10][10] = -8.959 ;
$s2[11][1] = -7.508 ;
$s2[11][2] = -7.508 ;
$s2[11][9] = -8.233 ;
$s2[14][1] =  0.099 ;
$s2[14][2] = -0.356 ;
$s2[14][3] = -0.635 ;
$s2[14][4] = -0.975 ;
$s2[14][5] =  2.061 ;
$s2[14][8] = -0.464 ;
$s2[14][9] = -0.592 ;
$s2[14][14] = -0.906 ;
$s2[14][15] = -0.483 ;
$s2[14][17] =  5.436 ;
$s2[14][25] =  0.810 ;
$s2[14][38] =  7.467 ;
$s2[14][41] = -8.233 ;
$s2[15][1] =  0.867 ;
$s2[15][2] =  0.166 ;
$s2[15][3] =  1.232 ;
$s2[15][4] = -0.101 ;
$s2[15][5] =  9.464 ;
$s2[15][8] = -0.059 ;
$s2[15][14] = -0.483 ;
$s2[15][15] = -0.377 ;
$s2[15][16] =  0.387 ;
$s2[15][17] = -0.128 ;
$s2[15][20] =  6.927 ;
$s2[15][23] = -0.048 ;
$s2[15][25] =  1.810 ;
$s2[15][26] =  0.596 ;
$s2[15][28] =  6.213 ;
$s2[15][31] =  6.213 ;
$s2[15][32] =  6.213 ;
$s2[15][36] =  6.213 ;
$s2[15][38] = -0.543 ;
$s2[15][40] =  8.217 ;
$s2[15][41] =  1.513 ;
$s2[15][42] = -0.174 ;
$s2[15][45] = -0.415 ;
$s2[15][46] =  6.213 ;
$s2[16][1] = -0.384 ;
$s2[16][2] =  0.783 ;
$s2[16][3] = -0.155 ;
$s2[16][8] = -0.127 ;
$s2[16][15] =  0.387 ;
$s2[16][17] = -2.020 ;
$s2[16][23] = -7.508 ;
$s2[16][25] =  5.645 ;
$s2[16][27] =  9.607 ;
$s2[16][34] = -8.233 ;
$s2[16][38] = -1.306 ;
$s2[16][42] =  6.213 ;
$s2[17][1] =  0.061 ;
$s2[17][2] =  0.576 ;
$s2[17][3] = -0.203 ;
$s2[17][4] = -3.170 ;
$s2[17][8] = -0.049 ;
$s2[17][14] =  5.436 ;
$s2[17][15] = -0.128 ;
$s2[17][16] = -2.020 ;
$s2[17][17] = -0.262 ;
$s2[17][19] = -1.295 ;
$s2[17][20] =  6.213 ;
$s2[17][23] =  6.213 ;
$s2[17][25] = -1.020 ;
$s2[17][26] = -0.656 ;
$s2[17][38] = -1.726 ;
$s2[17][41] =  6.213 ;
$s2[17][42] =  6.213 ;
$s2[17][45] = -1.122 ;
$s2[18][1] = -9.383 ;
$s2[19][1] =  0.253 ;
$s2[19][2] = -0.360 ;
$s2[19][3] =  0.539 ;
$s2[19][8] =  0.143 ;
$s2[19][17] = -1.295 ;
$s2[19][25] =  6.213 ;
$s2[20][2] =  8.059 ;
$s2[20][15] =  6.927 ;
$s2[20][17] =  6.213 ;
$s2[22][1] = -7.508 ;
$s2[23][1] =  1.049 ;
$s2[23][2] =  0.051 ;
$s2[23][3] =  1.635 ;
$s2[23][4] =  6.213 ;
$s2[23][5] =  6.213 ;
$s2[23][8] =  0.032 ;
$s2[23][15] = -0.048 ;
$s2[23][16] = -7.508 ;
$s2[23][17] =  6.213 ;
$s2[23][23] =  6.927 ;
$s2[23][25] = -1.564 ;
$s2[23][38] = -1.295 ;
$s2[25][1] =  0.282 ;
$s2[25][2] =  0.148 ;
$s2[25][3] =  0.677 ;
$s2[25][4] = -0.173 ;
$s2[25][5] =  0.153 ;
$s2[25][8] = -0.239 ;
$s2[25][9] = -8.959 ;
$s2[25][14] =  0.810 ;
$s2[25][15] =  1.810 ;
$s2[25][16] =  5.645 ;
$s2[25][17] = -1.020 ;
$s2[25][19] =  6.213 ;
$s2[25][23] = -1.564 ;
$s2[25][25] =  0.129 ;
$s2[25][26] = -0.361 ;
$s2[25][31] =  7.641 ;
$s2[25][36] =  6.927 ;
$s2[25][38] =  7.066 ;
$s2[25][41] = -6.090 ;
$s2[25][44] =  5.704 ;
$s2[25][45] = -7.227 ;
$s2[25][46] = -6.152 ;
$s2[26][1] = -0.002 ;
$s2[26][2] =  0.187 ;
$s2[26][3] = -0.357 ;
$s2[26][4] = -8.317 ;
$s2[26][5] = -0.791 ;
$s2[26][7] =  6.927 ;
$s2[26][8] = -0.497 ;
$s2[26][9] = -8.658 ;
$s2[26][15] =  0.596 ;
$s2[26][17] = -0.656 ;
$s2[26][25] = -0.361 ;
$s2[26][26] = -0.287 ;
$s2[26][28] = -1.847 ;
$s2[26][31] = -2.043 ;
$s2[26][36] =  6.213 ;
$s2[26][38] =  0.805 ;
$s2[26][44] = -9.054 ;
$s2[26][45] = -9.208 ;
$s2[26][46] = -6.415 ;
$s2[27][1] = -0.774 ;
$s2[27][2] =  1.009 ;
$s2[27][3] =  0.551 ;
$s2[27][8] =  7.112 ;
$s2[27][16] =  9.607 ;
$s2[28][1] = -9.684 ;
$s2[28][2] =  9.528 ;
$s2[28][3] =  6.927 ;
$s2[28][8] =  6.648 ;
$s2[28][15] =  6.213 ;
$s2[28][26] = -1.847 ;
$s2[28][28] = -0.622 ;
$s2[29][1] =  0.223 ;
$s2[29][2] = -0.174 ;
$s2[29][3] =  7.871 ;
$s2[29][9] = -7.508 ;
$s2[29][36] =  6.213 ;
$s2[31][1] = -2.032 ;
$s2[31][15] =  6.213 ;
$s2[31][25] =  7.641 ;
$s2[31][26] = -2.043 ;
$s2[32][1] = -3.735 ;
$s2[32][2] = -8.233 ;
$s2[32][8] = -8.554 ;
$s2[32][15] =  6.213 ;
$s2[32][32] = -9.383 ;
$s2[34][1] = -3.023 ;
$s2[34][2] = -8.559 ;
$s2[34][3] = -9.609 ;
$s2[34][4] = -8.233 ;
$s2[34][8] = -7.866 ;
$s2[34][16] = -8.233 ;
$s2[34][34] =-10.642 ;
$s2[34][38] = -8.233 ;
$s2[34][41] = -7.508 ;
$s2[34][45] = -7.508 ;
$s2[35][1] = -8.959 ;
$s2[35][8] =-10.409 ;
$s2[36][1] =  0.239 ;
$s2[36][2] = -1.295 ;
$s2[36][3] =  6.927 ;
$s2[36][4] = -7.508 ;
$s2[36][8] = -1.872 ;
$s2[36][15] =  6.213 ;
$s2[36][25] =  6.927 ;
$s2[36][26] =  6.213 ;
$s2[36][29] =  6.213 ;
$s2[36][36] = -2.456 ;
$s2[36][45] = -7.508 ;
$s2[38][1] =  0.463 ;
$s2[38][2] =  0.657 ;
$s2[38][3] =  1.034 ;
$s2[38][4] =  6.213 ;
$s2[38][5] =  6.213 ;
$s2[38][8] =  0.680 ;
$s2[38][14] =  7.467 ;
$s2[38][15] = -0.543 ;
$s2[38][16] = -1.306 ;
$s2[38][17] = -1.726 ;
$s2[38][23] = -1.295 ;
$s2[38][25] =  7.066 ;
$s2[38][26] =  0.805 ;
$s2[38][34] = -8.233 ;
$s2[38][38] = -1.418 ;
$s2[38][41] = -0.362 ;
$s2[38][42] =  6.213 ;
$s2[38][44] =  8.476 ;
$s2[38][45] = -1.313 ;
$s2[39][1] =  6.213 ;
$s2[40][1] = -1.748 ;
$s2[40][2] = -0.174 ;
$s2[40][3] =  6.927 ;
$s2[40][8] =  0.482 ;
$s2[40][15] =  8.217 ;
$s2[41][1] =  0.831 ;
$s2[41][2] =  1.188 ;
$s2[41][3] =  0.979 ;
$s2[41][4] = -7.508 ;
$s2[41][8] =  0.066 ;
$s2[41][14] = -8.233 ;
$s2[41][15] =  1.513 ;
$s2[41][17] =  6.213 ;
$s2[41][25] = -6.090 ;
$s2[41][34] = -7.508 ;
$s2[41][38] = -0.362 ;
$s2[41][41] = -8.233 ;
$s2[41][42] =  7.641 ;
$s2[41][44] = -1.313 ;
$s2[42][1] =  0.083 ;
$s2[42][2] =  0.577 ;
$s2[42][3] = -1.088 ;
$s2[42][4] = -7.508 ;
$s2[42][15] = -0.174 ;
$s2[42][16] =  6.213 ;
$s2[42][17] =  6.213 ;
$s2[42][38] =  6.213 ;
$s2[42][41] =  7.641 ;
$s2[42][45] =  7.345 ;
$s2[42][46] =  6.927 ;
$s2[42][47] =  6.213 ;
$s2[44][1] =  1.223 ;
$s2[44][2] = -0.034 ;
$s2[44][3] = -1.901 ;
$s2[44][5] =  6.927 ;
$s2[44][8] =  0.088 ;
$s2[44][25] =  5.704 ;
$s2[44][26] = -9.054 ;
$s2[44][38] =  8.476 ;
$s2[44][41] = -1.313 ;
$s2[44][44] =  0.810 ;
$s2[45][1] = -0.921 ;
$s2[45][2] = -0.744 ;
$s2[45][3] = -0.801 ;
$s2[45][5] =  7.641 ;
$s2[45][8] =  0.770 ;
$s2[45][15] = -0.415 ;
$s2[45][17] = -1.122 ;
$s2[45][25] = -7.227 ;
$s2[45][26] = -9.208 ;
$s2[45][34] = -7.508 ;
$s2[45][36] = -7.508 ;
$s2[45][38] = -1.313 ;
$s2[45][42] =  7.345 ;
$s2[45][45] = -0.866 ;
$s2[45][46] = -8.959 ;
$s2[46][1] =  0.215 ;
$s2[46][2] =  0.784 ;
$s2[46][3] = -9.742 ;
$s2[46][4] = -7.508 ;
$s2[46][8] =  0.503 ;
$s2[46][15] =  6.213 ;
$s2[46][25] = -6.152 ;
$s2[46][26] = -6.415 ;
$s2[46][42] =  6.927 ;
$s2[46][45] = -8.959 ;
$s2[46][46] =-10.269 ;
$s2[47][1] =  0.100 ;
$s2[47][2] = -1.731 ;
$s2[47][3] = -8.233 ;
$s2[47][4] =  6.213 ;
$s2[47][8] = -0.403 ;
$s2[47][42] =  6.213 ;

$s3[1][1] = -0.305 ;
$s3[1][2] =  0.145 ;
$s3[1][3] =  0.335 ;
$s3[1][4] =  0.414 ;
$s3[1][5] = -0.487 ;
$s3[1][8] =  0.372 ;
$s3[1][9] = -1.346 ;
$s3[1][10] = -2.165 ;
$s3[1][11] = -8.933 ;
$s3[1][14] = -0.047 ;
$s3[1][15] =  0.857 ;
$s3[1][16] =  0.030 ;
$s3[1][17] =  0.750 ;
$s3[1][18] = -8.933 ;
$s3[1][19] = -0.496 ;
$s3[1][20] =  7.746 ;
$s3[1][22] = -8.371 ;
$s3[1][23] =  0.800 ;
$s3[1][25] =  0.333 ;
$s3[1][26] =  0.343 ;
$s3[1][27] = -2.591 ;
$s3[1][28] =  5.433 ;
$s3[1][29] =  0.269 ;
$s3[1][31] = -0.942 ;
$s3[1][32] = -2.513 ;
$s3[1][34] = -3.068 ;
$s3[1][35] = -7.925 ;
$s3[1][36] =  0.191 ;
$s3[1][38] = -0.017 ;
$s3[1][40] = -1.949 ;
$s3[1][41] =  0.248 ;
$s3[1][42] = -0.238 ;
$s3[1][44] =  7.432 ;
$s3[1][45] = -0.718 ;
$s3[1][46] =  0.503 ;
$s3[1][47] =  0.555 ;
$s3[2][1] =  0.145 ;
$s3[2][2] = -0.215 ;
$s3[2][3] =  0.228 ;
$s3[2][4] = -0.667 ;
$s3[2][5] =  0.531 ;
$s3[2][8] =  0.602 ;
$s3[2][9] =-10.141 ;
$s3[2][10] =  6.595 ;
$s3[2][11] = -7.925 ;
$s3[2][14] =  0.157 ;
$s3[2][15] = -0.092 ;
$s3[2][16] =  0.844 ;
$s3[2][17] =  0.268 ;
$s3[2][19] = -0.591 ;
$s3[2][20] =  7.746 ;
$s3[2][23] = -0.131 ;
$s3[2][25] =  0.051 ;
$s3[2][26] =  0.208 ;
$s3[2][27] =  0.733 ;
$s3[2][28] =  9.199 ;
$s3[2][29] = -0.054 ;
$s3[2][31] =  5.869 ;
$s3[2][32] = -7.925 ;
$s3[2][34] = -7.706 ;
$s3[2][36] =  5.869 ;
$s3[2][38] =  0.508 ;
$s3[2][40] = -1.776 ;
$s3[2][41] =  0.615 ;
$s3[2][42] =  0.666 ;
$s3[2][44] =  8.460 ;
$s3[2][45] = -0.710 ;
$s3[2][46] =  7.321 ;
$s3[2][47] = -1.748 ;
$s3[3][1] =  0.335 ;
$s3[3][2] =  0.228 ;
$s3[3][3] =  0.075 ;
$s3[3][4] = -0.381 ;
$s3[3][5] = -0.058 ;
$s3[3][8] = -0.098 ;
$s3[3][9] = -9.896 ;
$s3[3][14] = -0.491 ;
$s3[3][15] =  0.003 ;
$s3[3][16] = -0.925 ;
$s3[3][17] =  0.575 ;
$s3[3][19] =  1.118 ;
$s3[3][23] =  8.681 ;
$s3[3][25] =  0.129 ;
$s3[3][26] =  0.047 ;
$s3[3][27] =  7.321 ;
$s3[3][28] = -9.896 ;
$s3[3][29] =  7.555 ;
$s3[3][32] =  5.869 ;
$s3[3][36] = -1.330 ;
$s3[3][38] =  1.602 ;
$s3[3][39] =  5.869 ;
$s3[3][40] =  0.158 ;
$s3[3][41] =  1.355 ;
$s3[3][42] = -0.132 ;
$s3[3][44] = -2.478 ;
$s3[3][45] = -1.380 ;
$s3[3][46] = -1.672 ;
$s3[4][1] =  0.414 ;
$s3[4][2] = -0.667 ;
$s3[4][3] = -0.381 ;
$s3[4][4] = -9.695 ;
$s3[4][5] =  7.321 ;
$s3[4][8] =  0.768 ;
$s3[4][14] = -0.781 ;
$s3[4][15] = -1.406 ;
$s3[4][16] = -7.925 ;
$s3[4][17] = -1.812 ;
$s3[4][25] =  0.620 ;
$s3[4][26] = -1.531 ;
$s3[4][36] = -7.163 ;
$s3[4][38] =  7.020 ;
$s3[4][41] = -7.163 ;
$s3[4][44] = -0.143 ;
$s3[4][45] = -8.933 ;
$s3[4][46] = -7.163 ;
$s3[5][1] = -0.487 ;
$s3[5][2] =  0.531 ;
$s3[5][3] = -0.058 ;
$s3[5][4] =  7.321 ;
$s3[5][5] = -2.564 ;
$s3[5][8] =  1.576 ;
$s3[5][14] =  1.860 ;
$s3[5][15] =  9.536 ;
$s3[5][17] =  7.321 ;
$s3[5][23] =  7.020 ;
$s3[5][25] =  0.316 ;
$s3[5][26] = -0.690 ;
$s3[5][38] =  6.595 ;
$s3[5][44] =  7.020 ;
$s3[5][45] =  7.746 ;
$s3[7][7] =  7.321 ;
$s3[7][8] =  8.281 ;
$s3[7][26] =  6.595 ;
$s3[8][1] =  0.372 ;
$s3[8][2] =  0.602 ;
$s3[8][3] = -0.098 ;
$s3[8][4] =  0.768 ;
$s3[8][5] =  1.576 ;
$s3[8][7] =  8.281 ;
$s3[8][8] = -0.127 ;
$s3[8][9] =  8.047 ;
$s3[8][14] = -0.362 ;
$s3[8][15] =  0.071 ;
$s3[8][16] = -0.069 ;
$s3[8][17] = -0.040 ;
$s3[8][19] =  0.467 ;
$s3[8][23] = -0.058 ;
$s3[8][25] = -0.315 ;
$s3[8][26] = -0.322 ;
$s3[8][27] =  7.275 ;
$s3[8][28] =  7.763 ;
$s3[8][32] = -8.578 ;
$s3[8][34] = -7.620 ;
$s3[8][35] =-10.342 ;
$s3[8][36] = -1.784 ;
$s3[8][38] =  0.327 ;
$s3[8][40] =  0.877 ;
$s3[8][41] = -0.136 ;
$s3[8][42] =  6.936 ;
$s3[8][44] =  0.067 ;
$s3[8][45] =  0.486 ;
$s3[8][46] =  0.516 ;
$s3[8][47] = -0.383 ;
$s3[9][1] = -1.346 ;
$s3[9][2] =-10.141 ;
$s3[9][3] = -9.896 ;
$s3[9][8] =  8.047 ;
$s3[9][9] = -1.328 ;
$s3[9][11] = -7.925 ;
$s3[9][14] =  0.122 ;
$s3[9][25] = -9.303 ;
$s3[9][26] = -9.450 ;
$s3[9][29] = -8.371 ;
$s3[10][1] = -2.165 ;
$s3[10][2] =  6.595 ;
$s3[10][10] = -3.617 ;
$s3[11][1] = -8.933 ;
$s3[11][2] = -7.925 ;
$s3[11][9] = -7.925 ;
$s3[14][1] = -0.047 ;
$s3[14][2] =  0.157 ;
$s3[14][3] = -0.491 ;
$s3[14][4] = -0.781 ;
$s3[14][5] =  1.860 ;
$s3[14][8] = -0.362 ;
$s3[14][9] =  0.122 ;
$s3[14][14] = -0.540 ;
$s3[14][15] =  0.401 ;
$s3[14][16] = -2.056 ;
$s3[14][17] =  6.954 ;
$s3[14][23] = -0.568 ;
$s3[14][25] =  0.432 ;
$s3[14][26] = -0.687 ;
$s3[14][27] =  6.595 ;
$s3[14][29] = -0.568 ;
$s3[14][32] = -9.133 ;
$s3[14][38] = -1.623 ;
$s3[14][40] = -7.163 ;
$s3[14][41] =  0.392 ;
$s3[14][42] = -0.604 ;
$s3[14][44] =  5.334 ;
$s3[14][45] = -3.579 ;
$s3[15][1] =  0.857 ;
$s3[15][2] = -0.092 ;
$s3[15][3] =  0.003 ;
$s3[15][4] = -1.406 ;
$s3[15][5] =  9.536 ;
$s3[15][8] =  0.071 ;
$s3[15][14] =  0.401 ;
$s3[15][15] =  0.359 ;
$s3[15][16] = -0.097 ;
$s3[15][17] = -0.329 ;
$s3[15][19] =  8.171 ;
$s3[15][20] =  7.020 ;
$s3[15][23] =  2.373 ;
$s3[15][25] =  1.136 ;
$s3[15][26] =  1.237 ;
$s3[15][27] =  1.082 ;
$s3[15][28] =  9.199 ;
$s3[15][29] =  5.869 ;
$s3[15][34] = -7.097 ;
$s3[15][36] =  7.020 ;
$s3[15][38] =  0.264 ;
$s3[15][40] = -0.781 ;
$s3[15][41] =  1.315 ;
$s3[15][42] =  0.631 ;
$s3[15][44] = -0.659 ;
$s3[15][45] =  0.345 ;
$s3[15][46] = -1.225 ;
$s3[15][47] = -8.688 ;
$s3[16][1] =  0.030 ;
$s3[16][2] =  0.844 ;
$s3[16][3] = -0.925 ;
$s3[16][4] = -7.925 ;
$s3[16][8] = -0.069 ;
$s3[16][14] = -2.056 ;
$s3[16][15] = -0.097 ;
$s3[16][16] =  7.321 ;
$s3[16][17] = -3.205 ;
$s3[16][23] = -1.294 ;
$s3[16][25] = -0.811 ;
$s3[16][26] = -0.035 ;
$s3[16][27] =  5.869 ;
$s3[16][34] = -7.163 ;
$s3[16][38] =  5.869 ;
$s3[16][42] = -7.163 ;
$s3[16][45] = -2.819 ;
$s3[16][46] = -7.163 ;
$s3[17][1] =  0.750 ;
$s3[17][2] =  0.268 ;
$s3[17][3] =  0.575 ;
$s3[17][4] = -1.812 ;
$s3[17][5] =  7.321 ;
$s3[17][8] = -0.040 ;
$s3[17][14] =  6.954 ;
$s3[17][15] = -0.329 ;
$s3[17][16] = -3.205 ;
$s3[17][17] = -0.734 ;
$s3[17][19] = -1.294 ;
$s3[17][25] =  0.020 ;
$s3[17][26] = -0.075 ;
$s3[17][27] =  7.020 ;
$s3[17][36] =  5.869 ;
$s3[17][38] =  1.019 ;
$s3[17][40] =  7.020 ;
$s3[17][41] =  1.734 ;
$s3[17][44] =  8.472 ;
$s3[17][45] = -0.318 ;
$s3[17][46] = -2.093 ;
$s3[18][1] = -8.933 ;
$s3[18][41] = -7.163 ;
$s3[19][1] = -0.496 ;
$s3[19][2] = -0.591 ;
$s3[19][3] =  1.118 ;
$s3[19][8] =  0.467 ;
$s3[19][15] =  8.171 ;
$s3[19][17] = -1.294 ;
$s3[19][25] =  0.583 ;
$s3[19][26] =  0.547 ;
$s3[19][44] =  5.869 ;
$s3[19][45] =  5.869 ;
$s3[20][1] =  7.746 ;
$s3[20][2] =  7.746 ;
$s3[20][15] =  7.020 ;
$s3[22][1] = -8.371 ;
$s3[23][1] =  0.800 ;
$s3[23][2] = -0.131 ;
$s3[23][3] =  8.681 ;
$s3[23][5] =  7.020 ;
$s3[23][8] = -0.058 ;
$s3[23][14] = -0.568 ;
$s3[23][15] =  2.373 ;
$s3[23][16] = -1.294 ;
$s3[23][23] = -0.604 ;
$s3[23][25] = -1.084 ;
$s3[23][26] =  0.989 ;
$s3[23][27] =  7.555 ;
$s3[23][36] = -7.163 ;
$s3[23][38] = -7.163 ;
$s3[23][41] =  7.321 ;
$s3[23][44] =  5.869 ;
$s3[25][1] =  0.333 ;
$s3[25][2] =  0.051 ;
$s3[25][3] =  0.129 ;
$s3[25][4] =  0.620 ;
$s3[25][5] =  0.316 ;
$s3[25][8] = -0.315 ;
$s3[25][9] = -9.303 ;
$s3[25][14] =  0.432 ;
$s3[25][15] =  1.136 ;
$s3[25][16] = -0.811 ;
$s3[25][17] =  0.020 ;
$s3[25][19] =  0.583 ;
$s3[25][23] = -1.084 ;
$s3[25][25] =  0.219 ;
$s3[25][26] =  0.757 ;
$s3[25][27] =  5.869 ;
$s3[25][29] = -0.018 ;
$s3[25][34] = -9.092 ;
$s3[25][36] = -0.942 ;
$s3[25][38] =  7.586 ;
$s3[25][40] =  6.595 ;
$s3[25][41] =  0.592 ;
$s3[25][42] = -0.179 ;
$s3[25][44] = -0.322 ;
$s3[25][45] = -0.364 ;
$s3[25][46] = -1.077 ;
$s3[25][47] = -1.133 ;
$s3[26][1] =  0.343 ;
$s3[26][2] =  0.208 ;
$s3[26][3] =  0.047 ;
$s3[26][4] = -1.531 ;
$s3[26][5] = -0.690 ;
$s3[26][7] =  6.595 ;
$s3[26][8] = -0.322 ;
$s3[26][9] = -9.450 ;
$s3[26][14] = -0.687 ;
$s3[26][15] =  1.237 ;
$s3[26][16] = -0.035 ;
$s3[26][17] = -0.075 ;
$s3[26][19] =  0.547 ;
$s3[26][23] =  0.989 ;
$s3[26][25] =  0.757 ;
$s3[26][26] = -0.185 ;
$s3[26][27] =  2.996 ;
$s3[26][28] =  6.595 ;
$s3[26][29] =  7.746 ;
$s3[26][34] = -9.273 ;
$s3[26][36] = -0.604 ;
$s3[26][38] =  0.139 ;
$s3[26][40] =  6.595 ;
$s3[26][41] =  0.022 ;
$s3[26][42] =  0.018 ;
$s3[26][44] = -1.977 ;
$s3[26][45] = -0.747 ;
$s3[26][46] = -9.192 ;
$s3[26][47] = -7.925 ;
$s3[27][1] = -2.591 ;
$s3[27][2] =  0.733 ;
$s3[27][3] =  7.321 ;
$s3[27][8] =  7.275 ;
$s3[27][14] =  6.595 ;
$s3[27][15] =  1.082 ;
$s3[27][16] =  5.869 ;
$s3[27][17] =  7.020 ;
$s3[27][23] =  7.555 ;
$s3[27][25] =  5.869 ;
$s3[27][26] =  2.996 ;
$s3[27][38] =  6.595 ;
$s3[28][1] =  5.433 ;
$s3[28][2] =  9.199 ;
$s3[28][3] = -9.896 ;
$s3[28][8] =  7.763 ;
$s3[28][15] =  9.199 ;
$s3[28][26] =  6.595 ;
$s3[28][44] = -8.371 ;
$s3[29][1] =  0.269 ;
$s3[29][2] = -0.054 ;
$s3[29][3] =  7.555 ;
$s3[29][9] = -8.371 ;
$s3[29][14] = -0.568 ;
$s3[29][15] =  5.869 ;
$s3[29][25] = -0.018 ;
$s3[29][26] =  7.746 ;
$s3[29][29] =  0.173 ;
$s3[29][36] =  5.869 ;
$s3[31][1] = -0.942 ;
$s3[31][2] =  5.869 ;
$s3[31][36] =  5.869 ;
$s3[32][1] = -2.513 ;
$s3[32][2] = -7.925 ;
$s3[32][3] =  5.869 ;
$s3[32][8] = -8.578 ;
$s3[32][14] = -9.133 ;
$s3[34][1] = -3.068 ;
$s3[34][2] = -7.706 ;
$s3[34][8] = -7.620 ;
$s3[34][15] = -7.097 ;
$s3[34][16] = -7.163 ;
$s3[34][25] = -9.092 ;
$s3[34][26] = -9.273 ;
$s3[34][44] = -7.925 ;
$s3[35][1] = -7.925 ;
$s3[35][8] =-10.342 ;
$s3[36][1] =  0.191 ;
$s3[36][2] =  5.869 ;
$s3[36][3] = -1.330 ;
$s3[36][4] = -7.163 ;
$s3[36][8] = -1.784 ;
$s3[36][15] =  7.020 ;
$s3[36][17] =  5.869 ;
$s3[36][23] = -7.163 ;
$s3[36][25] = -0.942 ;
$s3[36][26] = -0.604 ;
$s3[36][29] =  5.869 ;
$s3[36][31] =  5.869 ;
$s3[36][38] = -9.450 ;
$s3[36][45] =  7.020 ;
$s3[38][1] = -0.017 ;
$s3[38][2] =  0.508 ;
$s3[38][3] =  1.602 ;
$s3[38][4] =  7.020 ;
$s3[38][5] =  6.595 ;
$s3[38][8] =  0.327 ;
$s3[38][14] = -1.623 ;
$s3[38][15] =  0.264 ;
$s3[38][16] =  5.869 ;
$s3[38][17] =  1.019 ;
$s3[38][23] = -7.163 ;
$s3[38][25] =  7.586 ;
$s3[38][26] =  0.139 ;
$s3[38][27] =  6.595 ;
$s3[38][36] = -9.450 ;
$s3[38][38] =  8.472 ;
$s3[38][40] =  5.869 ;
$s3[38][41] =  7.020 ;
$s3[38][42] =  6.595 ;
$s3[38][44] =  7.907 ;
$s3[38][45] =  0.158 ;
$s3[38][47] =  5.869 ;
$s3[39][3] =  5.869 ;
$s3[40][1] = -1.949 ;
$s3[40][2] = -1.776 ;
$s3[40][3] =  0.158 ;
$s3[40][8] =  0.877 ;
$s3[40][14] = -7.163 ;
$s3[40][15] = -0.781 ;
$s3[40][17] =  7.020 ;
$s3[40][25] =  6.595 ;
$s3[40][26] =  6.595 ;
$s3[40][38] =  5.869 ;
$s3[41][1] =  0.248 ;
$s3[41][2] =  0.615 ;
$s3[41][3] =  1.355 ;
$s3[41][4] = -7.163 ;
$s3[41][8] = -0.136 ;
$s3[41][14] =  0.392 ;
$s3[41][15] =  1.315 ;
$s3[41][17] =  1.734 ;
$s3[41][18] = -7.163 ;
$s3[41][23] =  7.321 ;
$s3[41][25] =  0.592 ;
$s3[41][26] =  0.022 ;
$s3[41][38] =  7.020 ;
$s3[41][41] = -2.538 ;
$s3[41][42] =  6.595 ;
$s3[41][44] =  5.869 ;
$s3[41][45] =  9.532 ;
$s3[42][1] = -0.238 ;
$s3[42][2] =  0.666 ;
$s3[42][3] = -0.132 ;
$s3[42][8] =  6.936 ;
$s3[42][14] = -0.604 ;
$s3[42][15] =  0.631 ;
$s3[42][16] = -7.163 ;
$s3[42][25] = -0.179 ;
$s3[42][26] =  0.018 ;
$s3[42][38] =  6.595 ;
$s3[42][41] =  6.595 ;
$s3[42][45] =  5.869 ;
$s3[44][1] =  7.432 ;
$s3[44][2] =  8.460 ;
$s3[44][3] = -2.478 ;
$s3[44][4] = -0.143 ;
$s3[44][5] =  7.020 ;
$s3[44][8] =  0.067 ;
$s3[44][14] =  5.334 ;
$s3[44][15] = -0.659 ;
$s3[44][17] =  8.472 ;
$s3[44][19] =  5.869 ;
$s3[44][23] =  5.869 ;
$s3[44][25] = -0.322 ;
$s3[44][26] = -1.977 ;
$s3[44][28] = -8.371 ;
$s3[44][34] = -7.925 ;
$s3[44][38] =  7.907 ;
$s3[44][41] =  5.869 ;
$s3[44][44] = -3.278 ;
$s3[44][45] =  5.869 ;
$s3[45][1] = -0.718 ;
$s3[45][2] = -0.710 ;
$s3[45][3] = -1.380 ;
$s3[45][4] = -8.933 ;
$s3[45][5] =  7.746 ;
$s3[45][8] =  0.486 ;
$s3[45][14] = -3.579 ;
$s3[45][15] =  0.345 ;
$s3[45][16] = -2.819 ;
$s3[45][17] = -0.318 ;
$s3[45][19] =  5.869 ;
$s3[45][25] = -0.364 ;
$s3[45][26] = -0.747 ;
$s3[45][36] =  7.020 ;
$s3[45][38] =  0.158 ;
$s3[45][41] =  9.532 ;
$s3[45][42] =  5.869 ;
$s3[45][44] =  5.869 ;
$s3[45][45] = -0.894 ;
$s3[46][1] =  0.503 ;
$s3[46][2] =  7.321 ;
$s3[46][3] = -1.672 ;
$s3[46][4] = -7.163 ;
$s3[46][8] =  0.516 ;
$s3[46][15] = -1.225 ;
$s3[46][16] = -7.163 ;
$s3[46][17] = -2.093 ;
$s3[46][25] = -1.077 ;
$s3[46][26] = -9.192 ;
$s3[46][46] = -9.450 ;
$s3[47][1] =  0.555 ;
$s3[47][2] = -1.748 ;
$s3[47][8] = -0.383 ;
$s3[47][15] = -8.688 ;
$s3[47][25] = -1.133 ;
$s3[47][26] = -7.925 ;
$s3[47][38] =  5.869 ;

$s4[1][1] = -0.410 ;
$s4[1][2] =  0.081 ;
$s4[1][3] =  0.371 ;
$s4[1][4] =  0.512 ;
$s4[1][5] = -0.513 ;
$s4[1][8] =  0.414 ;
$s4[1][9] = -1.388 ;
$s4[1][10] = -2.507 ;
$s4[1][11] = -7.507 ;
$s4[1][14] =  0.035 ;
$s4[1][15] =  0.879 ;
$s4[1][16] =  1.137 ;
$s4[1][17] =  0.694 ;
$s4[1][18] = -9.386 ;
$s4[1][19] = -0.156 ;
$s4[1][20] =  6.672 ;
$s4[1][22] = -5.980 ;
$s4[1][23] =  1.452 ;
$s4[1][25] =  0.291 ;
$s4[1][26] =  0.441 ;
$s4[1][27] = -1.173 ;
$s4[1][28] = -0.526 ;
$s4[1][29] =  0.418 ;
$s4[1][31] = -0.705 ;
$s4[1][32] = -2.860 ;
$s4[1][34] = -2.211 ;
$s4[1][35] = -8.271 ;
$s4[1][36] =  0.147 ;
$s4[1][38] =  0.568 ;
$s4[1][39] =  5.960 ;
$s4[1][40] = -1.520 ;
$s4[1][41] =  0.940 ;
$s4[1][42] =  0.008 ;
$s4[1][44] =  1.048 ;
$s4[1][45] =  0.041 ;
$s4[1][46] =  0.246 ;
$s4[1][47] =  1.242 ;
$s4[2][1] =  0.041 ;
$s4[2][2] = -0.291 ;
$s4[2][3] =  0.079 ;
$s4[2][4] = -0.402 ;
$s4[2][5] =  0.534 ;
$s4[2][8] =  0.659 ;
$s4[2][9] =-10.222 ;
$s4[2][10] =  6.672 ;
$s4[2][11] = -7.507 ;
$s4[2][14] =  0.074 ;
$s4[2][15] =  0.170 ;
$s4[2][16] =  1.509 ;
$s4[2][17] =  0.308 ;
$s4[2][19] = -0.516 ;
$s4[2][20] =  6.830 ;
$s4[2][23] = -0.120 ;
$s4[2][25] =  0.115 ;
$s4[2][26] =  0.452 ;
$s4[2][27] =  0.833 ;
$s4[2][28] =  5.960 ;
$s4[2][29] = -0.268 ;
$s4[2][31] =  4.831 ;
$s4[2][32] = -7.954 ;
$s4[2][34] = -8.435 ;
$s4[2][36] = -0.518 ;
$s4[2][38] =  0.167 ;
$s4[2][40] = -0.211 ;
$s4[2][41] =  1.140 ;
$s4[2][42] =  0.528 ;
$s4[2][44] =  1.179 ;
$s4[2][45] =  0.409 ;
$s4[2][46] =  7.680 ;
$s4[2][47] = -1.792 ;
$s4[3][1] =  0.328 ;
$s4[3][2] =  0.079 ;
$s4[3][3] = -0.002 ;
$s4[3][4] = -0.866 ;
$s4[3][5] = -0.048 ;
$s4[3][7] =  5.543 ;
$s4[3][8] = -0.195 ;
$s4[3][9] = -9.833 ;
$s4[3][10] =  5.543 ;
$s4[3][14] = -0.211 ;
$s4[3][15] =  0.848 ;
$s4[3][16] = -0.393 ;
$s4[3][17] =  0.200 ;
$s4[3][19] =  0.731 ;
$s4[3][20] =  4.831 ;
$s4[3][23] =  1.975 ;
$s4[3][25] =  0.396 ;
$s4[3][26] =  0.227 ;
$s4[3][27] =  0.987 ;
$s4[3][28] =  8.178 ;
$s4[3][29] =  0.552 ;
$s4[3][34] = -9.403 ;
$s4[3][36] =  0.641 ;
$s4[3][38] =  1.720 ;
$s4[3][40] =  0.641 ;
$s4[3][41] =  1.241 ;
$s4[3][42] = -0.215 ;
$s4[3][44] =  0.397 ;
$s4[3][45] = -0.220 ;
$s4[3][46] = -0.292 ;
$s4[3][47] = -2.210 ;
$s4[4][1] =  0.237 ;
$s4[4][2] = -0.456 ;
$s4[4][3] = -1.176 ;
$s4[4][4] = -9.728 ;
$s4[4][5] =  7.466 ;
$s4[4][8] =  0.766 ;
$s4[4][14] = -0.661 ;
$s4[4][15] = -0.250 ;
$s4[4][16] = -6.743 ;
$s4[4][17] = -2.164 ;
$s4[4][23] =  4.831 ;
$s4[4][25] = -0.231 ;
$s4[4][26] = -0.501 ;
$s4[4][29] =  5.543 ;
$s4[4][34] = -7.954 ;
$s4[4][36] = -7.507 ;
$s4[4][38] =  4.831 ;
$s4[4][41] = -7.190 ;
$s4[4][42] = -5.980 ;
$s4[4][44] =  0.692 ;
$s4[4][45] = -2.210 ;
$s4[4][46] = -5.980 ;
$s4[5][1] = -0.443 ;
$s4[5][2] =  0.612 ;
$s4[5][3] =  0.404 ;
$s4[5][4] =  6.485 ;
$s4[5][5] = -3.467 ;
$s4[5][8] =  1.686 ;
$s4[5][14] =  1.900 ;
$s4[5][15] =  9.195 ;
$s4[5][17] =  0.799 ;
$s4[5][23] =  6.967 ;
$s4[5][25] =  0.605 ;
$s4[5][26] = -0.040 ;
$s4[5][38] =  6.255 ;
$s4[5][44] =  5.543 ;
$s4[5][45] =  7.197 ;
$s4[7][3] =  5.543 ;
$s4[7][7] =  5.543 ;
$s4[7][8] =  8.096 ;
$s4[8][1] =  0.367 ;
$s4[8][2] =  0.659 ;
$s4[8][3] = -0.195 ;
$s4[8][4] =  0.783 ;
$s4[8][5] =  1.782 ;
$s4[8][7] =  8.255 ;
$s4[8][8] = -0.134 ;
$s4[8][9] =  8.096 ;
$s4[8][14] =  0.012 ;
$s4[8][15] =  0.036 ;
$s4[8][16] = -0.133 ;
$s4[8][17] =  0.030 ;
$s4[8][19] =  0.549 ;
$s4[8][23] =  0.079 ;
$s4[8][25] = -0.217 ;
$s4[8][26] = -0.090 ;
$s4[8][27] =  7.295 ;
$s4[8][28] =  7.263 ;
$s4[8][31] =  4.831 ;
$s4[8][32] = -8.407 ;
$s4[8][34] = -7.905 ;
$s4[8][35] =-10.375 ;
$s4[8][36] = -1.973 ;
$s4[8][38] =  0.523 ;
$s4[8][40] =  0.612 ;
$s4[8][41] = -0.122 ;
$s4[8][42] =  7.054 ;
$s4[8][44] =  0.080 ;
$s4[8][45] =  0.560 ;
$s4[8][46] =  0.647 ;
$s4[8][47] = -0.085 ;
$s4[9][1] = -1.404 ;
$s4[9][2] = -9.988 ;
$s4[9][3] = -9.612 ;
$s4[9][8] =  7.909 ;
$s4[9][9] = -1.402 ;
$s4[9][11] = -7.507 ;
$s4[9][14] = -0.444 ;
$s4[9][25] = -9.102 ;
$s4[9][26] = -8.888 ;
$s4[9][29] = -7.954 ;
$s4[10][1] = -2.514 ;
$s4[10][2] =  6.672 ;
$s4[10][3] =  5.543 ;
$s4[10][10] = -3.397 ;
$s4[11][1] = -7.753 ;
$s4[11][2] = -7.507 ;
$s4[11][9] = -7.954 ;
$s4[11][11] = -6.743 ;
$s4[14][1] =  0.007 ;
$s4[14][2] =  0.061 ;
$s4[14][3] = -0.221 ;
$s4[14][4] = -0.369 ;
$s4[14][5] =  2.014 ;
$s4[14][8] =  0.010 ;
$s4[14][9] = -0.492 ;
$s4[14][14] = -1.350 ;
$s4[14][15] =  0.249 ;
$s4[14][16] =  0.275 ;
$s4[14][17] =  0.698 ;
$s4[14][19] =  4.831 ;
$s4[14][23] =  7.295 ;
$s4[14][25] =  0.700 ;
$s4[14][26] = -0.330 ;
$s4[14][29] = -1.200 ;
$s4[14][36] =  5.960 ;
$s4[14][38] =  0.720 ;
$s4[14][41] = -1.074 ;
$s4[14][42] =  0.505 ;
$s4[14][44] =  5.019 ;
$s4[14][45] =  6.369 ;
$s4[14][46] =  5.543 ;
$s4[15][1] =  0.847 ;
$s4[15][2] =  0.142 ;
$s4[15][3] =  0.857 ;
$s4[15][4] = -0.276 ;
$s4[15][5] =  9.660 ;
$s4[15][8] =  0.035 ;
$s4[15][14] =  0.211 ;
$s4[15][15] = -0.362 ;
$s4[15][16] = -0.085 ;
$s4[15][17] = -0.137 ;
$s4[15][19] = -0.041 ;
$s4[15][20] =  6.485 ;
$s4[15][23] = -0.052 ;
$s4[15][25] =  0.774 ;
$s4[15][26] =  0.176 ;
$s4[15][27] =  1.166 ;
$s4[15][28] =  4.831 ;
$s4[15][29] =  4.831 ;
$s4[15][32] =  5.543 ;
$s4[15][34] = -5.914 ;
$s4[15][36] =  6.830 ;
$s4[15][38] = -0.219 ;
$s4[15][39] =  4.831 ;
$s4[15][40] =  1.166 ;
$s4[15][41] =  0.503 ;
$s4[15][42] =  1.323 ;
$s4[15][44] = -0.863 ;
$s4[15][45] =  0.396 ;
$s4[15][46] =  6.255 ;
$s4[16][1] = -0.131 ;
$s4[16][2] =  1.086 ;
$s4[16][3] = -0.460 ;
$s4[16][4] = -6.743 ;
$s4[16][8] = -0.139 ;
$s4[16][14] =  0.275 ;
$s4[16][15] = -0.104 ;
$s4[16][16] = -9.035 ;
$s4[16][17] = -1.786 ;
$s4[16][20] =  4.831 ;
$s4[16][23] = -5.980 ;
$s4[16][25] =  1.954 ;
$s4[16][26] = -0.193 ;
$s4[16][27] =  8.255 ;
$s4[16][34] = -7.954 ;
$s4[16][38] = -0.705 ;
$s4[16][41] = -6.743 ;
$s4[16][42] =  4.831 ;
$s4[16][44] =  0.692 ;
$s4[16][45] = -3.123 ;
$s4[17][1] =  0.393 ;
$s4[17][2] =  0.309 ;
$s4[17][3] =  0.037 ;
$s4[17][4] = -2.164 ;
$s4[17][5] =  0.799 ;
$s4[17][8] = -0.052 ;
$s4[17][14] =  0.769 ;
$s4[17][15] = -0.233 ;
$s4[17][16] = -1.699 ;
$s4[17][17] = -0.191 ;
$s4[17][19] = -0.258 ;
$s4[17][20] =  5.960 ;
$s4[17][23] =  4.831 ;
$s4[17][25] = -0.425 ;
$s4[17][26] =  0.653 ;
$s4[17][27] =  4.831 ;
$s4[17][36] =  4.831 ;
$s4[17][38] = -0.282 ;
$s4[17][41] = -1.282 ;
$s4[17][42] =  4.831 ;
$s4[17][44] = -0.102 ;
$s4[17][45] = -1.278 ;
$s4[17][46] =  6.255 ;
$s4[17][47] =  5.960 ;
$s4[18][1] = -9.386 ;
$s4[18][26] = -6.743 ;
$s4[18][31] = -5.980 ;
$s4[18][41] = -5.980 ;
$s4[19][1] = -0.101 ;
$s4[19][2] = -0.516 ;
$s4[19][3] =  0.731 ;
$s4[19][8] =  0.272 ;
$s4[19][14] =  4.831 ;
$s4[19][15] = -0.041 ;
$s4[19][17] =  0.275 ;
$s4[19][19] = -6.743 ;
$s4[19][25] =  6.967 ;
$s4[19][26] =  1.217 ;
$s4[19][45] =  5.543 ;
$s4[20][1] =  6.672 ;
$s4[20][2] =  7.466 ;
$s4[20][3] =  4.831 ;
$s4[20][15] =  6.830 ;
$s4[20][16] =  4.831 ;
$s4[20][17] =  6.255 ;
$s4[20][25] =  5.960 ;
$s4[22][1] = -7.190 ;
$s4[23][1] =  1.254 ;
$s4[23][2] = -0.117 ;
$s4[23][3] =  1.810 ;
$s4[23][4] =  4.831 ;
$s4[23][5] =  7.197 ;
$s4[23][8] =  0.053 ;
$s4[23][14] =  7.295 ;
$s4[23][15] = -0.099 ;
$s4[23][16] = -6.743 ;
$s4[23][17] =  4.831 ;
$s4[23][23] =  0.799 ;
$s4[23][25] = -0.130 ;
$s4[23][26] =  1.039 ;
$s4[23][38] =  0.701 ;
$s4[23][41] =  4.831 ;
$s4[23][44] =  5.543 ;
$s4[23][45] =  6.672 ;
$s4[23][46] =  4.831 ;
$s4[25][1] =  0.236 ;
$s4[25][2] =  0.099 ;
$s4[25][3] =  0.432 ;
$s4[25][4] = -0.231 ;
$s4[25][5] =  0.447 ;
$s4[25][8] = -0.252 ;
$s4[25][9] = -9.527 ;
$s4[25][14] =  0.704 ;
$s4[25][15] =  0.913 ;
$s4[25][16] =  1.984 ;
$s4[25][17] = -0.407 ;
$s4[25][19] =  7.088 ;
$s4[25][20] =  5.960 ;
$s4[25][23] = -0.130 ;
$s4[25][25] =  0.090 ;
$s4[25][26] =  0.252 ;
$s4[25][27] = -0.488 ;
$s4[25][29] =  0.500 ;
$s4[25][31] =  5.960 ;
$s4[25][32] = -6.743 ;
$s4[25][36] =  0.007 ;
$s4[25][38] =  1.578 ;
$s4[25][40] = -0.835 ;
$s4[25][41] = -0.521 ;
$s4[25][42] = -1.230 ;
$s4[25][44] =  6.542 ;
$s4[25][45] =  0.834 ;
$s4[25][46] = -0.440 ;
$s4[25][47] =  4.831 ;
$s4[26][1] =  0.127 ;
$s4[26][2] =  0.230 ;
$s4[26][3] =  0.183 ;
$s4[26][4] = -1.466 ;
$s4[26][5] = -0.475 ;
$s4[26][7] =  6.255 ;
$s4[26][8] = -0.280 ;
$s4[26][9] = -9.435 ;
$s4[26][14] = -0.330 ;
$s4[26][15] =  0.438 ;
$s4[26][16] = -0.193 ;
$s4[26][17] =  0.447 ;
$s4[26][18] = -6.743 ;
$s4[26][19] =  1.217 ;
$s4[26][23] =  1.039 ;
$s4[26][25] = -0.208 ;
$s4[26][26] =  0.160 ;
$s4[26][27] =  0.505 ;
$s4[26][28] = -3.142 ;
$s4[26][29] =  6.672 ;
$s4[26][31] =  4.831 ;
$s4[26][32] =  4.831 ;
$s4[26][36] = -0.258 ;
$s4[26][38] =  1.688 ;
$s4[26][39] =  4.831 ;
$s4[26][40] =  0.275 ;
$s4[26][41] = -0.113 ;
$s4[26][42] =  0.235 ;
$s4[26][44] = -1.900 ;
$s4[26][45] =  0.265 ;
$s4[26][46] = -1.207 ;
$s4[27][1] = -1.089 ;
$s4[27][2] =  0.844 ;
$s4[27][3] =  0.799 ;
$s4[27][8] =  7.497 ;
$s4[27][15] =  1.166 ;
$s4[27][16] =  8.948 ;
$s4[27][17] =  4.831 ;
$s4[27][25] = -0.488 ;
$s4[27][26] =  0.505 ;
$s4[27][45] =  4.831 ;
$s4[28][1] = -1.174 ;
$s4[28][2] =  8.891 ;
$s4[28][3] =  8.255 ;
$s4[28][8] =  7.401 ;
$s4[28][15] =  5.543 ;
$s4[28][26] = -8.718 ;
$s4[28][28] = -8.718 ;
$s4[28][38] =  5.543 ;
$s4[29][1] =  0.225 ;
$s4[29][2] = -0.263 ;
$s4[29][3] =  0.937 ;
$s4[29][4] =  5.543 ;
$s4[29][9] = -8.271 ;
$s4[29][14] = -1.200 ;
$s4[29][15] =  4.831 ;
$s4[29][25] =  0.500 ;
$s4[29][26] =  6.672 ;
$s4[29][29] =  5.543 ;
$s4[29][36] =  5.960 ;
$s4[31][1] = -1.669 ;
$s4[31][2] =  4.831 ;
$s4[31][8] =  4.831 ;
$s4[31][15] =  4.831 ;
$s4[31][18] = -5.980 ;
$s4[31][25] =  6.830 ;
$s4[31][26] =  4.831 ;
$s4[32][1] = -3.110 ;
$s4[32][2] = -7.954 ;
$s4[32][8] = -8.407 ;
$s4[32][15] =  5.543 ;
$s4[32][25] = -6.743 ;
$s4[32][26] =  4.831 ;
$s4[32][32] = -9.165 ;
$s4[34][1] = -2.211 ;
$s4[34][2] = -8.052 ;
$s4[34][3] = -8.956 ;
$s4[34][4] = -6.743 ;
$s4[34][8] = -7.595 ;
$s4[34][15] = -5.914 ;
$s4[34][16] = -5.980 ;
$s4[34][34] =-10.661 ;
$s4[34][36] = -5.980 ;
$s4[34][38] = -7.190 ;
$s4[34][41] = -7.190 ;
$s4[34][44] = -8.622 ;
$s4[35][1] = -8.718 ;
$s4[35][8] =-10.375 ;
$s4[36][1] =  0.160 ;
$s4[36][2] = -0.258 ;
$s4[36][3] =  0.454 ;
$s4[36][4] = -5.980 ;
$s4[36][8] = -1.844 ;
$s4[36][14] =  7.088 ;
$s4[36][15] =  7.197 ;
$s4[36][17] =  4.831 ;
$s4[36][25] = -0.518 ;
$s4[36][26] = -1.200 ;
$s4[36][29] =  4.831 ;
$s4[36][34] = -5.980 ;
$s4[36][36] = -2.714 ;
$s4[36][45] =  6.672 ;
$s4[36][46] = -7.190 ;
$s4[38][1] =  0.477 ;
$s4[38][2] =  0.280 ;
$s4[38][3] =  1.666 ;
$s4[38][4] =  4.831 ;
$s4[38][5] =  6.672 ;
$s4[38][8] =  0.549 ;
$s4[38][14] =  1.059 ;
$s4[38][15] = -0.397 ;
$s4[38][16] = -0.360 ;
$s4[38][17] = -0.251 ;
$s4[38][23] =  0.701 ;
$s4[38][25] =  1.645 ;
$s4[38][26] =  1.768 ;
$s4[38][28] =  5.543 ;
$s4[38][34] = -8.124 ;
$s4[38][38] = -0.939 ;
$s4[38][40] =  4.831 ;
$s4[38][41] =  0.005 ;
$s4[38][42] =  5.543 ;
$s4[38][44] =  5.543 ;
$s4[38][45] =  0.723 ;
$s4[38][47] =  4.831 ;
$s4[39][1] =  5.960 ;
$s4[39][15] =  4.831 ;
$s4[39][26] =  4.831 ;
$s4[40][1] = -1.611 ;
$s4[40][2] = -0.211 ;
$s4[40][3] =  0.641 ;
$s4[40][8] =  0.612 ;
$s4[40][15] =  0.999 ;
$s4[40][25] = -0.835 ;
$s4[40][26] =  0.275 ;
$s4[40][38] =  4.831 ;
$s4[40][45] =  0.275 ;
$s4[41][1] =  0.940 ;
$s4[41][2] =  1.140 ;
$s4[41][3] =  1.341 ;
$s4[41][4] = -5.980 ;
$s4[41][8] = -0.130 ;
$s4[41][14] = -0.757 ;
$s4[41][15] =  0.423 ;
$s4[41][16] = -6.743 ;
$s4[41][17] = -1.699 ;
$s4[41][18] = -5.980 ;
$s4[41][23] =  4.831 ;
$s4[41][25] = -0.429 ;
$s4[41][26] = -0.113 ;
$s4[41][34] = -7.190 ;
$s4[41][38] =  0.552 ;
$s4[41][41] =  0.800 ;
$s4[41][42] =  6.672 ;
$s4[41][44] =  1.404 ;
$s4[41][45] =  1.486 ;
$s4[41][47] =  4.831 ;
$s4[42][1] =  0.004 ;
$s4[42][2] =  0.489 ;
$s4[42][3] = -0.484 ;
$s4[42][4] = -5.980 ;
$s4[42][8] =  7.054 ;
$s4[42][14] =  0.505 ;
$s4[42][15] =  0.918 ;
$s4[42][16] =  5.543 ;
$s4[42][17] =  4.831 ;
$s4[42][25] = -1.230 ;
$s4[42][26] =  0.235 ;
$s4[42][38] =  5.543 ;
$s4[42][41] =  7.542 ;
$s4[42][42] =  6.967 ;
$s4[42][47] =  5.543 ;
$s4[44][1] =  1.142 ;
$s4[44][2] =  0.595 ;
$s4[44][3] = -0.600 ;
$s4[44][4] =  0.692 ;
$s4[44][5] =  7.197 ;
$s4[44][8] =  0.023 ;
$s4[44][14] =  5.019 ;
$s4[44][15] = -0.863 ;
$s4[44][16] =  0.692 ;
$s4[44][17] = -0.102 ;
$s4[44][23] =  5.543 ;
$s4[44][25] =  6.654 ;
$s4[44][26] = -1.900 ;
$s4[44][34] = -8.622 ;
$s4[44][38] =  7.295 ;
$s4[44][41] = -0.558 ;
$s4[44][44] = -1.758 ;
$s4[44][45] =  0.454 ;
$s4[44][46] =  5.543 ;
$s4[45][1] = -0.653 ;
$s4[45][2] = -0.260 ;
$s4[45][3] = -0.572 ;
$s4[45][4] = -2.210 ;
$s4[45][5] =  8.178 ;
$s4[45][8] =  0.644 ;
$s4[45][14] =  6.369 ;
$s4[45][15] =  0.318 ;
$s4[45][16] = -3.123 ;
$s4[45][17] = -1.262 ;
$s4[45][19] =  5.543 ;
$s4[45][23] =  6.672 ;
$s4[45][25] =  0.316 ;
$s4[45][26] =  0.265 ;
$s4[45][27] =  4.831 ;
$s4[45][34] = -7.190 ;
$s4[45][36] = -0.835 ;
$s4[45][38] =  0.490 ;
$s4[45][40] =  0.275 ;
$s4[45][41] =  1.486 ;
$s4[45][42] =  5.960 ;
$s4[45][44] =  0.454 ;
$s4[45][45] = -0.125 ;
$s4[45][47] =  5.960 ;
$s4[46][1] =  0.177 ;
$s4[46][2] =  1.661 ;
$s4[46][3] = -2.436 ;
$s4[46][4] = -6.743 ;
$s4[46][8] =  0.510 ;
$s4[46][14] =  5.543 ;
$s4[46][15] =  6.672 ;
$s4[46][17] =  6.255 ;
$s4[46][23] =  4.831 ;
$s4[46][25] = -0.588 ;
$s4[46][26] = -1.207 ;
$s4[46][36] = -7.190 ;
$s4[46][42] =  5.543 ;
$s4[46][44] =  5.543 ;
$s4[46][46] = -3.345 ;
$s4[47][1] =  0.271 ;
$s4[47][2] = -2.014 ;
$s4[47][3] = -2.858 ;
$s4[47][4] =  4.831 ;
$s4[47][8] = -0.424 ;
$s4[47][17] =  5.960 ;
$s4[47][25] =  4.831 ;
$s4[47][38] =  4.831 ;
$s4[47][41] =  4.831 ;
$s4[47][42] =  4.831 ;
$s4[47][45] =  5.960 ;
$s4[47][47] = -3.175 ;

$s5[1][1] = -0.465 ;
$s5[1][2] =  0.197 ;
$s5[1][3] =  0.490 ;
$s5[1][4] =  1.132 ;
$s5[1][5] = -0.176 ;
$s5[1][8] =  0.638 ;
$s5[1][9] = -1.276 ;
$s5[1][10] = -2.410 ;
$s5[1][11] = -8.333 ;
$s5[1][14] = -0.048 ;
$s5[1][15] =  0.980 ;
$s5[1][16] =  1.222 ;
$s5[1][17] =  1.125 ;
$s5[1][18] = -8.896 ;
$s5[1][19] = -0.235 ;
$s5[1][20] =  6.539 ;
$s5[1][22] = -6.635 ;
$s5[1][23] =  1.592 ;
$s5[1][25] =  0.465 ;
$s5[1][26] =  0.902 ;
$s5[1][27] = -1.231 ;
$s5[1][28] =  8.276 ;
$s5[1][29] =  0.354 ;
$s5[1][31] =  0.673 ;
$s5[1][32] = -1.831 ;
$s5[1][34] = -2.650 ;
$s5[1][35] = -6.988 ;
$s5[1][36] =  0.168 ;
$s5[1][38] =  1.143 ;
$s5[1][40] = -1.176 ;
$s5[1][41] =  0.621 ;
$s5[1][42] =  0.239 ;
$s5[1][44] =  0.869 ;
$s5[1][45] =  0.712 ;
$s5[1][46] =  1.229 ;
$s5[1][47] =  0.837 ;
$s5[2][1] =  0.140 ;
$s5[2][2] = -0.340 ;
$s5[2][3] =  0.324 ;
$s5[2][4] =  0.189 ;
$s5[2][5] =  0.665 ;
$s5[2][8] =  0.786 ;
$s5[2][9] =-10.440 ;
$s5[2][10] =  6.908 ;
$s5[2][11] = -7.484 ;
$s5[2][14] =  0.297 ;
$s5[2][15] =  0.159 ;
$s5[2][16] =  1.741 ;
$s5[2][17] =  0.185 ;
$s5[2][19] = -0.887 ;
$s5[2][20] =  6.539 ;
$s5[2][23] = -0.250 ;
$s5[2][25] =  0.115 ;
$s5[2][26] =  0.531 ;
$s5[2][27] =  0.846 ;
$s5[2][28] =  5.246 ;
$s5[2][29] = -0.467 ;
$s5[2][32] = -7.484 ;
$s5[2][34] = -7.486 ;
$s5[2][36] =  6.192 ;
$s5[2][38] =  0.595 ;
$s5[2][40] = -0.158 ;
$s5[2][41] =  0.950 ;
$s5[2][42] =  1.015 ;
$s5[2][44] =  2.103 ;
$s5[2][45] =  0.315 ;
$s5[2][46] =  7.623 ;
$s5[2][47] = -1.311 ;
$s5[3][1] =  0.411 ;
$s5[3][2] =  0.324 ;
$s5[3][3] =  0.076 ;
$s5[3][4] =  0.990 ;
$s5[3][5] =  1.182 ;
$s5[3][8] = -0.021 ;
$s5[3][9] =-10.176 ;
$s5[3][10] =  6.380 ;
$s5[3][14] =  0.219 ;
$s5[3][15] =  0.095 ;
$s5[3][16] =  0.179 ;
$s5[3][17] =  0.805 ;
$s5[3][19] =  0.552 ;
$s5[3][20] =  4.531 ;
$s5[3][23] =  8.631 ;
$s5[3][25] =  0.019 ;
$s5[3][26] =  0.013 ;
$s5[3][27] =  2.644 ;
$s5[3][28] = -8.333 ;
$s5[3][29] =  1.375 ;
$s5[3][32] =  5.246 ;
$s5[3][36] = -0.576 ;
$s5[3][38] =  1.830 ;
$s5[3][39] =  5.665 ;
$s5[3][40] =  0.267 ;
$s5[3][41] =  1.101 ;
$s5[3][42] =  0.134 ;
$s5[3][44] =  1.311 ;
$s5[3][45] = -0.076 ;
$s5[3][46] = -0.676 ;
$s5[3][47] = -6.139 ;
$s5[4][1] =  0.756 ;
$s5[4][2] = -0.077 ;
$s5[4][3] =  0.670 ;
$s5[4][4] = -9.456 ;
$s5[4][5] =  7.096 ;
$s5[4][8] =  0.920 ;
$s5[4][14] = -0.267 ;
$s5[4][15] = -1.001 ;
$s5[4][16] = -7.484 ;
$s5[4][17] = -1.645 ;
$s5[4][25] =  0.780 ;
$s5[4][26] =  1.181 ;
$s5[4][29] =  5.962 ;
$s5[4][36] = -6.988 ;
$s5[4][38] =  6.799 ;
$s5[4][41] = -6.635 ;
$s5[4][42] = -5.289 ;
$s5[4][44] = -6.139 ;
$s5[4][45] = -1.026 ;
$s5[5][1] = -0.242 ;
$s5[5][2] =  0.852 ;
$s5[5][3] =  1.057 ;
$s5[5][4] =  5.962 ;
$s5[5][5] = -7.105 ;
$s5[5][8] =  1.842 ;
$s5[5][14] =  2.188 ;
$s5[5][15] =  9.148 ;
$s5[5][17] =  0.906 ;
$s5[5][23] =  7.006 ;
$s5[5][25] =  1.766 ;
$s5[5][26] =  8.122 ;
$s5[5][38] =  5.962 ;
$s5[5][44] =  6.908 ;
$s5[5][45] =  7.768 ;
$s5[7][8] =  8.108 ;
$s5[8][1] =  0.565 ;
$s5[8][2] =  0.786 ;
$s5[8][3] = -0.021 ;
$s5[8][4] =  0.995 ;
$s5[8][5] =  1.730 ;
$s5[8][7] =  8.286 ;
$s5[8][8] = -0.311 ;
$s5[8][9] =  8.413 ;
$s5[8][14] =  0.299 ;
$s5[8][15] =  0.177 ;
$s5[8][16] = -0.008 ;
$s5[8][17] =  0.068 ;
$s5[8][19] =  0.991 ;
$s5[8][23] =  0.217 ;
$s5[8][25] = -0.139 ;
$s5[8][26] =  0.055 ;
$s5[8][27] =  7.498 ;
$s5[8][28] =  7.343 ;
$s5[8][29] = -5.289 ;
$s5[8][31] =  5.246 ;
$s5[8][32] = -8.212 ;
$s5[8][34] = -7.561 ;
$s5[8][35] =-10.384 ;
$s5[8][36] = -1.839 ;
$s5[8][38] =  0.536 ;
$s5[8][39] =  4.531 ;
$s5[8][40] =  0.732 ;
$s5[8][41] = -0.310 ;
$s5[8][42] =  7.480 ;
$s5[8][44] =  0.069 ;
$s5[8][45] =  0.786 ;
$s5[8][46] =  0.699 ;
$s5[8][47] =  0.164 ;
$s5[9][1] = -1.235 ;
$s5[9][2] = -9.980 ;
$s5[9][3] = -9.808 ;
$s5[9][8] =  8.140 ;
$s5[9][9] = -1.634 ;
$s5[9][11] = -6.988 ;
$s5[9][14] = -0.026 ;
$s5[9][25] = -4.365 ;
$s5[9][26] = -8.607 ;
$s5[9][29] = -7.981 ;
$s5[9][45] =  5.246 ;
$s5[10][1] = -2.420 ;
$s5[10][2] =  6.908 ;
$s5[10][3] =  6.380 ;
$s5[10][10] = -4.433 ;
$s5[10][15] =  5.246 ;
$s5[10][26] =  5.246 ;
$s5[11][1] = -8.760 ;
$s5[11][2] = -7.484 ;
$s5[11][9] = -7.484 ;
$s5[14][1] = -0.104 ;
$s5[14][2] =  0.289 ;
$s5[14][3] =  0.210 ;
$s5[14][4] =  0.108 ;
$s5[14][5] =  1.894 ;
$s5[14][8] =  0.296 ;
$s5[14][9] = -0.225 ;
$s5[14][14] = -0.717 ;
$s5[14][15] =  0.441 ;
$s5[14][16] =  0.903 ;
$s5[14][17] =  8.122 ;
$s5[14][19] =  4.531 ;
$s5[14][23] =  2.104 ;
$s5[14][25] =  0.645 ;
$s5[14][26] =  1.125 ;
$s5[14][27] =  5.246 ;
$s5[14][28] =  6.192 ;
$s5[14][29] = -0.893 ;
$s5[14][31] = -0.043 ;
$s5[14][32] = -8.830 ;
$s5[14][36] =  4.531 ;
$s5[14][38] = -0.423 ;
$s5[14][40] = -6.139 ;
$s5[14][41] =  1.133 ;
$s5[14][42] =  1.187 ;
$s5[14][44] =  4.004 ;
$s5[14][45] =  6.347 ;
$s5[14][46] =  4.531 ;
$s5[15][1] =  0.947 ;
$s5[15][2] =  0.086 ;
$s5[15][3] =  0.089 ;
$s5[15][4] = -0.882 ;
$s5[15][5] =  9.524 ;
$s5[15][8] =  0.179 ;
$s5[15][10] =  5.246 ;
$s5[15][14] =  0.349 ;
$s5[15][15] = -0.150 ;
$s5[15][16] =  0.332 ;
$s5[15][17] = -0.588 ;
$s5[15][19] =  7.722 ;
$s5[15][20] =  6.192 ;
$s5[15][23] =  1.590 ;
$s5[15][25] =  1.031 ;
$s5[15][26] =  0.687 ;
$s5[15][27] =  1.715 ;
$s5[15][28] =  5.962 ;
$s5[15][29] =  4.531 ;
$s5[15][31] = -6.139 ;
$s5[15][32] =  4.531 ;
$s5[15][34] = -6.562 ;
$s5[15][36] =  7.255 ;
$s5[15][38] = -0.073 ;
$s5[15][40] = -0.299 ;
$s5[15][41] =  1.120 ;
$s5[15][42] =  1.539 ;
$s5[15][44] =  1.310 ;
$s5[15][45] =  1.109 ;
$s5[15][46] =  0.053 ;
$s5[15][47] = -6.988 ;
$s5[16][1] =  0.255 ;
$s5[16][2] =  1.091 ;
$s5[16][3] = -0.117 ;
$s5[16][4] = -7.837 ;
$s5[16][8] = -0.038 ;
$s5[16][14] = -0.722 ;
$s5[16][15] =  0.240 ;
$s5[16][16] =  0.734 ;
$s5[16][17] = -2.853 ;
$s5[16][19] =  4.531 ;
$s5[16][23] = -0.043 ;
$s5[16][25] =  0.492 ;
$s5[16][26] = -1.746 ;
$s5[16][27] =  7.971 ;
$s5[16][29] =  5.246 ;
$s5[16][34] = -6.635 ;
$s5[16][38] =  0.376 ;
$s5[16][41] = -2.104 ;
$s5[16][42] = -0.758 ;
$s5[16][44] =  0.376 ;
$s5[16][45] = -1.918 ;
$s5[17][1] =  0.952 ;
$s5[17][2] =  0.126 ;
$s5[17][3] =  0.714 ;
$s5[17][4] = -1.875 ;
$s5[17][5] =  1.183 ;
$s5[17][8] = -0.044 ;
$s5[17][14] =  8.208 ;
$s5[17][15] = -0.619 ;
$s5[17][16] = -2.506 ;
$s5[17][17] = -0.114 ;
$s5[17][19] = -0.673 ;
$s5[17][20] =  4.531 ;
$s5[17][23] =  5.962 ;
$s5[17][25] =  0.879 ;
$s5[17][26] =  0.370 ;
$s5[17][27] =  6.908 ;
$s5[17][28] =  7.096 ;
$s5[17][31] =  4.531 ;
$s5[17][36] =  5.665 ;
$s5[17][38] =  0.867 ;
$s5[17][40] =  6.380 ;
$s5[17][41] =  0.470 ;
$s5[17][42] =  5.246 ;
$s5[17][44] =  7.393 ;
$s5[17][45] =  1.329 ;
$s5[17][46] = -2.457 ;
$s5[18][1] = -8.896 ;
$s5[18][25] = -5.289 ;
$s5[18][26] = -6.139 ;
$s5[18][41] = -6.635 ;
$s5[19][1] = -0.365 ;
$s5[19][2] = -0.887 ;
$s5[19][3] =  0.552 ;
$s5[19][8] =  0.658 ;
$s5[19][14] =  4.531 ;
$s5[19][15] =  7.722 ;
$s5[19][16] =  4.531 ;
$s5[19][17] = -1.323 ;
$s5[19][19] = -0.181 ;
$s5[19][25] =  1.966 ;
$s5[19][26] =  1.717 ;
$s5[20][1] =  7.570 ;
$s5[20][2] =  7.179 ;
$s5[20][3] =  4.531 ;
$s5[20][15] =  6.192 ;
$s5[20][17] =  4.531 ;
$s5[20][20] =  6.380 ;
$s5[20][25] =  6.539 ;
$s5[20][26] =  4.531 ;
$s5[20][31] =  4.531 ;
$s5[22][1] = -7.981 ;
$s5[23][1] =  1.203 ;
$s5[23][2] = -0.294 ;
$s5[23][3] =  8.883 ;
$s5[23][5] =  7.570 ;
$s5[23][8] =  0.110 ;
$s5[23][14] =  1.375 ;
$s5[23][15] =  1.874 ;
$s5[23][16] = -1.608 ;
$s5[23][17] =  5.962 ;
$s5[23][23] =  0.138 ;
$s5[23][25] =  0.704 ;
$s5[23][26] =  6.780 ;
$s5[23][27] =  6.192 ;
$s5[23][36] = -6.635 ;
$s5[23][38] =  1.595 ;
$s5[23][41] =  7.179 ;
$s5[23][44] =  5.962 ;
$s5[23][45] =  5.246 ;
$s5[23][46] =  5.665 ;
$s5[25][1] =  0.346 ;
$s5[25][2] =  0.064 ;
$s5[25][3] =  0.015 ;
$s5[25][4] =  0.780 ;
$s5[25][5] =  0.748 ;
$s5[25][8] = -0.230 ;
$s5[25][9] = -5.182 ;
$s5[25][14] =  0.530 ;
$s5[25][15] =  1.077 ;
$s5[25][16] =  0.503 ;
$s5[25][17] =  0.930 ;
$s5[25][18] = -5.289 ;
$s5[25][19] =  1.484 ;
$s5[25][20] =  6.539 ;
$s5[25][23] =  0.704 ;
$s5[25][25] =  0.236 ;
$s5[25][26] =  0.505 ;
$s5[25][27] = -0.177 ;
$s5[25][29] =  0.802 ;
$s5[25][31] =  4.531 ;
$s5[25][34] = -9.058 ;
$s5[25][36] = -0.889 ;
$s5[25][38] =  2.300 ;
$s5[25][40] =  7.006 ;
$s5[25][41] =  1.125 ;
$s5[25][42] =  0.267 ;
$s5[25][44] =  6.614 ;
$s5[25][45] =  0.910 ;
$s5[25][46] = -0.319 ;
$s5[25][47] = -3.142 ;
$s5[26][1] =  0.469 ;
$s5[26][2] =  0.274 ;
$s5[26][3] =  0.174 ;
$s5[26][4] = -0.486 ;
$s5[26][5] = -0.070 ;
$s5[26][7] =  5.246 ;
$s5[26][8] = -0.196 ;
$s5[26][9] = -9.777 ;
$s5[26][10] =  5.246 ;
$s5[26][11] = -6.139 ;
$s5[26][14] = -0.107 ;
$s5[26][15] =  0.935 ;
$s5[26][16] = -1.433 ;
$s5[26][17] =  0.258 ;
$s5[26][18] = -6.139 ;
$s5[26][19] =  0.940 ;
$s5[26][20] =  4.531 ;
$s5[26][23] =  1.619 ;
$s5[26][25] =  0.552 ;
$s5[26][26] =  0.401 ;
$s5[26][27] =  2.685 ;
$s5[26][28] =  7.894 ;
$s5[26][29] =  0.957 ;
$s5[26][34] = -8.986 ;
$s5[26][36] = -0.323 ;
$s5[26][38] =  1.659 ;
$s5[26][40] =  0.957 ;
$s5[26][41] =  1.079 ;
$s5[26][42] =  0.498 ;
$s5[26][44] =  0.905 ;
$s5[26][45] =  0.330 ;
$s5[26][46] = -0.162 ;
$s5[26][47] = -1.299 ;
$s5[27][1] = -1.587 ;
$s5[27][2] =  0.732 ;
$s5[27][3] =  2.787 ;
$s5[27][8] =  7.700 ;
$s5[27][14] =  5.962 ;
$s5[27][15] =  1.336 ;
$s5[27][16] =  8.007 ;
$s5[27][17] =  6.908 ;
$s5[27][23] =  6.192 ;
$s5[27][25] = -0.177 ;
$s5[27][26] =  0.673 ;
$s5[27][27] =  0.796 ;
$s5[27][38] =  5.246 ;
$s5[28][1] =  8.311 ;
$s5[28][2] =  8.569 ;
$s5[28][3] = -9.183 ;
$s5[28][8] =  7.690 ;
$s5[28][14] =  6.192 ;
$s5[28][15] =  8.610 ;
$s5[28][17] =  7.096 ;
$s5[28][26] =  7.811 ;
$s5[28][38] =  8.437 ;
$s5[29][1] =  0.143 ;
$s5[29][2] = -0.321 ;
$s5[29][3] =  1.832 ;
$s5[29][4] =  5.962 ;
$s5[29][8] = -5.289 ;
$s5[29][9] = -8.686 ;
$s5[29][14] = -0.673 ;
$s5[29][15] =  5.246 ;
$s5[29][16] =  5.246 ;
$s5[29][25] =  0.802 ;
$s5[29][26] =  0.241 ;
$s5[29][29] =  5.962 ;
$s5[29][36] =  5.665 ;
$s5[31][1] = -0.710 ;
$s5[31][2] =  5.246 ;
$s5[31][8] =  5.246 ;
$s5[31][14] = -0.043 ;
$s5[31][15] = -6.139 ;
$s5[31][17] =  4.531 ;
$s5[31][20] =  4.531 ;
$s5[31][25] =  4.531 ;
$s5[31][36] =  5.665 ;
$s5[32][1] = -2.031 ;
$s5[32][2] = -7.484 ;
$s5[32][3] =  5.246 ;
$s5[32][8] = -8.212 ;
$s5[32][14] = -8.830 ;
$s5[32][15] =  4.531 ;
$s5[32][32] = -6.139 ;
$s5[32][46] = -6.139 ;
$s5[34][1] = -2.650 ;
$s5[34][2] = -7.096 ;
$s5[34][8] = -7.229 ;
$s5[34][15] = -5.216 ;
$s5[34][25] = -8.124 ;
$s5[34][44] = -6.139 ;
$s5[34][46] = -5.289 ;
$s5[35][1] = -7.484 ;
$s5[35][8] =-10.384 ;
$s5[35][35] = -6.139 ;
$s5[36][1] =  0.087 ;
$s5[36][2] =  5.962 ;
$s5[36][3] = -0.311 ;
$s5[36][4] = -5.289 ;
$s5[36][8] = -1.641 ;
$s5[36][14] =  4.531 ;
$s5[36][15] =  7.006 ;
$s5[36][17] =  4.531 ;
$s5[36][23] = -5.289 ;
$s5[36][25] = -0.772 ;
$s5[36][26] =  0.241 ;
$s5[36][41] =  4.531 ;
$s5[36][45] =  4.531 ;
$s5[36][46] = -6.635 ;
$s5[38][1] =  0.640 ;
$s5[38][2] =  0.522 ;
$s5[38][3] =  1.592 ;
$s5[38][4] =  6.799 ;
$s5[38][5] =  6.908 ;
$s5[38][8] =  0.450 ;
$s5[38][14] = -0.892 ;
$s5[38][15] =  0.163 ;
$s5[38][16] = -0.043 ;
$s5[38][17] =  1.752 ;
$s5[38][23] =  1.595 ;
$s5[38][25] =  2.300 ;
$s5[38][26] =  1.989 ;
$s5[38][27] =  5.246 ;
$s5[38][28] =  8.437 ;
$s5[38][36] = -9.183 ;
$s5[38][38] =  0.082 ;
$s5[38][40] =  6.192 ;
$s5[38][41] =  1.715 ;
$s5[38][42] =  6.192 ;
$s5[38][44] = -2.104 ;
$s5[38][45] =  2.941 ;
$s5[39][3] =  5.246 ;
$s5[39][8] =  4.531 ;
$s5[40][1] = -1.351 ;
$s5[40][2] = -0.158 ;
$s5[40][3] =  0.267 ;
$s5[40][8] =  0.732 ;
$s5[40][14] = -6.139 ;
$s5[40][15] = -0.299 ;
$s5[40][17] =  5.665 ;
$s5[40][25] =  6.799 ;
$s5[40][26] =  0.538 ;
$s5[40][38] =  5.962 ;
$s5[40][40] = -1.742 ;
$s5[40][45] = -6.139 ;
$s5[41][1] =  0.637 ;
$s5[41][2] =  1.077 ;
$s5[41][3] =  1.037 ;
$s5[41][4] = -5.289 ;
$s5[41][8] = -0.282 ;
$s5[41][14] =  1.375 ;
$s5[41][15] =  1.217 ;
$s5[41][16] = -2.104 ;
$s5[41][17] = -0.087 ;
$s5[41][18] = -6.635 ;
$s5[41][23] =  6.192 ;
$s5[41][25] =  1.303 ;
$s5[41][26] =  2.409 ;
$s5[41][36] =  4.531 ;
$s5[41][38] =  1.431 ;
$s5[41][41] = -1.567 ;
$s5[41][42] =  6.677 ;
$s5[41][44] =  1.388 ;
$s5[41][45] =  1.729 ;
$s5[41][46] = -0.893 ;
$s5[42][1] =  0.037 ;
$s5[42][2] =  0.764 ;
$s5[42][3] =  0.031 ;
$s5[42][4] = -5.289 ;
$s5[42][8] =  7.633 ;
$s5[42][14] =  0.582 ;
$s5[42][15] =  1.256 ;
$s5[42][16] =  4.531 ;
$s5[42][17] =  5.246 ;
$s5[42][25] =  0.267 ;
$s5[42][26] =  0.935 ;
$s5[42][38] =  6.192 ;
$s5[42][41] =  7.096 ;
$s5[42][42] =  0.495 ;
$s5[42][45] =  5.665 ;
$s5[42][46] =  5.246 ;
$s5[42][47] =  4.531 ;
$s5[44][1] =  1.416 ;
$s5[44][2] =  2.757 ;
$s5[44][3] = -0.079 ;
$s5[44][4] = -0.970 ;
$s5[44][5] =  7.455 ;
$s5[44][8] = -0.010 ;
$s5[44][14] =  5.138 ;
$s5[44][15] =  0.196 ;
$s5[44][16] =  0.376 ;
$s5[44][17] =  7.971 ;
$s5[44][19] =  5.246 ;
$s5[44][23] =  6.192 ;
$s5[44][25] =  0.887 ;
$s5[44][26] =  0.905 ;
$s5[44][28] = -6.635 ;
$s5[44][34] = -7.837 ;
$s5[44][38] =  0.042 ;
$s5[44][41] =  1.717 ;
$s5[44][44] = -2.410 ;
$s5[44][45] = -2.238 ;
$s5[45][1] = -0.115 ;
$s5[45][2] = -0.265 ;
$s5[45][3] = -0.492 ;
$s5[45][4] = -2.019 ;
$s5[45][5] =  8.258 ;
$s5[45][8] =  0.557 ;
$s5[45][9] =  5.246 ;
$s5[45][14] = -1.248 ;
$s5[45][15] =  0.736 ;
$s5[45][16] = -2.291 ;
$s5[45][17] =  0.751 ;
$s5[45][19] =  5.246 ;
$s5[45][23] =  5.246 ;
$s5[45][25] =  0.402 ;
$s5[45][26] =  0.330 ;
$s5[45][36] =  6.908 ;
$s5[45][38] =  2.200 ;
$s5[45][40] = -6.139 ;
$s5[45][41] =  2.997 ;
$s5[45][42] =  5.962 ;
$s5[45][44] = -2.238 ;
$s5[45][45] = -1.372 ;
$s5[46][1] =  0.728 ;
$s5[46][2] =  7.971 ;
$s5[46][3] = -1.005 ;
$s5[46][4] = -5.289 ;
$s5[46][8] =  0.546 ;
$s5[46][14] =  4.531 ;
$s5[46][15] = -0.358 ;
$s5[46][17] = -2.172 ;
$s5[46][23] =  5.665 ;
$s5[46][25] = -0.721 ;
$s5[46][26] = -0.162 ;
$s5[46][32] = -6.139 ;
$s5[46][34] = -5.289 ;
$s5[46][36] = -6.635 ;
$s5[46][41] = -0.893 ;
$s5[46][42] =  5.246 ;
$s5[47][1] =  0.715 ;
$s5[47][2] = -1.857 ;
$s5[47][3] = -6.139 ;
$s5[47][8] = -0.254 ;
$s5[47][15] = -8.333 ;
$s5[47][25] = -1.847 ;
$s5[47][26] = -1.299 ;
$s5[47][38] =  4.531 ;
$s5[47][42] =  4.531 ;


$ntat  = 0 ; # sum of all atoms (atom types)
$ntpq  = 0 ; # sum of all atom 1-2 pairs (bonds)
$ntpq2 = 0 ; # sum of all atom 1-3 pairs (angles)
$ntpq3 = 0 ; # sum of all atom 1-4 pairs (torsions)
$ntpq4 = 0 ; # sum of all atom 1-5 pairs
$ntpq5 = 0 ; # sum of all atom 1-6 pairs
$sco0  = 0 ; # score of atom types
$sco   = 0 ; # score of 1-2 pairs
$sco2  = 0 ; # score of 1-3 pairs
$sco3  = 0 ; # score of 1-4 pairs
$sco4  = 0 ; # score of 1-5 pairs 
$sco5  = 0 ; # score of 1-6 pairs


for ($ll=1; $ll<=$natoms; $ll++) {

# assign numeric atomtype and neighbour list

    $el = $elem[$ll] ;  # element
    $t0 = $typ[$ll] ;  # given atom type

    if ($el eq "H") {
      $t0 = "H?" ;
    }
    else {
      $ntat++ ; # count all non-hydrogen atoms
    }

    if ($t0 eq "**") {
      if ($el eq "C")  { $t0 = "C?" }
      if ($el eq "N")  { $t0 = "N?" }
      if ($el eq "O")  { $t0 = "O?" }
      if ($el eq "Si") { $t0 = "SI" }
      if ($el eq "P")  { $t0 = "P?" }
      if ($el eq "S")  { $t0 = "S?" }
      if ($el eq "S")  { $t0 = "S?" }
      if ($el eq "F")  { $t0 = "F" }
      if ($el eq "Cl") { $t0 = "CL" }
      if ($el eq "Br") { $t0 = "BR" }
      if ($el eq "I")  { $t0 = "I" }
    }

    $tt[$ll] = 0 ; # default is hydrogen
    for ($j=1 ; $j <= $maxt; $j++) {
      if ($t0 eq $t[$j]) {
        $sco0 = $sco0 + $p0[$j] ; # count frequency
        $tt[$ll] = $j ; # assign numeric atom type
      }
    }


    $nn[$ll] = $nsubst[$ll] ; # number of neighbours of atom 

    if ($nn[$ll] == 0) {
      $nb[$ll][1] = 0 ;
    }
    else {
      for ($j=1 ; $j <= $nn[$ll]; $j++) {
        $nb[$ll][$j] = $neib[$ll][$j] ; # neighbour matrix of atom 
      }
    }

} 


# determine frequencies 

for ($j=1 ; $j <= $natoms; $j++) {
  if ($tt[$j] != 0) {   
    if ($nn[$j] > 0) {
      for ($kk=1 ; $kk <= $nn[$j] ; $kk++) {
        $q1 = $nb[$j][$kk] ; # first neighbour of atom i
        if ($tt[$q1] != 0 and $q1 != $j) {
          $sco = $sco + $s[$tt[$j]][$tt[$q1]] ;
          $ntpq++ ;
        } 
        if ($nn[$q1] != 0 and $q1 != $j) {
          for ($ll=1 ; $ll <= $nn[$q1] ; $ll++) {
            $q2 = $nb[$q1][$ll] ; # second neighbour of atom i
            if ($tt[$q2] != 0 and $q2 != $j and $q2 != $q1) {
              $sco2 = $sco2 + $s2[$tt[$j]][$tt[$q2]] ;
              $ntpq2++ ;
            }
            if ($nn[$q2] != 0 and $q2 != $j) {
              for ($mm=1 ; $mm <= $nn[$q2] ; $mm++) {
                $q3 = $nb[$q2][$mm] ; # third neighbour of atom i
                if ($tt[$q3] != 0 and $q3 != $q1) {
                  $sco3 = $sco3 + $s3[$tt[$j]][$tt[$q3]] ;
                  $ntpq3++ ;
                }
                if ($nn[$q3]!=0 and $q3!=$j and $q3!=$q1 and $q3!=$q2) {
                  for ($oo=1 ; $oo <= $nn[$q3] ; $oo++) {
                    $q4 = $nb[$q3][$oo] ; # fourth neighbour
                    if ($tt[$q4]!=0 and $q4!=$q3) { 
                      $sco4 = $sco4 + $s4[$tt[$j]][$tt[$q4]] ;
                      $ntpq4++ ;
                    }
                    if ($nn[$q4]!=0 and $q4!=$j and $q4!=$q1 and $q4!=$q2
                        and $q4!=$q3) {
                      for ($pp=1 ; $pp <= $nn[$q4] ; $pp++) {
                        $q5 = $nb[$q4][$pp] ; # fifth neighbour
                        if ($tt[$q5]!=0 and $q5!=$q4) {
                          $sco5 = $sco5 + $s5[$tt[$j]][$tt[$q5]] ;
                          $ntpq5++
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
  }
}

if ($ntat  > 0) {$sco0 = $sco0 / $ntat}
if ($ntpq  > 0) {$sco  = $sco / $ntpq}
if ($ntpq2 > 0) {$sco2 = $sco2 / $ntpq2}
if ($ntpq3 > 0) {$sco3 = $sco3 / $ntpq3}
if ($ntpq4 > 0) {$sco4 = $sco4 / $ntpq4}
if ($ntpq5 > 0) {$sco5 = $sco5 / $ntpq5}

$scsum = $sco0 + $sco + $sco2 + $sco3 + $sco4 + $sco5 ;

print "$hinname $scsum\n" ;

exit(0)

