#!/bin/bash

##    Alberto Roldan --> 11/2013

# input -->  ibrion-1/charge/charge_SUM/
# input -->  OUTCAR
# input -->  POSCAR


#---------------------------------------------------------------------------
    if [ $@ ]; then
        surface=$@
    else
       echo " --- Introduce the surface --- "
    fi
#---------------------------------------------------------------------------
    if [ $surface == 001 ]; then
 A=(6 7 8)                    # Tetrahedral
 B=(21 22 23 24)              # Octahedral
 S=(49 50 51 52 53 54 55 56)  # Sulphurs
    elif [ $surface == 011 ]; then
 A=(5 6 8)                    # Tetrahedral
 B=(21 22 23 24)              # Octahedral
 S=(41 42 47 48 49 50 51 52 53 54 55 56)  # Sulphurs
    elif [ $surface == 111 ]; then ###################################################
 A=(7 8)                    # Tetrahedral
 B=(24)              # Octahedral
 S=(49 50 51 52 53 54 55 56)  # Sulphurs
    fi
echo " --- Surface $surface ---"

#---------------------------------------------------------------------------
   folder="./ibrion-1/charge/charge_SUM/"
   DOSfolder="./ibrion-1/DOS/LDOS/"
 # Q = electrons 
 # Ms = spin density
 # Z = cartesian value of Z
 # G = Grimme energy eV
 # dup = d-band centre UP
 # ddo = d-band centre DOWM

     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- A
  for i in "${A[@]}"
    do
     QA=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsA=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZA=$(grep " $i " $folder/ACF.dat | awk '{print $4}')
     dupA=$(grep "D-up-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ddoA=$(grep "D-down-centre =" $DOSfolder/$i.xy | awk '{print $4}')

     q=$(echo "$q $QA "); ms=$(echo "$ms $MsA "); z=$(echo "$z $ZA ");
     dup=$(echo "$dup $dupA "); ddo=$(echo "$ddo $ddoA ");
    done
   qa=$q; msa=$ms; za=$z; dupa=$dup; ddoa=$ddo;
     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- B
  for i in "${B[@]}"
    do
     QB=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsB=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZB=$(grep " $i " $folder/ACF.dat | awk '{print $4}')
     dupB=$(grep "D-up-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ddoB=$(grep "D-down-centre =" $DOSfolder/$i.xy | awk '{print $4}')

     q=$(echo "$q $QB "); ms=$(echo "$ms $MsB "); z=$(echo "$z $ZB ");
     dup=$(echo "$dup $dupB "); ddo=$(echo "$ddo $ddoB ");
    done
   qb=$q; msb=$ms; zb=$z; dupb=$dup; ddob=$ddo;
     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- S
  for i in "${S[@]}"
    do
     QS=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsS=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZS=$(grep " $i " $folder/ACF.dat | awk '{print $4}')

     q=$(echo "$q $QS "); ms=$(echo "$ms $MsS "); z=$(echo "$z $ZS ");
    done
   qs=$q; mss=$ms; zs=$z; 
     q=''; ms=''; z=''; 

#--------------------------------------------------------------------------- Grimme
   G=$(grep -A1 "Grimme" OUTCAR |tail -1 | awk '{print $5}')

#---------------------------------------------------------------------------
 nO=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $3}')
 nH=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $4}')

     q=''; z=''
#--------------------------------------------------------------------------- nO
   plusnO=$(echo "56 + $nO" | bc)
   for (( i=57; i<=$plusnO; i++ ))
    do
     QO=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     ZO=$(grep " $i " $folder/ACF.dat | awk '{print $4}')

     q=$(echo "$q $QO ");  z=$(echo "$z $ZO ");
    done
   qo=$q; zo=$z
     q=''; z=''
#--------------------------------------------------------------------------- nH
   plusnH=$(echo "$plusnO + $nH" | bc)
   for (( i=$plusnO+1; i<=$plusnH; i++ ))
    do
     QH=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     ZH=$(grep " $i " $folder/ACF.dat | awk '{print $4}')

     q=$(echo "$q $QH ");  z=$(echo "$z $ZH ");
    done
   qh=$q; zh=$z
     q=''; z=''

#---------------------------------------------------------------------------
 echo "Grimme"  >> Data.dat ;               echo " Q(A) Q(B) Ms(A) Ms(B)"  >> Data.dat ;
 echo " d(A)-down d(B)-down"  >> Data.dat ; echo " d(A)-up d(B)-up"  >> Data.dat ;
 echo " z(A) z(B)"  >> Data.dat ;           echo " Z(S)"  >> Data.dat;      
 echo " Z(O) Q(O)"  >> Data.dat;            echo " Z(H) Q(H)"  >> Data.dat;
 echo "  $G"  >> Data.dat
 echo " $qa $qb $msa $msb"  >> Data.dat;
 echo " $ddoa $ddob"  >> Data.dat;          echo " $dupa $dupb"  >> Data.dat;
 echo " $za $zb"  >> Data.dat;              echo " $zs"  >> Data.dat
 echo " $zo $qo"  >> Data.dat;              echo " $zh $qh"  >> Data.dat


