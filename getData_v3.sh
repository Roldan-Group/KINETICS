#!/bin/bash

##    Alberto Roldan --> 11/2013

# input -->  ibrion-1/charge/charge_SUM/
# input -->  OUTCAR
# input -->  POSCAR


#---------------------------------------------------------------------------
    if [ $@ ]; then
        TM=$@
    else
       TM=Fe
    fi
#---------------------------------------------------------------------------
    if [ $TM == Fe ]; then
 M=(8)
 A=(7)                        # Tetrahedral
 B=(24)                       # Octahedral
 S=(49 50 51 52 53 54 55 56)  # Sulphurs
    else
 M=(56)
 A=(7)                        # Tetrahedral
 B=(23)                       # Octahedral
 S=(48 49 50 51 52 53 54 55)  # Sulphurs
    fi
echo " --- Transition Metal is $TM (atom n=$M) ---"

#---------------------------------------------------------------------------
   folder="./ibrion-1/charge/charge_SUM/"
   DOSfolder="./ibrion-1/DOS/LDOS/"
   Workfunction="./ibrion-1/workfunction/Workfunction"
 vasp2xyz.pl CONTCAR

 # Q = electrons 
 # Ms = spin density
 # Z = cartesian value of Z
 # G = Grimme energy eV
 # dup = d-band centre UP
 # ddo = d-band centre DOWM

     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- TM
  for i in "${M[@]}"
    do
      line=$(($i + 2))
     QM=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsM=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZM=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')
     xM=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $2}')
     yM=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $3}')

     dupM=$(grep "D-up-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ddoM=$(grep "D-down-centre =" $DOSfolder/$i.xy | awk '{print $4}')

     q=$(echo "$q $QM "); ms=$(echo "$ms $MsM "); z=$(echo "$z $ZM ");
     dup=$(echo "$dup $dupM "); ddo=$(echo "$ddo $ddoM ");
    done
   qm=$q; msm=$ms; zm=$z; dupm=$dup; ddom=$ddo;
     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- A
  for i in "${A[@]}"
    do
      line=$(($i + 2))
     QA=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsA=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     dupA=$(grep "D-up-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ddoA=$(grep "D-down-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ZA=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')

     q=$(echo "$q $QA "); ms=$(echo "$ms $MsA "); z=$(echo "$z $ZA ");
     dup=$(echo "$dup $dupA "); ddo=$(echo "$ddo $ddoA ");
    done
   qa=$q; msa=$ms; za=$z; dupa=$dup; ddoa=$ddo;
     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- B
  for i in "${B[@]}"
    do
      line=$(($i + 2))
     QB=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsB=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZB=$(grep " $i " $folder/ACF.dat | awk '{print $4}')
     dupB=$(grep "D-up-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ddoB=$(grep "D-down-centre =" $DOSfolder/$i.xy | awk '{print $4}')
     ZB=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')

     q=$(echo "$q $QB "); ms=$(echo "$ms $MsB "); z=$(echo "$z $ZB ");
     dup=$(echo "$dup $dupB "); ddo=$(echo "$ddo $ddoB ");
    done
   qb=$q; msb=$ms; zb=$z; dupb=$dup; ddob=$ddo;
     q=''; ms=''; z=''; ddo=''; dup='';
#--------------------------------------------------------------------------- S
  for i in "${S[@]}"
    do
      line=$(($i + 2))
     QS=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     MsS=$(grep " $i " $folder/SpinDensity/ACF.dat | awk '{print $5}')
     ZS=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')

     q=$(echo "$q $QS "); ms=$(echo "$ms $MsS "); z=$(echo "$z $ZS ");
    done
   qs=$q; mss=$ms; zs=$z; 
     q=''; ms=''; z=''; 

#---------------------------------------------------------------------------
  if [ $TM == Fe ]; then
     n1=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $3}')
     n2=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $4}')
     n3=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $5}')
    else
     n1=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $4}')
     n2=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $5}')
     n3=$(grep -B1 "Selective dynamics" POSCAR | head -1 | awk '{print $6}')
    fi
     q=''; z=''; d1=''; j=0;
#--------------------------------------------------------------------------- n1
  if [ $n1 ]; then
   plusn1=$(echo "56 + $n1" | bc)
   for (( i=57; i<=$plusn1; i++ ))
    do
      line=$(($i + 2))
     j=$(echo "1 + $j" | bc)
     Q1=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     Z1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')
     x1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $2}')
     y1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $3}')
     D1M=$(echo "scale = 5; sqrt(($x1 - $xM)*($x1 - $xM) + ($y1 - $yM)*($y1 - $yM) + ($Z1 - $ZM)*($Z1 - $ZM) )" | bc)

      x1=$x1; y1=$y1; z1=$Z1;

     q=$(echo "$q $Q1 ");  z=$(echo "$z $Z1 "); d1=$(echo "$d1 $D1M "); 
    done
   q1=$q; z1=$z; d1M=$d1;
  fi
     q=''; z=''; d2=''; 
#--------------------------------------------------------------------------- n2
  if [ $n2 ]; then
   plusn2=$(echo "$plusn1 + $n2" | bc)
   for (( i=$plusn1+1; i<=$plusn2; i++ ))
    do
      line=$(($i + 2))
     Q2=$(grep " $i " $folder/ACF.dat | awk '{print $5}')
     Z2=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')
     x2=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $2}')
     y2=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $3}')
      for (( j=57; j<=$plusn1; j++ ))
       do
      line=$(($j + 2))
        Z1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $4}')
        x1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $2}')
        y1=$(sed -n "${line}p" CONTCAR.xyz | awk '{print $3}')
         D21=$(echo "scale = 5; sqrt(($x1 - $x2)*($x1 - $x2) + ($y1 - $y2)*($y1 - $y2) + ($z1 - $Z2)*($z1 - $Z2))" | bc)
         d2=$(echo "$d2 $D21 ");
       done

     q=$(echo "$q $Q2 ");  z=$(echo "$z $Z2 ");
    done
   q2=$q; z2=$z; d21=$d2;
  fi
     q=''; z=''


 rm CONTCAR.xyz

#---------------------------------------------------------------------------
 echo "" >> Data.dat ;
 pwd >> Data.dat ;
 echo " Energy (eV) Grimme (eV) WorkFunction (eV)"  >> Data.dat;
 echo " Q(TM) Ms(TM) d(TM)-down d(TM)-up Z(TM)"  >> Data.dat;
 echo " Q(A) Q(B) Ms(A) Ms(B) d(A)-down d(B)-down d(A)-up d(B)-up"  >> Data.dat;
 echo " z(A) z(B) Z(S)"  >> Data.dat;
 if [ $n1 ]; then
   echo " Q(1) Z(1) d[1-TM]"   >> Data.dat;
 fi
 if [ $n2 ]; then
   echo " Q(2) Z(2) d[2-1]"   >> Data.dat;
 fi

energy=$(pp OUTCAR | awk '{print $2}')
G=$(grep "Estimated vdW energy (eV):" OUTCAR | tail -1 | awk '{print $5}')
workfunction=$(grep -A1 WorkFunction $Workfunction |tail -1 | awk '{print $3}')
 echo "$energy $G $workfunction"  >> Data.dat
     q=''; ms=''; z=''; ddo=''; dup='';
     q=''; ms=''; z=''; ddo=''; dup='';
 echo "$qm $msm $ddom $dupm $zm"  >> Data.dat;
 echo "$qa $qb $msa $msb $ddoa $ddob $dupa $dupb"  >> Data.dat;
 echo " $za $zb $zs"  >> Data.dat;
 if [ $n1 ]; then
   echo "$q1 $z1 $d1M"  >> Data.dat;
 fi
 if [ $n2 ]; then
   echo "$q2 $z2 $d21"  >> Data.dat;
 fi
