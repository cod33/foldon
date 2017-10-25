#!/usr/bin/sh

#removes junk at beginning
sed '1,121d' 1rfo_oriented.pdb > temp.pdb
#removes chain ID
sed -i -r 's/\S+//5' temp.pdb
#removes junk at end of each row
sed -i -r 's/\S+//9' temp.pdb
sed -i -r 's/\S+//9' temp.pdb
sed -i -r 's/\S+//9' temp.pdb
sed -i -r 's/\S+//9' temp.pdb
#removes junk in last row
sed -i '$ d' temp.pdb
sed -i '$ d' temp.pdb
mv temp.pdb 1rfo_21_renum.pdb

#renumbers chain B but with shitty formatting
#field: Atom, atomno, space, atom name, alt conf indic, resname, space, chain indent, res seq no, spaces, x, y, z, occup, temp fact, spaces, seg id
# $1, $2, " ", $3, " ", $4, " ", " ", " ", $5, "    ", $6, $7, $8

#awk '{printf"%6s %5d %1s %-4s %1s %3s %1s %1s %4s %4s %8s %8s %8s %6s\n",$1, $2, " ", $3, " ", $4, " ", " ", " ", $5, "    ", $6, $7, $8;}' temp.pdb > temp2.pdb
#awk 'NR>=431{$5 += 27; printf"%6s %5d %1s %-4s %1s %3s %1s %1s %4s %4s %8s %8s %8s %6s\n",
#                              $1, $2, " ", $3, " ", $4, " ", " ", " ", $5, "    ", $6, $7, $8; }' temp2.pdb > temp3.pdb
##you need to move temp2 into temp before you do the thing
##ok so you're gonna have to replace the last lines of temp with temp2, and be careful about formatting. 
##cp temp2.pdb temp2.pdb~
##awk 'NR == FNR{if(FNR >= 434) {patch = patch $0 ORS} next } FNR == 0{ $0 = patch $0} 1' temp.pdb temp2.pdb~
##sed -i "$(sed 1,434!d temp.pdb),\$d" temp.pdb > temp2.pdb~
#
#sed '431,$d' temp2.pdb > temp2.pdb~
#cat temp2.pdb~ > temp4.pdb
#tail -n +2 temp3.pdb >> temp4.pdb
#
##[%4f ATOM][%6f atomnum] [%5s atomtype] [%3f residue] [%8s resid] [9 x] [8 y] [8 z]
##renumbers chain C but with shitty formatting
#awk 'NR>=864{$5 += 27; printf"%6s %5d %1s %-4s %1s %3s %1s %1s %4s %4s %8s %8s %8s %6s\n",
#                              $1, $2, " ", $3, " ", $4, " ", " ", " ", $5, "    ", $6, $7, $8;}' temp4.pdb > temp5.pdb
#sed '864,$d' temp4.pdb >temp4.pdb~
#cat temp4.pdb~ > 1rfo_30_renum.pdb
#tail -n +2 temp5.pdb >> 1rfo_30_renum.pdb
#
#rm temp*
