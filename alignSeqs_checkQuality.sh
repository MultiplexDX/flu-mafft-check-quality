## © 2020, Viera Kovacova for MultiplexDX. All rights reserved. 
# sequences downloaded from gisaid.org 
# filtering in gisaid: 1.1.2018-24.6.2020, human, 4 flu types (IAV H1N1, IAV H3N2, IBV Victori, IBV Yamagata); 3 segments (PB1, PB2, PA)
# bash code;  environment: unix Os, ubuntu v 20.04 LTS; hardware: Dell XPS 15
# for every single downloaded flu record, we have sequences of 3 segments

# separate sequences of three segments (PB1, PB2 and PA) into new fasta files
cd ~/Documents/flu
awk ' BEGIN{RS=">"}; /\|PB1\|2/ { print ">"$0 } ' IAV_H1N1/gisaid_epiflu_sequence.fasta > IAV_H1N1/H1N1_PB1_segm2.fa
awk ' BEGIN{RS=">"}; /\|PB2\|1/ { print ">"$0 } ' IAV_H1N1/gisaid_epiflu_sequence.fasta > IAV_H1N1/H1N1_PB2_segm1.fa
awk ' BEGIN{RS=">"}; /\|PA\|3/ { print ">"$0 } ' IAV_H1N1/gisaid_epiflu_sequence.fasta > IAV_H1N1/H1N1_PA_segm3.fa

awk ' BEGIN{RS=">"}; /\|PB1\|2/ { print ">"$0 } ' IAV_H3N2/gisaid_epiflu_sequence.fasta > IAV_H3N2/H3N2_PB1_segm2.fa
awk ' BEGIN{RS=">"}; /\|PB2\|1/ { print ">"$0 } ' IAV_H3N2/gisaid_epiflu_sequence.fasta > IAV_H3N2/H3N2_PB2_segm1.fa
awk ' BEGIN{RS=">"}; /\|PA\|3/ { print ">"$0 } ' IAV_H3N2/gisaid_epiflu_sequence.fasta > IAV_H3N2/H3N2_PA_segm3.fa

awk ' BEGIN{RS=">"}; /\|PB1\|2/ { print ">"$0 } ' IBV_Victoria/gisaid_epiflu_sequence.fasta > IBV_Victoria/Vict_PB1_segm2.fa
awk ' BEGIN{RS=">"}; /\|PB2\|1/ { print ">"$0 } ' IBV_Victoria/gisaid_epiflu_sequence.fasta > IBV_Victoria/Vict_PB2_segm1.fa
awk ' BEGIN{RS=">"}; /\|PA\|3/ { print ">"$0 } ' IBV_Victoria/gisaid_epiflu_sequence.fasta > IBV_Victoria/Vict_PA_segm3.fa

awk ' BEGIN{RS=">"}; /\|PB1\|2/ { print ">"$0 } ' IBV_Yamagata/gisaid_epiflu_sequence.fasta > IBV_Yamagata/Yama_PB1_segm2.fa
awk ' BEGIN{RS=">"}; /\|PB2\|1/ { print ">"$0 } ' IBV_Yamagata/gisaid_epiflu_sequence.fasta > IBV_Yamagata/Yama_PB2_segm1.fa
awk ' BEGIN{RS=">"}; /\|PA\|3/ { print ">"$0 } ' IBV_Yamagata/gisaid_epiflu_sequence.fasta > IBV_Yamagata/Yama_PA_segm3.fa

# align sequences of the given segments, per flu strain, using mafft software (version)
mafft --thread 4 --auto IAV_H1N1/H1N1_PB1_segm2.fa > IAV_H1N1/H1N1_PB1_segm2_mafft.fa
mafft --thread 4 --auto IAV_H1N1/H1N1_PB2_segm1.fa > IAV_H1N1/H1N1_PB2_segm1_mafft.fa
mafft --thread 4 --auto IAV_H1N1/H1N1_PA_segm3.fa > IAV_H1N1/H1N1_PA_segm3_mafft.fa

mafft --thread 4 --auto IAV_H3N2/H3N2_PB1_segm2.fa > IAV_H3N2/H3N2_PB1_segm2_mafft.fa
mafft --thread 4 --auto IAV_H3N2/H3N2_PB2_segm1.fa > IAV_H3N2/H3N2_PB2_segm1_mafft.fa
mafft --thread 4 --auto IAV_H3N2/H3N2_PA_segm3.fa > IAV_H3N2/H3N2_PA_segm3_mafft.fa

mafft --thread 4 --auto IBV_Victoria/Vict_PB1_segm2.fa > IBV_Victoria/Vict_PB1_segm2_mafft.fa
mafft --thread 4 --auto IBV_Victoria/Vict_PB2_segm1.fa > IBV_Victoria/Vict_PB2_segm1_mafft.fa
mafft --thread 4 --auto IBV_Victoria/Vict_PA_segm3.fa > IBV_Victoria/Vict_PA_segm3_mafft.fa

mafft --thread 4 --auto IBV_Yamagata/Yama_PB1_segm2.fa > IBV_Yamagata/Yama_PB1_segm2_mafft.fa
mafft --thread 4 --auto IBV_Yamagata/Yama_PB2_segm1.fa > IBV_Yamagata/Yama_PB2_segm1_mafft.fa
mafft --thread 4 --auto IBV_Yamagata/Yama_PA_segm3.fa > IBV_Yamagata/Yama_PA_segm3_mafft.fa

# check per sequence quality - to discard all outliers

###### H1N1 ###
### reference seq: EPI_ISL_390379; Guandong-Maonan; 1536; 2019; selected for vaccine 28.2.2020
### PB1
# how many seqs in the fasta file
count=$(grep -c ">" IAV_H1N1/H1N1_PB1_segm2_mafft.fa | cut -d" " -f1)
# read the reference sequence into 'Ref' array    
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 }' IAV_H1N1/H1N1_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
# for each sequence in the fasta file get a new array
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 } ' IAV_H1N1/H1N1_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H1N1/H1N1_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H1N1/H1N1_PB1_segm2_mafft.fa |  sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
# compare the reference with other sequences
# check for alternative bases (SNPs)
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
# count all gaps (dashes) and bases
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H1N1/H1N1_PB1_segm2_mafft_overview.csv
done

### PB2
count=$(grep -c ">" IAV_H1N1/H1N1_PB2_segm1_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 }' IAV_H1N1/H1N1_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 } ' IAV_H1N1/H1N1_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H1N1/H1N1_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H1N1/H1N1_PB2_segm1_mafft.fa |  sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H1N1/H1N1_PB2_segm1_mafft_overview.csv
done

### PA
count=$(grep -c ">" IAV_H1N1/H1N1_PA_segm3_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 }' IAV_H1N1/H1N1_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_390379\|/ { print ">"$0 } ' IAV_H1N1/H1N1_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H1N1/H1N1_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H1N1/H1N1_PA_segm3_mafft.fa |  sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H1N1/H1N1_PA_segm3_mafft_overview.csv
done

###### H3N2 ###
### reference seq: EPI_ISL_391201; Hong Kong; 2671; 2019; selected for vaccine 28.2.2020
# PB1
count=$(grep -c ">" IAV_H3N2/H3N2_PB1_segm2_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 }' IAV_H3N2/H3N2_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 } ' IAV_H3N2/H3N2_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H3N2/H3N2_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H3N2/H3N2_PB1_segm2_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H3N2/H3N2_PB1_segm2_mafft_overview.csv
done

# PB2
count=$(grep -c ">" IAV_H3N2/H3N2_PB2_segm1_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 }' IAV_H3N2/H3N2_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 } ' IAV_H3N2/H3N2_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H3N2/H3N2_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H3N2/H3N2_PB2_segm1_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H3N2/H3N2_PB2_segm1_mafft_overview.csv
done

# PA
count=$(grep -c ">" IAV_H3N2/H3N2_PA_segm3_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 }' IAV_H3N2/H3N2_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_391201\|/ { print ">"$0 } ' IAV_H3N2/H3N2_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"};  NR==ord { print ">"$0 }' IAV_H3N2/H3N2_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IAV_H3N2/H3N2_PA_segm3_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IAV_H3N2/H3N2_PA_segm3_mafft_overview.csv
done


###### IBV Victoria ###
# reference seq: EPI_ISL_362540; Washington; 02; 2019; selected for vaccine 28.2.2020
# PB1
count=$(grep -c ">" IBV_Victoria/Vict_PB1_segm2_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_362540\|/ { print ">"$0 }' IBV_Victoria/Vict_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_362540\|/ { print ">"$0 } ' IBV_Victoria/Vict_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; NR==ord { print ">"$0 }' IBV_Victoria/Vict_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Victoria/Vict_PB1_segm2_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Victoria/Vict_PB1_segm2_mafft_overview.csv
done
# PB2
count=$(grep -c ">" IBV_Victoria/Vict_PB2_segm1_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_362540\|/ { print ">"$0 }' IBV_Victoria/Vict_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_362540\|/ { print ">"$0 } ' IBV_Victoria/Vict_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; NR==ord { print ">"$0 }' IBV_Victoria/Vict_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Victoria/Vict_PB2_segm1_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Victoria/Vict_PB2_segm1_mafft_overview.csv
done

# PA
count=$(grep -c ">" IBV_Victoria/Vict_PA_segm3_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_3625401\|/ { print ">"$0 }' IBV_Victoria/Vict_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_362540\|/ { print ">"$0 } ' IBV_Victoria/Vict_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; NR==ord { print ">"$0 }' IBV_Victoria/Vict_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Victoria/Vict_PA_segm3_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Victoria/Vict_PA_segm3_mafft_overview.csv
done

###### IBV Yamagata ###
# reference seq: from 2013; here as the ref: EPI_ISL_296613; Kyiv,9, 2018
# PB1
count=$(grep -c ">" IBV_Yamagata/Yama_PB1_segm2_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 }' IBV_Yamagata/Yama_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 } ' IBV_Yamagata/Yama_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; NR==ord { print ">"$0 }' IBV_Yamagata/Yama_PB1_segm2_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Yamagata/Yama_PB1_segm2_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Yamagata/Yama_PB1_segm2_mafft_overview.csv
done

# PB2
count=$(grep -c ">" IBV_Yamagata/Yama_PB2_segm1_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 }' IBV_Yamagata/Yama_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 } ' IBV_Yamagata/Yama_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; !/EPI_ISL_2966130\|/ && NR==ord { print ">"$0 }' IBV_Yamagata/Yama_PB2_segm1_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Yamagata/Yama_PB2_segm1_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Yamagata/Yama_PB2_segm1_mafft_overview.csv
done

# PA
count=$(grep -c ">" IBV_Yamagata/Yama_PA_segm3_mafft.fa | cut -d" " -f1)
eval Ref=( $( awk 'BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 }' IBV_Yamagata/Yama_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
len=$(awk ' BEGIN{RS=">"}; /EPI_ISL_296613\|/ { print ">"$0 } ' IBV_Yamagata/Yama_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | wc -c | cut -d" " -f1 )
for ((i=2; i<$(($count+2)); i++))
do
eval Ar1=( $( awk -v ord=$i 'BEGIN{RS=">"}; NR==ord { print ">"$0 }' IBV_Yamagata/Yama_PA_segm3_mafft.fa | grep -v ">" | tr -d '\n' | awk -v len=$len 'BEGIN{FS=""}; { for (i=1; i<len+1; i++) { print $i }} ' ))
ID=$(grep ">" IBV_Yamagata/Yama_PA_segm3_mafft.fa | sed -n ''$i'p' | cut -d"|" -f2)
Das=0
As=0
Ts=0
Cs=0
Gs=0
Ns=0
asRef=0
for ((m=0; m<$len; m++))
do
if [[ ${Ar1[$m]} == ${Ref[$m]} ]]
then
asRef=$(($asRef+1))
fi
done
Das=$(echo ${Ar1[@]} | grep -o "-" | wc -l)
As=$(echo ${Ar1[@]} | grep -Po "[aA]" | wc -l)
Ts=$(echo ${Ar1[@]} | grep -Po "[tT]" | wc -l)
Cs=$(echo ${Ar1[@]} | grep -Po "[cC]" | wc -l)
Gs=$(echo ${Ar1[@]} | grep -Po "[gG]" | wc -l)
Ns=$(echo ${Ar1[@]} | grep -Po "[nN]" | wc -l)
echo -e $ID"\t"$Das"\t"$As"\t"$Ts"\t"$Cs"\t"$Gs"\t"$Ns"\t"$asRef"\t"$count >> IBV_Yamagata/Yama_PA_segm3_mafft_overview.csv
done

