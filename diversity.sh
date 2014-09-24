# sh diversity.sh /home/javier/Desktop/virus h5tyturkey105_complete.bfa 20 60     

for file in `ls *1_001.fastq.gz`

do

echo $file

nb=`basename $file "_R1_001.fastq.gz"`
nf=$(echo $nb | cut -d '_' -f1)

echo "Creating directory "$nf

mkdir $1/$nf

gunzip -c $1/$nb"_R1_001.fastq.gz" > $1/$nf/$nf"_1.fastq"
gunzip -c $1/$nb"_R2_001.fastq.gz" > $1/$nf/$nf"_2.fastq"

cd $1/$nf

cp $1/$2 .

java -jar $1/trimmomatic-0.30.jar PE -phred33 $nf"_1.fastq" $nf"_2.fastq" $nf"_fp1.fastq" $nf"_f1.fastq" $nf"_fp2.fastq" $nf"_f2.fastq" CROP:127

python $1/cutTo127.py $1/$nf $nf"_fp1.fastq" $nf"_fp2.fastq" $nf"_127_fp1.fastq" $nf"_127_fp2.fastq"

maq fastq2bfq $nf"_127_fp1.fastq" $nf"_1.bfq"
maq fastq2bfq $nf"_127_fp2.fastq" $nf"_2.bfq"
maq match $nf".map" $2 $nf"_1.bfq" $nf"_2.bfq"
maq pileup -v $2 $nf".map" > $nf."pileup"

python $1/pileupCount.py $1/$nf $nf."pileup" $3 $4

rm -f *.fastq
rm -f *.bfq
rm -f *.bfa
rm -f *.map


done
