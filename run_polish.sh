SPLIT="/media/k/SSDdata/220701splitsam"
PREVIOUS="/media/k/HDD4/220628hypohisat"
MAPPING="/media/k/HDD4/220701hypohisat"
REFFILE="220627yewgenome_polished.fasta"
RESULT="220701yewgenome_polished.fasta"
STORE="/media/k/HDD5/220608yewgenome_polishing"
SCORE="-0.3"
SPLITNUM="7"

LOGFILE=$SPLIT/log.txt

if [ ! -e $PREVIOUS/$REFFILE ];then 
echo previous reffile not found
exit
fi

echo polishing log $LOGFILE >> $LOGFILE
echo >> $LOGFILE
echo previous working folder: $PREVIOUS >> $LOGFILE
echo main mapping folder: $MAPPING >> $LOGFILE
echo split working folder: $SPLIT >> $LOGFILE
echo >> $LOGFILE
echo reference file: $PREVIOUS/$REFFILE >> $LOGFILE
echo result file: $MAPPING/$RESULT >> $LOGFILE
echo >> $LOGFILE

mkdir $MAPPING
cd $MAPPING

echo hisat2 start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE

hisat2-build -p 64 $PREVIOUS/$REFFILE index

hisat2 \
-x index \
-1 /media/k/HDD3/kusano/201221_novaseq_yew_unzip/abyssin_yewgenome_1.fastq \
-2 /media/k/HDD3/kusano/201221_novaseq_yew_unzip/abyssin_yewgenome_2.fastq \
-p 64 \
--rdg 1,1 \
--rfg 1,1 \
--score-min L,0,$SCORE \
--no-spliced-alignment -X 4000 \
--summary-file hisat2report.txt \
-S mapping.sam

rm index.*.ht2l

echo hisat2 finish >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE
ls -lh mapping.sam >> $LOGFILE
ls -l mapping.sam >> $LOGFILE
echo >> $LOGFILE

cd $SPLIT

echo splitsam start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE

perl splitsam.pl $MAPPING/mapping.sam $SPLITNUM 100
rm $MAPPING/mapping.sam

echo splitsam end >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE
echo >> $LOGFILE

cd $MAPPING
cp $PREVIOUS/il_names.txt .

SEQNUM=$(($SPLITNUM-1))
for i in `seq 0 $SEQNUM`
do
echo samtools $i start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE
samtools view -@64 -bo splitbam_$i.bam $SPLIT/splitsam_$i.sam
samtools sort -@64 -o splitbam_sorted$i.bam splitbam_$i.bam
samtools index -@64 splitbam_sorted$i.bam
echo hypo $i start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE
hypo -d $PREVIOUS/splithypo_$i.fa \
-r @il_names.txt -s 4g -c 30 \
-b splitbam_sorted$i.bam -t 64 -o splithypo_$i.fa
done

echo hypo end >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE
echo >> $LOGFILE

rm splitbam*
rm -r aux
rm $SPLIT/splitsam_*.sam

echo make result file start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE

for i in `seq 0 $SEQNUM`
do
cat splithypo_$i.fa >> $RESULT
done
ls -lh $RESULT >> $LOGFILE
ls -l $RESULT >> $LOGFILE

if [ ! -e $RESULT ];then
echo result file not found error >> $LOGFILE
fi

if [ -e $RESULT ]; then

echo delete old files start >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE

cd $PREVIOUS
rm splithypo_*.fa
cp $REFFILE $STORE
rm $REFFILE

echo done >> $LOGFILE
echo `date '+%y/%m/%d %H:%M:%S'` >> $LOGFILE

fi

