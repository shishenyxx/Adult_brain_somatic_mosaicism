##code for collecting insertion size using GATK

cat input_files.txt |cut -f1,4|grep -v tumor|while read name bam;do echo "java -jar -Xmx14G -Djava.io.tmpdir=temp picard.jar CollectInsertSizeMetrics I="${bam}" O=/Collect_insertion_size/"${name}".insert_size_metrics.txt H=/Collect_insertion_size/"${name}".insert_size_histogram.pdf M=0.5"|qsub -q home -l nodes=1:ppn=4,mem=15g -l walltime=48:00:00;done
