# NGS
for NGS PJ

## План действий

1. Скачать хромосому 3 (chr3.fasta) - референс
2. FastQC - анализ качества 
3. Trimmomatic - тримминг 
4. BWA - картирование 
5. Samtools + bcftools - получить VCF 
6. Фильтрация SNP 
7. Написать отчет

# FastQC - анализ качества
(fastq2026) MacBook-Pro-Konstantin-2:~ konstantingeneralov$ cd PJ_NGS
(fastq2026) MacBook-Pro-Konstantin-2:PJ_NGS konstantingeneralov$ fastqc data.fastq
null
Started analysis of data.fastq
Approx 5% complete for data.fastq
Approx 10% complete for data.fastq
Approx 15% complete for data.fastq
Approx 20% complete for data.fastq
Approx 25% complete for data.fastq
Approx 30% complete for data.fastq
Approx 35% complete for data.fastq
Approx 40% complete for data.fastq
Approx 45% complete for data.fastq
Approx 50% complete for data.fastq
Approx 55% complete for data.fastq
Approx 60% complete for data.fastq
Approx 65% complete for data.fastq
Approx 70% complete for data.fastq
Approx 75% complete for data.fastq
Approx 80% complete for data.fastq
Approx 85% complete for data.fastq
Approx 90% complete for data.fastq
Approx 95% complete for data.fastq
Analysis complete for data.fastq
2026-04-21 16:42:39.744 java[13213:35307968] [JRSAppKitAWT markAppIsDaemon] failed. SetApplicationIsDaemon returned -50

