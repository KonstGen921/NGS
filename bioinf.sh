#!/bin/bash

echo "CHROM	POS	REF	ALT	QUAL	DP	INFO" > snp_table.txt

grep -v "^#" snp_filtered.vcf | while read line; do
  CHROM=$(echo "$line" | cut -f1)
  POS=$(echo "$line" | cut -f2)
  REF=$(echo "$line" | cut -f4)
  ALT=$(echo "$line" | cut -f5)
  QUAL=$(echo "$line" | cut -f6)
  INFO=$(echo "$line" | cut -f8)
  
  # Извлекаем DP из INFO
  DP=$(echo "$INFO" | grep -o "DP=[0-9]*" | cut -d= -f2)
  
  echo -e "$CHROM\t$POS\t$REF\t$ALT\t$QUAL\t$DP\t$INFO" >> snp_table.txt
done

cat snp_table.txt