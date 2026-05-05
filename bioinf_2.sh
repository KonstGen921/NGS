#!/bin/bash

# Заголовок со ВСЕМИ колонками
echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAF\tMQ\tFQ\tINDEL\tIDV\tIMF\tVDB\tSGB\tRPB\tMQB\tMQSB\tBQB\tSCB\tMQ0F\tAF1\tAC1\tDP4\tPV4" > snp_table_2.txt

grep -v "^#" snp_filtered.vcf | while read line; do
  CHROM=$(echo "$line" | cut -f1)
  POS=$(echo "$line" | cut -f2)
  REF=$(echo "$line" | cut -f4)
  ALT=$(echo "$line" | cut -f5)
  QUAL=$(echo "$line" | cut -f6)
  INFO=$(echo "$line" | cut -f8)
  
  # Функция для извлечения значения из INFO
  get_value() {
    echo "$INFO" | grep -o "$1=[^ ;]*" | cut -d= -f2
  }
  
  # Извлекаем все параметры
  DP=$(get_value "DP")
  AF=$(get_value "AF")
  MQ=$(get_value "MQ")
  FQ=$(get_value "FQ")
  INDEL=$(get_value "INDEL")
  IDV=$(get_value "IDV")
  IMF=$(get_value "IMF")
  VDB=$(get_value "VDB")
  SGB=$(get_value "SGB")
  RPB=$(get_value "RPB")
  MQB=$(get_value "MQB")
  MQSB=$(get_value "MQSB")
  BQB=$(get_value "BQB")
  SCB=$(get_value "SCB")
  MQ0F=$(get_value "MQ0F")
  AF1=$(get_value "AF1")
  AC1=$(get_value "AC1")
  DP4=$(get_value "DP4")
  PV4=$(get_value "PV4")

  # Выводим ВСЕ колонки
  echo -e "$CHROM\t$POS\t$REF\t$ALT\t$QUAL\t$DP\t$AF\t$MQ\t$FQ\t$INDEL\t$IDV\t$IMF\t$VDB\t$SGB\t$RPB\t$MQB\t$MQSB\t$BQB\t$SCB\t$MQ0F\t$AF1\t$AC1\t$DP4\t$PV4" >> snp_table_2.txt
done

cat snp_table_2.txt