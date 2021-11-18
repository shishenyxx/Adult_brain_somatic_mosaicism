#! /bin/bash
## random simulation of each group for 10000 times based on number of variants detected

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..132};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_01_sorted_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..452};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_02_03_04_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..95};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_02_bulk_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..226};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_03_bulk_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..131};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_04_bulk_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..187};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_01_4tissue_bulk_simulation.txt;done

for i in {1..10000};do declare -i A;for i in {1..25};do A[i]="0";done;for i in {1..78};do i=$RANDOM%25+1;A[i]+=1;done;for i in {1..25};do echo ${A[i]};done|tr '\n' ','|awk '{print $0}'|sed -z 's/,\n/\n/g' >>ID_02_sorted_simulation.txt;done
