#! /bin/bash

cp ../*txt ./ #将上层文件夹的结果复制到当前文件夹

files_list=`ls -d err_prd*.txt` #输入文件
for CPTlim in '0.08' '0.5'           #色标范围
do
for infile in ${files_list[@]}  #循环
do
./ploterr.sh $infile 0 $CPTlim 0.1
done
done

files_list=`ls -d V_prd*.txt`  #输入文件
for CPTlim in '0.18'           #色标范围
do
for infile in ${files_list[@]}
do

./plot.sh $infile $CPTlim 0.06 

done
done

rm *.cpt *.grd gmt* clip

