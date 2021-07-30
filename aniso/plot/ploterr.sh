#! /bin/bash

gmt set MAP_FRAME_TYPE=plain
gmt set PS_MEDIA=a4  #图纸大小 
R=99.6/106/25.9/32.5 #画图区域：最小经度/最大经度/最小纬度/最大纬度
J=M102.5/29/7i       #图形中心经度/图形中心纬度/图形大小
FT=CN-faults.dat     #断层数据
BD1=block.dat		 #板块边界数据1
BD2=block_add.dat	 #板块边界数据2


infile=$1
v1=$2
v2=$3
ban=$4

echo "$infile"
gawk '{print $2,$1,$3}' $infile > tmp.tmp
pp1=`echo $infile | awk -F '.txt' '{print $1}'`
PS=${pp1}.ps

echo $v1 $v2
gmt psxy -R$R -J$J -T -Ba01WesN -P -K  > $PS
gmt makecpt -Cseis -T$v1/$v2/0.01 -Z > velocity.cpt
#makecpt -Cjet -T$v2/$v1/0.05 -Z> velocity.cpt

line=$(grep -n "B" velocity.cpt | cut -d ":" -f 1)
rlp=`echo "B 170/0/0"`
sed -i "${line}c $rlp" velocity.cpt

rlp=`echo "F mediumblue"`
line=$(grep -n "F" velocity.cpt | cut -d ":" -f 1)
sed -i "${line}c $rlp" velocity.cpt

if [ -f velocity.grd ]
then
	rm velocity.grd
fi 

gmt surface tmp.tmp -R -I10m -Gvelocity.grd
########剪切区域
gmt psmask tmp.tmp -I5m -R -J -F -G50 -S20m -Dclip -V -O -K >> $PS
gmt psclip clip -R -J -B -K -O>> $PS
gmt grdview velocity.grd -Gvelocity.grd -Cvelocity.cpt -R -J -B -Qi720 -K -O >>$PS
gmt psmask -C -K -O >> $PS

echo "plot faults..."
gmt psxy $FT -R$R -J$J  -V -W0.1p,blue -K -O  >> $PS
gmt psxy $BD1  -R$R -J$J  -V -W0.8p,black -K -O  >> $PS
gmt psxy $BD2  -R$R -J$J  -V -W0.8p,black -K -O  >> $PS

echo "plot err"
cat $infile | gawk '{print $2,$1,0,$5}' | gmt psxy -R$R -J$J -SW0.5c -Gwhite -K -O >> $PS
cat $infile | gawk '{print $2,$1,0,0.15*$4}' | gmt psxy -R$R -J$J -SV0.5p -W1p,black -Gblack -V -K -O >> $PS

echo "plot err"
echo 104.25 26.25 0 0 | gmt psxy -R$R -J$J -SW0.25c -Ggray -K -O >> $PS
echo 104.5 26.5 0 0 | gmt psxy -R$R -J$J -SW0.5c -Ggray -K -O >> $PS
echo 104.75 26.75 0 0 | gmt psxy -R$R -J$J -SW0.75c -Ggray -K -O >> $PS
echo 105 26.5 90 0.3  | gmt psxy -R$R -J$J -SV0.5p -W2p,black -Gblack -V -K -O >> $PS
echo 104.7 26.5 "1% " | gmt pstext -J$J -R$R -F+f14p,2,black -K -O >> $PS
gmt psscale -Y-0.4c -Cvelocity.cpt -D1.3i/0i/2.5i/0.15ih -B$ban:"":/:"km/s": -S -O -K >>$PS
gmt psxy -R$R -J$J -T -O >> $PS
ps2pdf $PS ${pp1}_$v2.pdf
rm *.cpt *.grd gmt* clip temp.temp *.ps



