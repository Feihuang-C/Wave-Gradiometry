#! /bin/bash
gmt set MAP_FRAME_TYPE=plain
gmt set PS_MEDIA=a4  #图纸大小 
R=99.6/106/25.9/32.5 #画图区域：最小经度/最大经度/最小纬度/最大纬度
J=M102.5/29/7i       #图形中心经度/图形中心纬度/图形大小
FT=CN-faults.dat     #断层数据
BD1=block.dat		 #板块边界数据1
BD2=block_add.dat	 #板块边界数据2

infile=$1
b=$2
ban=$3

echo "$infile"

#find mean velocity
gawk '{print $2,$1,$3}' $infile > tmp.tmp
cac=`awk '{v[NR]=$3;sum+=$3}{if($3>max)max=$3;if(min>$3)min=$3}NR==1{max=min=$3}END{avg=sum/NR;{dv=max-avg-0.0};printf("Max=%f\nMin=%f\nAvg=%f\ndv=%f\n",max,min,avg,dv)}' tmp.tmp`
vavg=`echo $cac | awk -F ' ' '{print $3}'| awk -F '=' '{print $2}'`
v1=`echo "$vavg - $b"|bc`
v2=`echo "$vavg + $b"|bc`

#################
pp1=`echo $infile | awk -F '.dat' '{print $1}'`
PS=${pp1}.ps

echo $v1 $v2
#########################################################################
gmt psxy -R$R -J$J -T -Ba01WesN -P -K  > $PS
gmt makecpt -Cseis -T$v1/$v2/0.01 -Z > velocity.cpt

line=$(grep -n "B" velocity.cpt | cut -d ":" -f 1)
rlp=`echo "B 170/0/0"`
sed -i "${line}c $rlp" velocity.cpt

rlp=`echo "F mediumblue"`
line=$(grep -n "F" velocity.cpt | cut -d ":" -f 1)
sed -i "${line}c $rlp" velocity.cpt

########################################################################
if [ -f velocity.grd ]
then
	rm velocity.grd
fi 
gmt surface tmp.tmp -R -I5m -Gvelocity.grd

########
gmt psmask tmp.tmp -I5m -R -J -F -G50 -S20m -Dclip -V -O -K >> $PS
gmt psclip clip -R -J -B -K -O>> $PS
gmt grdview velocity.grd -Gvelocity.grd -Cvelocity.cpt -R -J -B -Qi720 -K -O >>$PS
gmt psmask -C -K -O >> $PS

echo "plot faults..."
gmt psxy $FT -R$R -J$J  -V -W0.1p,blue -K -O  >> $PS
gmt psxy $BD1  -R$R -J$J  -V -W0.8p,black -K -O  >> $PS
gmt psxy $BD2  -R$R -J$J  -V -W0.8p,black -K -O  >> $PS


#plot aniso
cat $infile | gawk '{print $2,$1,$5,0.25*$4}' | gmt psxy -R$R -J$J -SV3p -W1p,black -Gblack -V -K -O >> $PS
echo 104.8 26.25 90 0.5 |gmt psxy -R$R -J$J -SV4p -W2p,black -Gblack -V -K -O >> $PS
echo 104.8 26.5 90 1 |gmt psxy -R$R -J$J -SV4p -W2p,black -Gblack -V -K -O >> $PS
echo 104.6 26.25 "2% " | gmt pstext -J$J -R$R -F+f18p,2,black -K -O >> $PS
echo 104.6 26.5 "4% " | gmt pstext -J$J -R$R -F+f18p,2,black -K -O >> $PS

gmt psscale -Y-0.4c -Cvelocity.cpt -D1.3i/0i/2.5i/0.15ih -B$ban:"":/:"km/s": -S -O -K >>$PS
echo "end ..."
gmt psxy -R$R -J$J -T -O >> $PS
ps2pdf $PS ${pp1}_$b.pdf
rm *.cpt *.grd gmt* clip tmp.tmp *.ps

