%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   - by Liang Chuantao, 2009                  %%%%%%%%%%%%%%%
%%%                                              %%%%%%%%%%%%%%%
%%%   - liangct@cdut.edu.cn                      %%%%%%%%%%%%%%%
%%%   - before running this program, also check   %%%%%%%%%%%%%%%
%%%     all parameters.                          %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function GradAnaReg
clc
clear;
close all;
if strcmp(computer,'GLNX86')
    fs='/';
else
    fs='\';
end
pathroot='.\SACdata';
datout='.\WG2';
wparamters=[pathroot fs 'evg.txt'];
cmp='BHZ';
sacftyp=['T1*' cmp '.sac'];
par.dt=0.1;
par.reg=[25.8 32.5 ; 99.8 105.2;];%range of the longitude and latitude of array
par.dsmax=50;  %dsmax: the maximum distance of station
par.dsmin=5;  %dsmin: the minimum distance of station

par.datout=datout;
par.kwvmodel='sph';
kphase='R'; %%P/R/S/L
par.pathroot=pathroot;
par.cmp=cmp;
par.sacftyp=sacftyp;
par.omark='nz00pct.';
par.tsac=[];   %%%T0-9 .or. [];
par.knz=0; %noise level: 0.0, no noise added to the real data
par.kredu='vrphs'; %% 'vr???/vrphs/vrgrp/vrnan';
par.dvmax=0.025;
par.kazm=0;  %%%always be 0, other values for test purpose
par.kgauss=[]; %%don't change this parameter in this file
par.kwi=1;  %%if kwi<=0, no weight applied
par.fs=fs;

[evfdd,vg,pband(:,1),pband(:,2)]=textread(wparamters,'%s %f %f %f');
% loops for all earthquakes 
for i=264:length(evfdd)
    par.pband=pband(i,:);
    par.evtime=evfdd(i); 
    par.vg=vg(i);
    par.vredu=par.vg+0.2; %%%starting reducing velocity
    pstr=[kphase '.p' num2str2(mean(par.pband),3,0)];
    par.tprd=pstr;
    disp(['loops=' num2str(i) ', ev=' cell2mat(par.evtime) ', pb=' num2str(par.pband)]);
    GradAna1eq(par);
end

%-------------------------------------------------------
%   WG for one earthquake
%-------------------------------------------------------
function GradAna1eq(par)
kwvmodel=par.kwvmodel;
sacftyp=par.sacftyp;
fs=par.fs;
kwi=par.kwi;
kdat=cell2mat(par.evtime);
reg=par.reg;
pband=par.pband;
tprd=par.tprd;
dseis=[par.pathroot fs ];                       
cmp=par.cmp;

fwin=[1/pband(2) 1/pband(1)];
dout=[par.datout fs kdat fs];
dseis=[dseis kdat fs ];
omark=par.omark; 

if ~exist(dout,'dir')
    mkdir(dout)
end

if kwi==-1
    fout=[dout 'wga.' tprd '.' cmp '.' kdat '.wtn.' par.kredu omark 'txt'];
elseif kwi==0
    fout=[dout 'wga.' tprd '.' cmp '.' kdat '.wt0.' par.kredu omark 'txt'];
elseif kwi==1
    fout=[dout 'wga.' tprd '.' cmp '.' kdat '.wt1.' par.kredu omark 'txt'];
elseif kwi==2
    fout=[dout 'wga.' tprd '.' cmp '.' kdat '.wt2.' par.kredu omark 'txt'];
end
par.kwi=kwi;
par.ermax=[];
par.clr='b';
par.wlen=max(pband)/2;
par.vavg=par.vredu;
par.dvg=0.5;
par.kdat=kdat;

%--find all stations in the given region--
sta=allsta(dseis,sacftyp,reg);
nst=length(sta);
fido=fopen(fout,'w');

for ist=1:nst
    %--find stations around the central element--
    stn=sta(ist).stn;
    arr=ChsSta(sta,ist,par);
    while length(arr)<5 && par.dsmax<80
        par.dsmax=par.dsmax+10;
        arr=ChsSta(sta,ist,par);
    end
    st=[];
    if ~isempty(arr) && length(arr) >= 3
        [st,ev]=readwaves(dseis,arr);
        %%%add random errors:
        if par.knz>0
            st=addnoise(st,par.knz);
        end
        %[wf]=gaussfilt(st(1).dt,0.5,st(1).dis,mean(pband),st(1).dat,kplot);
        wf=bpfilt(st(1).dat,st(1).dt,1/pband(2),1/pband(1));
        vgmax=[];
        if ~isempty(par.tsac)
            vgmax=find_tpick(par.tsac,st(1).tphs,st(1).dis);
        else
            vgmax=find_vgmax(wf,st(1).dis,st(1).tbeg,st(1).dt,par.vg,par.dvg);
        end
        if ~isempty(vgmax) &&  ~isnan(vgmax)
            %narrow band filter is applied in stFilter
            st=stFilter(st,par.tsac,vgmax,pband,200);
            par.vgmax=vgmax;
            if ~isempty(st) && length(st)>=3
                %---match times---
                [st]=timmatch(st,par.dt);  
                %----WG for one station----
                GradAna0(st,ev,par,fwin,fido,kwvmodel,par.kredu) 
            end
        else
            warning('...vgmax not set!!');
        end 
    end
end
% end
fclose(fido);

%%-------------------------------------------
%%      add noise                           %
%%-------------------------------------------
function st=addnoise(st,knz)
if isempty(knz)
    return
end

nst=length(st);


for ii=1:nst
    figure(1)
    plot(st(ii).dat)
    npt=length(st(ii).dat);
    wmax=max(abs(st(ii).dat));
    wfnz=rand(1,npt)*2-1.0;
    wfnz=wfnz*(knz*wmax);
    st(ii).dat=st(ii).dat+wfnz;
    figure(2)

    plot(st(ii).dat)
end



%%
%%
%%
function vpk=find_tpick(tsac,tpk,dis)
ktim=str2num(tsac(2))+1;
tpk=tpk(ktim);
vpk=dis/tpk;

%
% station selections
%
function st=stFilter(st,tsac,vgmax,pband,twin)

nst=length(st);
amax(1)=100000;
for ii=1:nst

    dt=st(ii).dt;
    tg=st(ii).dis/vgmax-st(ii).tbeg;
    tw=[tg-twin tg+twin];
    nt=ceil(tw/dt);
    if nt(1)<1
        nt(1)=1;
    end
    if nt(2)>length(st(ii).dat)
        nt(2)=length(st(ii).dat);
    end
    nt=[nt(1):nt(2)];
    if ~isempty(nt)
        %[wf]=gaussfilt(st(ii).dt,0.5,st(ii).dis,mean(pband),st(ii).dat,0);

        wf=bpfilt(st(ii).dat,dt,1/pband(2),1/pband(1));
        st(ii).dat=wf;
        amax(ii)=max(abs(wf(nt)));
    else
        amax(ii)=amax(1)*100;
    end

    %%%If no value picked: discard
    if ~isempty(tsac)
        ktim=str2num(tsac(2))+1;
        tpk=st(ii).tphs(ktim);
        if isnan(tpk)
            amax(ii)=amax(1)*100;
        end
    end
end
ist=[1:nst];
amax=amax/amax(1);
%amplitude variation should be less than 20%
ist(amax>1.3)=[];
amax=amax(ist);
ist(amax<0.7)=[];
st=st(ist);


%
%--Find all station in the region
%
function st=allsta(dseis,sacftyp,reg)
sti=[dseis sacftyp];
fl=dir(sti);
nst=length(fl);
kst=0;
for ist=1:nst
    fli=fl(ist).name;    
    sti=readseis(dseis,fli,'');  %read 'EW' component Seismogram 
    if ~isempty(sti)
        st0=sti.st;
        if st0(1)>=reg(1) & st0(1,1)<=reg(1,2) & st0(2)>=reg(2,1) & st0(2)<=reg(2,2)        
            kst=kst+1;
            sti.fl=fli;
            st(kst)=sti;
        end
    end      
end


%
%--Search neighboring stations---
%
function stn=ChsSta(st,kst,par)
%--calculate distances and azimuths between stc and other stations
stc=st(kst);
nst=length(st);
nst=[1:nst];
nst(kst)=[];
for ist=1:length(nst);
    ii=nst(ist);
    sti=st(ii);
    azm(ist)= azimuth0(stc.st(1),stc.st(2),sti.st(1),sti.st(2));
    dst(ist)=distance0(stc.st(1),stc.st(2),sti.st(1),sti.st(2))*111.199;
end

idx=nst(dst>=par.dsmin & dst<=par.dsmax);
npth=length(idx);
if npth<2
    %stn=[];
    %msg='path too few'
    %return
end
stn(1)=stc;
stn(2:npth+1)=st(idx);

%
%--read seismogram-----------------------
%
function dat=readseis(dir0,fl,kdat)
if ~exist([dir0 fl],'file')
    dat=[];
    return
end
[A B C D]=readsac0([dir0 fl]);
%%%%%%%%%%%%%%%
%if isnan(A(4))
 %return
%end
%%%%%%%%%%%%%%
dat.st(1:3)=A(32:34);
dat.stn=cell2mat(C(1));
dd=B(17);
if ~isempty(kdat) && kdat==1
    dat.dt=A(1);
    %dat.scale=A(4);
    dat.scale=1;
    dat.ev=[A(36) A(37) A(39)];
    dat.dat=D'/dat.scale;
    tbeg=[A(6)];
    tend=[A(7)];
    dat.dis=A(51);
    dat.azm=A(52);

    %t0=[str2num(sdate(1:4)) str2num(sdate(6:8)) str2num(sdate(10:11)) str2num(sdate(13:14)) str2num(sdate(16:22))];
    t0=[B(1:5)];  %start time of the time series
    t0(5)=t0(5)+B(6)/1000+tbeg;
    dat.t0=t0;
    t1=t0;
    t1(5)=t1(5)+tend-tbeg;
    dat.t1=t1;
    dat.tphs=A(11:20);
    dat.tbeg=tbeg;
end



%
%
%
function [st,ev]=readwaves(dseis,arr)
nst=length(arr);
kst=0;
for ist=1:nst
    fl=arr(ist).fl;    
    sti=readseis(dseis,fl,1);         %read 'EW' component Seismogram
    if ~isempty(sti)
        kst=kst+1;
        st(kst)=sti;
        ev=st(kst).ev;
    end
end


%-------------------------------------------------
%---match the original time and sample rate-------
%-------------------------------------------------
function [st]=timmatch(st,dtm)
nst=length(st);
for is=1:nst
    tm2=st(is).t0;
    t0(is)=tm2(2)*24*3600+tm2(3)*3600+tm2(4)*60+tm2(5);
    dt(is)=st(is).dt;
    t1(is)=t0(is)+dt(is)*length(st(is).dat)-dt(is);
end

dt=round(dt*100000)/100000;
if isempty(dtm)
    %--Find the sampling rate that is integer times of all sampling rates--
    nn=1;
    dr=mod(max(dt)*nn./dt,1.0);
    while ~isempty(find(dr~=0.0))
        nn=nn+1;
        dr=mod(max(dt)*nn./dt,1.0);
    end
    dtm=max(dt)*nn;
end

%--decimate all time series----
for is=1:nst
    r=dtm./dt(is);
    if r>1
        fnq=1.0/2/dtm;
        st(is).dat=decimate(st(is).dat,r);
        st(is).dt=dtm;
        dt(is)=dtm;
    end
end

%--match all time series with same the same starting and ending time----
t0max=max(t0);
t1min=min(t1);
tlen=(t1min-t0max);
if tlen<0
    st=[];
    return
end
nt0=round((t0max-t0)./dt)+1;
nt1=nt0+round(tlen./dt)-1;

for is=1:nst
    n0=nt0(is);
    n1=nt1(is);
    st(is).dat=st(is).dat(n0:n1);
    st(is).t0(5)=st(is).t0(5)+t0max-t0(is);
    st(is).tbeg=st(is).tbeg+t0max-t0(is);
end
