%%
%% --------------Wave Gradiometry-----------------
%%
function GradAna0(st,ev,par,fwin,fido,kwvmod,kredu)         
nst=length(st);
stc=st(1);
wfc=stc.dat;
tbeg=stc.tbeg;
dt=stc.dt;
dis0=distance0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))*111.199;
azm0=azimuth0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))+180; %st(1)台站经度，st(2)台站维度；ev（1）事件经度
if azm0>360
    azm0=azm0-360;
end
nsti=length(st);
    %find the group velocity
    %vgmax=find_vgmax(wfc,dis0,tbeg,dt,par.vg,par.dvg);
    vgmax=par.vgmax;
    if strcmp(kredu,'vrphs')
            %find the global phase velocity
        dv=100;
        nloop=0;
        vphs0=par.vredu;
        while dv>par.dvmax
            nloop=nloop+1;
            if nloop>25
                return
            end
            [ab vphs stc]=GradAna8redu(st,ev,par,fwin,fido,kwvmod,kredu);
            dv0=dv;
            dv=abs(vphs-par.vredu);
            if isnan(vphs)||median([st.dis])/vphs>length(st(1).dat)
                return
            end
            par.vredu=vphs;
        end
        nloop=nloop;
    elseif strcmp(kredu,'vrgrp')
        par.vredu=vgmax;
        [ab vphs stc]=GradAna8redu(st,ev,par,fwin,fido,kwvmod,kredu)
    elseif strcmp(kredu,'vrnan')
        par.vredu=-1;
        [ab vphs stc]=GradAna8redu(st,ev,par,fwin,fido,kwvmod,kredu)
    else
        par.vredu=par.vredu;
        [ab vphs stc]=GradAna8redu(st,ev,par,fwin,fido,kwvmod,kredu)
    end
    OutPutVelAzm(fido,ab,azm0,stc,nsti,dt,vgmax,par.vg,par.wlen,dis0)

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WG: par.vredu >0 <0, w/o reducing velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ab vphs stc]=GradAna8redu(st,ev,par,fwin,fido,kwvmod,kredu)         
nst=length(st);
stc=st(1);
wfc=stc.dat;
tbeg=stc.tbeg;
dt=stc.dt;
dis0=distance0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))*111.199;

    vgmax=par.vgmax;
    if par.vredu>0
        [st]=timeshift(st,ev,par.vredu,dt);
    end

    [dudx,dudy,stc]=WGA_dxdy(st,ev,fwin,par);
    stc.dst=distance0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))*111.199;
    nsti=length(st);
    clear st;
    %[ab,wfc,dt]=wga_spc(stc,dudx,dudy,fwin);
    [ab,wfc,dt]=WGA_AB(stc,dudx,dudy,par.kazm,kwvmod);
    azm0=azimuth0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))+180;
    if azm0>360
        azm0=azm0-360;
    end
    ab=WGA_PhyParam(par.vredu,ab,azm0,stc.dst,par.kazm);
    vphs=find_vphase(ab,stc,dt,vgmax,par.wlen);
    

%%-----------------------------------------------
%% add back the effects of the reducing velocity\
%%-----------------------------------------------
function ab=WGA_PhyParam(vredu,ab,azm0,dst,kazm)
    if vredu>0
        predu=1/vredu;
        ab.Bx=ab.Bx-predu*sin(azm0/180*pi);
        ab.By=ab.By-predu*cos(azm0/180*pi);
    end
            %-----Azimuth Estimate---------------
        if kazm==0
            ab.theta=atan2(-ab.Bx,-ab.By);
        else
            ab.theta=atan2(ab.dudx,ab.dudy);
        end
            %-----Ray Parameter Estimate-----------
        ab.pr=sqrt(ab.Bx.*ab.Bx+ab.By.*ab.By);
        am=azm0/180*pi;
        ab.radiation=(ab.Ax.*cos(am)-ab.Ay.*sin(am))*dst;
        ab.gsprd=(ab.Ax.*sin(am)+ab.Ay.*cos(am));
    

%
%       -----reduce time-----
%
function [st]=timeshift(st,ev,vredu,dt)
nst=length(st);
%shift the waveforms.
for is=1:nst
    dis0=distance0(st(is).st(1),st(is).st(2),ev(1),ev(2))*111.199;
    tsh=dis0/vredu;
    nts(is)=ceil(tsh/dt);
end
nts=nts-min(nts);

%make sure all waveforms are in same length
[nmax idx]=max(nts);
np=length(st(1).dat)-nmax;
for is=1:nst
    st(is).tsh=nts(is)*dt;
    st(is).tbeg=st(is).tbeg+st(is).tsh;
    if nts(is)>0
        st(is).dat(1:nts(is))=[];
    end
    npi=length(st(is).dat);
    if npi>np
        st(is).dat(np+1:npi)=[];
    end
end



%-----------------------------------
%------Output results-------------
%-----------------------------------
function vphs=find_vphase(ab,stc,dt,vg,twin)
kmax=ceil((stc.dst/vg-stc.tbeg)/dt); %--point number of estimated velocity

nt=round(twin/dt);                  %time window to be averaged.
nt=[-nt:1:nt]+kmax;
nt(nt<1)=[];
nt(nt>length(ab.pr))=[];

pr=ab.pr(nt);                      %ray-parameter in the time window
vapp=1.0./pr;                 %velocity in the average-window
vphs=mean(vapp);

return
%---test----
subplot(2,1,1)
plot(1.0./ab.pr)
hold on;
plot(nt,vapp,'r')
ylim([0 10])

%-----------------------------------
%------Output results-------------
%-----------------------------------
function OutPutVelAzm(fido,ab,azm0,stc,nsti,dt,vgpick,vg0,twin,dis0)
stn=stc.stn;

st=stc.st;

[envu,phau]=envelope(stc.dat);
latitude=stc.ev(1);
longitude=stc.ev(2);
depth=stc.ev(3);
tpeak=stc.dst/vgpick-stc.tbeg;
nt=round(tpeak/dt);
dnt=round(twin/dt);
nt1=nt-dnt;
nt2=nt+dnt;

if nt2>length(ab.pr); nt2=length(ab.pr); end
if nt1<1; nt1=1; end
        
if nt2>nt1
    amax=max(envu(nt1:1:nt2));
else
    amax=[];
    return
end


qual=1;

%---apparant velocity,azimuth and ---
vel=1.0./ab.pr(nt1:nt2);
azmo=ab.theta(nt1:nt2)*180/pi;
azm=azmo-azm0;
azm(azm<-180)=azm(azm<-180)+360;
azm(azm>180)=azm(azm>180)-360;
gsp=ab.gsprd(nt1:nt2)*1000;          %Geometrical Spreading * 1000
rad=ab.radiation(nt1:nt2);          %Geometrical Spreading  * 1000

azmo=mean(azmo);
v0=mean(vel);
dv=sqrt(var(vel));
a0=mean(azm);
da=sqrt(var(azm));
g0=mean(gsp);
dg=sqrt(var(gsp));
r0=mean(rad);
dr=sqrt(var(rad));
a000=azm0;
dis=dis0;

fprintf(fido,'%7.3f %8.3f %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %8.3f %2g %2g %5s %8.3f %7.3f %8.3f %7.3f %7.3f\n', ...
    st(1),st(2),st(3),v0,dv,a0,da,r0,dr,g0,dg,vgpick,azmo,nsti,amax,stn,a000,latitude,longitude,depth,dis);
msg=['station:' stn ' saved!!'];

return
figure(22)
subplot(3,1,1)
plot(stc.dat);
subplot(3,1,2)
plot(1.0./ab.pr);
ylim([0 7])
subplot(3,1,1)
plot(ab.theta);


%
%--Find all station in the region
%
function st=allsta(dseis,ev,fs,reg,cmp)
s0=[num2str(ev(1)) num2str2(ev(2),2,0) ...
    num2str2(ev(3),2,0)  num2str2(ev(4),2,0)];
dseis=[dseis s0 fs];

sti=[dseis s0 '*' cmp '*'];
fl=dir(sti);
nst=length(fl);

kst=0;
for ist=1:nst
    fli=fl(ist).name;    
    sti=readseis(dseis,fli,'');  %read 'EW' component Seismogram        
    st0=sti.st;
    if st0(1)>=reg(1) & st0(1,1)<=reg(1,2) & st0(2)>=reg(2,1) & st0(2)<=reg(2,2)        
        kst=kst+1;
        sti.fl=fli;
        st(kst)=sti;
    end
end


%
%--read seismogram-----------------------
%
function dat=readseis(dir0,fl,kdat)
if ~exist([dir0 fl],'file')
    dat=[];
    return
end
[A B C D]=readsac0([dir0 fl]);

dat.st(1:3)=A(32:34);
dat.stn=cell2mat(C(1));
dd=B(17);
if ~isempty(kdat) & kdat==1
    dat.dt=A(1);
    dat.scale=A(4);
    dat.ev=[A(36) A(37) A(39)];
    dat.dat=D';
    tbeg=[A(6)];
    tend=[A(7)];

    %t0=[str2num(sdate(1:4)) str2num(sdate(6:8)) str2num(sdate(10:11)) str2num(sdate(13:14)) str2num(sdate(16:22))];
    t0=[B(1:5)];  %start time of the time series
    t0(5)=t0(5)+B(6)/1000+tbeg;
    dat.t0=t0;
    t1=t0;
    t1(5)=t1(5)+tend-tbeg;
    dat.t1=t1;
    dat.tphs=A(12:17);
    dat.tbeg=tbeg;
end

if isnan(A(4))
    dat=[];
end


%
%
%
function [st,ev]=readwaves(dseis,ev,fs,arr)
damp=0.01;
s0=[num2str(ev(1)) num2str2(ev(2),2,0) ...
    num2str2(ev(3),2,0)  num2str2(ev(4),2,0)];
dseis=[dseis s0 fs];
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


%-----------------------------------
%-------calculate ray parameter-----
%-----------------------------------
function [pr]=plotVelAzm(ab,azm0,stc,dt,twin,par,ifig,ifig2)
wfc=stc.dat;
tp=stc.tphs(1)-stc.tbeg;
ts=stc.tphs(2)-stc.tbeg;
np=length(wfc);
nh=np/2;
tim=[0:np-1]*dt;
[emax,kmax]=max(abs(hilbert(wfc)));

px=ab.Bx;
py=ab.By;

pr=ab.pr;
%pr=mean2(pr,2);

clr=par.clr;
if ~isempty(ifig2)
    figure(ifig2)
    hold on;
    plot(tim,1.0./pr);
    text(tim(kmax),1.0/pr(kmax),num2str(ifig))
    %xlim([twin]);
    ylim([0 20]);
    xlabel('Time (sec)');
    ylabel('Apprarent Velocity');
    title('Waveform ');
    grid on;
    return
end

figure(ifig);
subplot(3,1,1)
plot(tim,wfc,'color',clr);
hold on;
plot(tim,abs(hilbert(wfc)),'r');
plot([tp tp],[min(wfc) max(wfc)],'r');
plot([ts ts],[min(wfc) max(wfc)],'r');
xlim(twin);
xlabel('Time (sec)');
ylabel('Amplitude');
title(par.ttl);
grid on;

subplot(3,1,2)
hold on;
%set(gca,'xtick',[0:200:1000])
plot(tim,1.0./pr,'color',clr);
xlim([twin]);
ylim([0 25]);
xlabel('Time (sec)');
ylabel('Apprarent Velocity');
title('Waveform ');
grid on;

subplot(3,1,3)
theta=ab.theta*180/pi;
plot(tim,theta,'.','color',clr);
hold on;
grid on;

plot([tim(1) tim(np)],[azm0 azm0],'g');
plot([tim(1) tim(np)],-[azm0 azm0],'g');
%plot(tim,abs(hilbert(theta)),'g');
title('Azimuth');
%dang=acal-avg
xlim(twin);
ylim([0 180]);
xlabel('Time (sec)');
ylabel('Azimuth');
grid on;

%------------------------------------------------
%-------calculate Ax Bx Ay and By values---------
%------------------------------------------------
function [AB,Dc,dt]=WGA_AB(stc,dudx,dudy,kazm,kmodel)
Dc=stc.dat;
dt=stc.dt;
np=length(Dc);
nh=np/2;
df=1.0/dt/np;

[envdx,phadx]=envelope(dudx);
[envdy,phady]=envelope(dudy);
[envu,phau]=envelope(Dc);

wt=InsantFrq(Dc,dt);
envudt=diff(envu)/dt;
envudt(np)=envudt(np-1);

if strfind(kmodel,'pla')
    AB.kmodel=0;
    %Bx=envdx./envu./wt.*sin(phadx-phau);
    %By=envdy./envu./wt.*sin(phady-phau);

    dudt=diff(Dc)/dt;
    dudt(np)=dudt(np-1);
    [envdudt,phadudt]=envelope(dudt);
    Bx=envdx./envdudt;
    By=envdy./envdudt;

    Ax=0.0;
    Ay=0.0;

    AB.Ax=Ax;
    AB.Ay=Ay;
    AB.Bx=Bx;
    AB.By=By;

    AB.xc=stc.st(2);
    AB.yc=stc.st(1);
        %-----Azimuth Estimate---------------
    if kazm~=0
        theta=atan2(dudx,dudy);
    else
        theta=atan2(-AB.Bx,-AB.By);
    end
    AB.theta=theta;
        %-----Ray Parameter Estimate-----------
    %AB.pr=-AB.Bx.*sin((theta))-AB.By.*cos((theta));
    AB.pr=sqrt(AB.Bx.*AB.Bx+AB.By.*AB.By);
    %AB.radiation=(AB.Ax.*cos(theta)-AB.Ay.*sin(theta))*stc.dst;
    %AB.gsprd=(AB.Ax.*sin(theta)+AB.Ay.*cos(theta));
elseif strfind(kmodel,'sph')
    AB.kmodel=1;
    Bx=envdx./envu./wt.*sin(phadx-phau);
    By=envdy./envu./wt.*sin(phady-phau);

    Ax=envdx./envu.*cos(phadx-phau)-envdx./envu./envu./wt.*(envudt).*sin(phadx-phau);
    Ay=envdy./envu.*cos(phady-phau)-envdy./envu./envu./wt.*(envudt).*sin(phady-phau);

    AB.Ax=Ax;
    AB.Ay=Ay;
    AB.Bx=Bx;
    AB.By=By;
    
    AB.xc=stc.st(2);
    AB.yc=stc.st(1);
        %-----Azimuth Estimate---------------
    if kazm~=0
        theta=atan2(dudx,dudy);
    else
        theta=atan2(-AB.Bx,-AB.By);
    end
    AB.theta=theta;
        %-----Ray Parameter Estimate-----------
    %AB.pr=-AB.Bx.*sin((theta))-AB.By.*cos((theta));
    AB.pr=sqrt(AB.Bx.*AB.Bx+AB.By.*AB.By);
    %AB.radiation=(AB.Ax.*cos(theta)-AB.Ay.*sin(theta))*stc.dst;
    %AB.gsprd=(AB.Ax.*sin(theta)+AB.Ay.*cos(theta));

    return
    subplot(4,1,1)
    plot(Dc)
    grid on
    
    subplot(4,1,2)
    plot(dudx)
    grid on

    subplot(4,1,3)
    plot(dudy)
    grid on

    subplot(4,1,4)
    plot(AB.Bx)
    grid on
    
end


%
% calculate the instantenous frequency--瞬时频率
%
function wt=InsantFrq(u,dt)
np=length(u);
ut=diff(u)/dt;
ut(np)=ut(np-1);
uh=-imag(hilbert(u));
uth=-imag(hilbert(ut));
uu=abs(hilbert(u));
uu=uu.*uu;
wt=(ut.*uh-u.*uth)./uu;

%
%-----find the peak value of one envelope----------
%
function [kmax]=findmaxima(dat,mmax)
diff0=diff(dat);
ndf=length(diff0);
nmax=0;
for id=2:ndf-1
    if diff0(id)==0 & diff0(id-1)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id;
        %[blw(nmax),bup(nmax)]=errbars(dat,id,0.95);
    elseif diff0(id)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id+1;
        %[blw(nmax),bup(nmax)]=errbars(dat,id+1,0.95);
    end
end 

if nmax==0
    kmax=[];
    return
end

mmax=min(length(pmax),mmax);
for ii=1:mmax
    [fmax,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end
    
