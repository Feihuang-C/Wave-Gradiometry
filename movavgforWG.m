function movavgforWG
indir = './WG/'; 
outdir = './WGmovavrg/';
evid = dir([indir,'200*']);
wgfileType = 'wga.R.*nz00pct.txt';
fs = '/';
vlim = 0.3;% upper and lower limits of velocity in each periods
smradi=50;%%%smooth radius 
par.knst=5; % minimum number of support stations
isFixedSmrad = 1; %=1, smooth radius is a fixed number of smradi, =0,smooth radius is half wave length
MaxSmrad = deg2km(0.6); %the maxmum smooth radius,only isFixedSmrad=0 is availabel

for kevi = 1:length(evid) %loops for events
   disp(['kevi:  ' num2str(kevi)])
   evi = evid(kevi).name;
   wgid = dir([indir,evi,fs,wgfileType]);
for wi=1:length(wgid);    %½«loops for files
    
    disp(['wi:  ' num2str(wi)])
    cntrprd = strsplit(wgid(wi).name,'.');
    cntrprd = cell2mat(cntrprd(3));
    cntrprd = str2double(cntrprd(2:4));
    
    wginpth = [indir,evi,fs,wgid(wi).name]; 
    wga=wgaread(wginpth);
    
    if isempty(wga.st);
        continue
    end
    
    p.cntrprd = cntrprd;
    p.vavg=median(wga.v);
    kv=p.vavg+[-vlim vlim];
    
	if isFixedSmrad==1
	   p.smrad=smradi;%%%smooth range 
	elseif isFixedSmrad==0
	   p.smrad=min(kv)*0.5*cntrprd;%%%smooth range 
       if p.smrad>MaxSmrad
           p.smrad=MaxSmrad;
       end
    end
    
	prd=p;
    par.kv=kv;
    smrad=prd.smrad;
    wga=datselect(wga,par);
    wga=smoothavg(wga,smrad);
    
    if isempty(wga)
        continue
    end
        
    outdir2 = [outdir evi fs];
    if ~exist(outdir2,'dir');
        mkdir(outdir2)
    end    
    flo=[outdir2  wgid(wi).name num2str2(smrad,4,1) 'km'];
    OutPutVelAzm(flo,wga);    
end
end

%
%---Read the WGA files---
%
function wga=wgaread(fl)

[stla,stlo,stel,v,dv,az,da,rd,dr,gs,dg,vg,ao,ns,amp,stn,~,evla,evlo,evdpth,~]= ...
textread(fl,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f'); %#ok<DTXTRD>

badID = v>5 | v<2;
stn(badID) = [];
stla(badID) = [];
stlo(badID) = [];
stel(badID) = [];
v(badID) = [];
dv(badID) = [];
az(badID) = [];
da(badID) = [];
rd(badID) = [];
dr(badID) = [];
gs(badID) = [];
dg(badID) = [];
vg(badID) = [];
ao(badID) = [];
ns(badID) = [];
amp(badID) = [];
evla(badID) = [];
evlo(badID) = [];
evdpth(badID) = [];


wga.st=[stla stlo stel];
wga.v=v;
wga.dv=abs(dv);
az(az<-180)=az(az<-180)+360;
az(az>180)=az(az>180)-360;
wga.az=az;
wga.da=abs(da);
wga.rd=rd;
wga.dr=abs(dr);
wga.gs=gs;
wga.dg=abs(dg);
wga.vg=vg;
wga.ao=ao;
wga.nst=ns;
wga.amp=amp;
wga.stn=stn;
wga.evla=evla;
wga.evlo=evlo;
wga.evdpth=evdpth;



%%-------------------------------
%% data selection
%%-------------------------------
function wga=datselect(wga,par)
if isfield(par,'kv')
    idx=wga.v>par.kv(1) & wga.v<par.kv(2);
end
if isfield(par,'knst')
    idx=idx.*(wga.nst>par.knst);
end
idx=find(idx);
wga.st=wga.st(idx,1:3);
wga.v=wga.v(idx);
wga.dv(idx)=wga.dv(idx);
wga.az=wga.az(idx);
wga.da=wga.da(idx);
wga.rd=wga.rd(idx);
wga.dr=wga.dr(idx);
wga.gs=wga.gs(idx);
wga.dg=wga.dg(idx);
wga.vg=wga.vg(idx);
wga.ao=wga.ao(idx);
wga.nst=wga.nst(idx);
wga.amp=wga.amp(idx);
wga.stn=wga.stn(idx);
wga.evla=wga.evla(idx);
wga.evlo=wga.evlo(idx);
wga.evdpth=wga.evdpth(idx);



%%
%%average around a station
%%
function bk=smoothavg(wga,dsmax)
if isempty(dsmax)
    bk=wga;
    return
end
st=wga.st;
[nst , ~]=size(st);
bk=[];

for is=1:nst
    la0=st(is,1);
    lo0=st(is,2);
    nn=0;
    idx=[];
    for isi=1:nst
        lai=st(isi,1);
        loi=st(isi,2);    
        ds=distance0(la0,lo0,lai,loi)*111.199;
        if ds<dsmax
            nn=nn+1;
            idx(nn)=isi;
        end
    end
    bk.st(is,1:3)=st(is,1:3);
    bk.v(is)=mean(wga.v(idx));
    bk.dv(is)=mean(wga.dv(idx));
    bk.az(is)=mean(wga.az(idx));
    bk.da(is)=mean(wga.da(idx));
    bk.rd(is)=mean(wga.rd(idx));
    bk.dr(is)=mean(wga.dr(idx));
    bk.gs(is)=mean(wga.gs(idx));
    bk.dg(is)=mean(wga.dg(idx));
    bk.vg(is)=mean(wga.vg(idx));
    bk.ao(is)=mean(wga.ao(idx));
    bk.nst(is)=mean(wga.nst(idx));
    bk.amp(is)=mean(wga.amp(idx));
    bk.stn(is)=wga.stn(is);
    bk.evla(is)=wga.evla(is);
    bk.evdpth(is)=wga.evdpth(is);
    bk.evlo(is)=wga.evlo(is);
end


%-----------------------------------
%------Output results-------------
%-----------------------------------
function OutPutVelAzm(fout,wga)
nst=length(wga.v);
if isempty(nst)
    return
end
fido=fopen(fout,'w');
for is=1:nst
    st=wga.st(is,1:3);
    v0=wga.v(is);
    dv=wga.dv(is);
    a0=wga.az(is);
    da=wga.da(is);
    r0=wga.rd(is);
    dr=wga.dr(is);
    g0=wga.gs(is);
    dg=wga.dg(is);
    vg=wga.vg(is);
    ao=wga.ao(is);
    nst=round(wga.nst(is));
    amp=wga.amp(is);
    stn=cell2mat(wga.stn(is));
    evla=wga.evla(is);
    evlo=wga.evlo(is);
    evdpth=wga.evdpth(is);
    fprintf(fido,'%5s %7.3f %8.3f %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %8.3f %2g %2g %7.3f %7.3f %7.3f\n', ...
    stn,st(1),st(2),st(3),v0,dv,a0,da,r0,dr,g0,dg,vg,ao,nst,amp,evla,evlo,evdpth);
end
fclose(fido);