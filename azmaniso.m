
function azmaniso
global MaxAzm
infile = './wga.R.BHZ.mat';
outpath = './aniso/';
periods = [20 40]; %periods
bootNum = 500; % bootstrap times
vlims = 0.5; % velocity limit (km/s)
azmban = 10;
MaxAzms = 180; 
isbootstrap = 1;
isoverwrite = 0;
WG = load(infile);
WG = WG.st;

if ~exist(outpath,'dir')
    mkdir(outpath)
end

for li = 1:length(vlims) %loops for velocity-limit parameter
    vlim = vlims(li);
    for MxAzmi = 1:length(MaxAzms) %loops for anisotropic periods' parameter
        MaxAzm = MaxAzms(MxAzmi);
        clear Aniso
        outfile = [outpath '2ChuanxiAniso' num2str(MaxAzm) 'ban' num2str(azmban) 'vlim' num2str(vlim) '.mat'];
        
        if exist(outfile,'file')&&isoverwrite==0
            load(outfile)
            kist = length(Aniso.st)-1;
        else
            kist = 1;
        end
        
        for ist = kist:length(WG)
            sti = WG(ist);
            Aniso.stn(ist,:) = cell2mat(sti.stn); Aniso.st(ist,:) = [sti.stla,sti.stlo];
            for iprd = 1:length(periods)
                for iev = 1:length(sti.evla)
                    v(iev) = sti.v0(iev,iprd);
                    dv(iev) = sti.dv(iev,iprd);
                    azmo(iev) = sti.azmo(iev,iprd);
                end
                azmo(azmo<0)=azmo(azmo<0)+360; 
                
                %% remove the nan or zeros values
                dv(isnan(v)) = [];
                azmo(isnan(v)) = [];
                v(isnan(v)) = [];
                dv(v==0) = [];
                azmo(v==0) = [];
                v(v==0) = [];
                stdv2=std(v);
               
                
                if ~isempty(v)
                    nloops=0;
                    cmv=mean(v);
                    dvlim=100;
                    while dvlim>0.001
                        nloops=nloops+1;
                        rmid= abs(v-median(v))>vlims;
                        v(rmid)=[];
                        dv(rmid)=[];
                        azmo(rmid)=[];                        
                        dvlim=abs(cmv-mean(v));
                        cmv=mean(v);
                        if nloops>100
                            v=cmv;
                            break
                        end
                    end
                end
                
                
                %%
                if isempty(v)
                    
                    Aniso.viso(ist,iprd) = nan; Aniso.stdv(ist,iprd) = nan;
                    Aniso.Fai(ist,iprd) = nan; Aniso.M(ist,iprd) = nan;
                    Aniso.a(ist,iprd) = nan; Aniso.b(ist,iprd) = nan;Aniso.fitRsrm(ist,iprd) = nan;
                    
                    if isbootstrap==1
                        Aniso.stdFai(ist,iprd) = nan; Aniso.stdM(ist,iprd) = nan;
                        Aniso.stdA(ist,iprd) = nan; Aniso.stdB(ist,iprd) = nan;
                        if MaxAzm ==360
                            Aniso.stdFai1(ist,iprd) = nan; Aniso.stdM1(ist,iprd) = nan;
                            Aniso.stdA1(ist,iprd) = nan; Aniso.stdB1(ist,iprd) = nan;
                        end
                    end
                    continue
                end
                
                %%
                %%
                if length(v)<9&&~isempty(v)
                    Aniso.viso(ist,iprd) = median(v);Aniso.stdv(ist,iprd) = median(dv);
                    Aniso.Fai(ist,iprd) = nan;Aniso.M(ist,iprd) = nan;
                    Aniso.a(ist,iprd) = nan;Aniso.b(ist,iprd) = nan; Aniso.fitRsrm(ist,iprd) = nan;
                    
                    if MaxAzm ==360
                        Aniso.Fai1(ist,iprd) = nan; Aniso.M1(ist,iprd) = nan;
                        Aniso.a1(ist,iprd) = nan; Aniso.b1(ist,iprd) = nan;
                    end
                    
                    if isbootstrap==1
                        Aniso.stdFai(ist,iprd) = nan;Aniso.stdM(ist,iprd) = nan;
                        Aniso.stdA(ist,iprd) = nan; Aniso.stdB(ist,iprd) = nan;
                        if MaxAzm ==360
                            Aniso.stdFai1(ist,iprd) = nan;Aniso.stdM1(ist,iprd) = nan;
                            Aniso.stdA1(ist,iprd) = nan; Aniso.stdB1(ist,iprd) = nan;
                        end
                    end
                    continue
                end
                %% band averege
                if MaxAzm==180
                    azmo(azmo>180)= azmo(azmo>180)-180;
                end
                bandNum = MaxAzm/azmban;
                k=0;
                for ibn = 1:bandNum
                    bd1 = (ibn-1)*azmban; bd2 = ibn*azmban;
                    vID = find(azmo>=bd1 & azmo<=bd2);
                    if isempty(vID)
                        continue
                    else
                        k = k+1;
                        vt2(k) = median(v(vID));
                    end
                end
                %%
                if ~exist('vt2','var')
                    Aniso.viso(ist,iprd) = median(v);Aniso.stdv(ist,iprd) = median(dv);
                    Aniso.Fai(ist,iprd) = nan; Aniso.M(ist,iprd) = nan;
                    Aniso.a(ist,iprd) = nan; Aniso.b(ist,iprd) = nan;Aniso.fitRsrm(ist,iprd) = nan;
                    
                    if isbootstrap==1
                        Aniso.stdFai(ist,iprd) = nan; Aniso.stdM(ist,iprd) = nan;
                        Aniso.stdA(ist,iprd) = nan; Aniso.stdB(ist,iprd) = nan;
                        if MaxAzm ==360
                            Aniso.stdFai1(ist,iprd) = nan; Aniso.stdM1(ist,iprd) = nan;
                            Aniso.stdA1(ist,iprd) = nan; Aniso.stdB1(ist,iprd) = nan;
                        end
                    end
                    
                    continue
                end
                %%
                if length(vt2)<6&&~isempty(vt2)
                    Aniso.viso(ist,iprd) = median(v);Aniso.stdv(ist,iprd) = median(dv);
                    Aniso.Fai(ist,iprd) = nan;Aniso.M(ist,iprd) = nan;
                    Aniso.a(ist,iprd) = nan;Aniso.b(ist,iprd) = nan; Aniso.fitRsrm(ist,iprd) = nan;
                    
                    if MaxAzm ==360
                        Aniso.Fai1(ist,iprd) = nan; Aniso.M1(ist,iprd) = nan;
                        Aniso.a1(ist,iprd) = nan; Aniso.b1(ist,iprd) = nan;
                    end
                    
                    if isbootstrap==1
                        Aniso.stdFai(ist,iprd) = nan;Aniso.stdM(ist,iprd) = nan;
                        Aniso.stdA(ist,iprd) = nan; Aniso.stdB(ist,iprd) = nan;
                        if MaxAzm ==360
                            Aniso.stdFai1(ist,iprd) = nan;Aniso.stdM1(ist,iprd) = nan;
                            Aniso.stdA1(ist,iprd) = nan; Aniso.stdB1(ist,iprd) = nan;
                        end
                    end
                    continue
                end
                %%
                fitNum = length(v);
                if MaxAzm==180
                    azmo(azmo>180)= azmo(azmo>180)-180;
                end
                
                
                
                An = CompAniso(azmban,v,azmo,MaxAzm,vlims);
                
                
                Aniso.viso(ist,iprd) = An.viso;Aniso.Fai(ist,iprd) = An.Fai;
                Aniso.M(ist,iprd) = An.M;Aniso.fitRsrm(ist,iprd) = An.fitRsrm;
                Aniso.a(ist,iprd) = An.a;Aniso.b(ist,iprd) = An.b;
                Aniso.stdv(ist,iprd) = median(dv);
                Aniso.stdv2(ist,iprd) = stdv2;

                
                
                if isbootstrap==1
                    
                    boot = bootaniso(fitNum,v,azmo,azmban,bootNum,MaxAzm,vlims);
                    
                    Aniso.stdFai(ist,iprd) = boot.stdFai;Aniso.stdA(ist,iprd) = boot.stdA;
                    Aniso.stdM(ist,iprd) = boot.stdM;Aniso.stdB(ist,iprd) = boot.stdB;
                    if MaxAzm == 360
                        Aniso.stdM1(ist,iprd) = boot.stdM1;Aniso.stdFai1(ist,iprd) = boot.stdFai1;
                        Aniso.stdB1(ist,iprd) = boot.stdB1;Aniso.stdA1(ist,iprd) = boot.stdA1;
                    end
                end
                if MaxAzm == 360
                    Aniso.Fai1(ist,iprd) = An.Fai1;Aniso.M1(ist,iprd) = An.M1;
                    Aniso.a1(ist,iprd) = An.a1;Aniso.b1(ist,iprd) = An.b1;
                end
                Aniso.period(iprd)=periods(iprd);
            end
            save(outfile,'Aniso')
        end
        save(outfile,'Aniso')
    end    
end
Anisomat2gmt(outfile,outpath,isbootstrap)




function boot = bootaniso(fitNum,v,azmo,azmban,bootNum,MaxAzm,vlims)
global bootviso
%bootstrap method
bootID = 1:fitNum;
bootID = bootstrp(bootNum, @bootline,bootID);
for iboot = 1:bootNum
    bootV = v(bootID(iboot,:));
    bootazm = azmo(bootID(iboot,:));
    bootdat = [bootV',bootazm'];
    bootdat=sortrows(bootdat,2);
    bootV = bootdat(:,1);
    bootazm = bootdat(:,2);
    
    %          band averege
    bandNum = MaxAzm/azmban;
    k=0;
    for ibn = 1:bandNum
        bd1 = (ibn-1)*azmban;
        bd2 = ibn*azmban;
        vID = find(bootazm>=bd1 & bootazm<=bd2);
        if isempty(vID)
            continue
        else
            k = k+1;
            bootv2(k) = median(bootV(vID));
            vi=bootV(vID);
            
            
            vi(abs(vi-bootv2(k))>vlims)=[];
            
            nloops=0;
            dvlim=100;
            while dvlim>0.001
                nloops=nloops+1;
                rmid= abs(bootv2(k)-vi)>vlims;
                vi(rmid)=[];
                dvlim=abs(bootv2(k)-mean(vi));
                bootv2(k)=mean(vi);
                if nloops>100
                    break
                end
            end
            
            
            bootazmi = bootazm(vID);
            ddv = abs(bootV(vID)-bootv2(k));
            azmID = ddv==min(ddv);
            if length(azmID)>1
                bootazm2 (k) = mean(bootazmi(azmID));
            elseif length(azmID)==1
                bootazm2 (k) = bootazmi(azmID);
            end
        end
    end
    
    if ~exist('bootv2','var')
        anisoPr(iboot,:) = [nan,nan];
        bootM(iboot) = nan;
        bootFai(iboot)=nan;
        continue
    end
    bootviso = mean(bootv2);
    bootazm2 = bootazm2/180*pi;
    
    if MaxAzm ==360
        co = [0,0,0,0];
    elseif MaxAzm ==180;
        co = [0,0];
    end
    [anisoPr(iboot,:)] = lsqcurvefit(@ansioboot,co,bootazm2,bootv2);
    bootM(iboot) = sqrt(anisoPr(iboot,1).^2+anisoPr(iboot,2).^2)/bootviso*100*2;
    bootFai(iboot) = atan2(anisoPr(iboot,2),anisoPr(iboot,1))*90/pi;
    
    if MaxAzm ==360
        bootM1(iboot) = sqrt(anisoPr(iboot,3).^2+anisoPr(iboot,4).^2)/bootviso*100*2;
        bootFai1(iboot) = atan2(anisoPr(iboot,4),anisoPr(iboot,3))*180/pi;
    end
end
boot.mdbootA = median(anisoPr(:,1));
boot.mdbootB = median(anisoPr(:,2));

boot.mdbootM = median(bootM);
boot.mdbootFai = median(bootFai);
boot.stdM = std(bootM);
boot.stdFai = std(bootFai);
boot.stdA = std(anisoPr(:,1));
boot.stdB = std(anisoPr(:,2));

if MaxAzm ==360
    boot.mdbootA1 = median(anisoPr(:,3));
    boot.mdbootB1 = median(anisoPr(:,4));
    boot.mdbootM1 = median(bootM1);
    boot.stdM1 = std(bootM1);
    boot.mdbootFai1 = median(bootFai1);
    boot.stdFai1 = std(bootFai1);
    boot.stdA1 = std(anisoPr(:,3));
    boot.stdB1 = std(anisoPr(:,4));
end

function Aniso = CompAniso(azmban,v,azm,MaxAzm,vlims)
global viso

%band averege
bandNum = MaxAzm/azmban;
k=0;
for ibn = 1:bandNum
    bd1 = (ibn-1)*azmban;
    bd2 = ibn*azmban;
    vID = find(azm>=bd1 & azm<=bd2);
    if isempty(vID)
        continue
    else
        k = k+1;
        vi=v(vID);
        v2(k) = median(vi);
        vi(abs(vi-v2(k))>vlims)=[];
        nloops=0;
        dvlim=100;     
        while dvlim>0.001
            nloops=nloops+1;
            rmid= abs(v2(k)-vi)>vlims;
            vi(rmid)=[];
            dvlim=abs(v2(k)-median(vi));
            v2(k)=median(vi);
            if nloops>100
                break
            end
        end
        
        azmi = azm(vID);
        ddv = abs(v(vID)-v2(k));
        azmID = ddv==min(ddv);
               
        if length(azmID)>1
            azm2 (k) = mean(azmi(azmID));
        elseif length(azmID)==1
            azm2 (k) = azmi(azmID);
        end
    end
end


viso = mean(v2);
if MaxAzm ==360
    co = [0,0,0,0];
elseif MaxAzm ==180;
    co = [0,0];
end
azm2 = azm2/180*pi;
[anisoPr,anisoRsrm] = lsqcurvefit(@ansiofit,co,azm2,v2);
Aniso.viso = viso;
Aniso.Fai = atan2(anisoPr(2),anisoPr(1))*90/pi;
Aniso.M = sqrt(anisoPr(1).^2+anisoPr(2).^2)/viso*100*2;
Aniso.a = anisoPr(1);
Aniso.b = anisoPr(2);
Aniso.fitRsrm = anisoRsrm;
if MaxAzm ==360
    Aniso.Fai1 = atan2(anisoPr(4),anisoPr(3))*180/pi;
    Aniso.M1 = sqrt(anisoPr(3)^2+anisoPr(4)^2)/viso*100*2;
    Aniso.a1 = anisoPr(3);
    Aniso.b1 = anisoPr(4);
end



function y1=bootline(y)
y1=y;

function fy1=ansioboot(abc,fx)
global bootviso  MaxAzm
if MaxAzm==360
    fy1=bootviso+abc(1)*cos(2*fx)+abc(2)*sin(2*fx)+abc(3)*cos(fx)+abc(4)*sin(fx);
end
if MaxAzm==180
    fy1=bootviso+abc(1)*cos(2*fx)+abc(2)*sin(2*fx);
end


function fy1=ansiofit(abc,fx)
global viso MaxAzm
if MaxAzm==360
    fy1=viso+abc(1)*cos(2*fx)+abc(2)*sin(2*fx)+abc(3)*cos(fx)+abc(4)*sin(fx);
end
if MaxAzm==180
    fy1=viso+abc(1)*cos(2*fx)+abc(2)*sin(2*fx);
end




