


function WGoupt2mat
fs = '/';
dseis  = './WGmovavrg/';
dataF = dir([dseis '200*.*']);
fout = './wga.R.BHZ.mat';
stinf='./stalist.txt';
periods = [20 40];
stnall=textread(stinf,'%s'); %#ok<DTXTRD>
st = intstrt(periods,stnall,dataF);


for ei =1:length(dataF)
    eloc = dataF(ei).name;
    for  pi= 1:length(periods)
        tprd = num2str2(periods(pi),3,0);
        wgfile  = dir([dseis eloc fs 'wga.R.p' tprd '.BHZ.' eloc '*']);
        if isempty([wgfile.name])
            continue
        end
        
        wgfile  = [dseis eloc fs wgfile.name];
        
        [stn,stla,stlo,~,v0,dv,a0,da,r0,dr,g0,dg,vg,azmo,nsti,~,evla,evlo,evdpth]= ...
         textread(wgfile,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'); %#ok<DTXTRD>
     
        if isempty(stn)
            continue
        end
        
        for si = 1:length(stnall)
            ks = find(ismember(stn,stnall(si)), 1);
            if isempty(ks)
                continue
            end
            
            st(si).stn = stn(ks);
            st(si).stla = stla(ks);
            st(si).stlo = stlo(ks);
            st(si).evla(ei) = evla(ks);
            st(si).evlo(ei) = evlo(ks);
            st(si).dpth(ei) = evdpth(ks);
            
            st(si).periods(ei,pi) = periods(pi);
            st(si).v0(ei,pi) = v0(ks);
            st(si).dv(ei,pi) = dv(ks);
            st(si).a0(ei,pi) = a0(ks);
            st(si).da(ei,pi) = da(ks);
            st(si).r0(ei,pi) = r0(ks);
            st(si).dr(ei,pi) = dr(ks);
            st(si).g0(ei,pi) = g0(ks);
            st(si).dg(ei,pi) = dg(ks);
            st(si).azmo(ei,pi) = azmo(ks);
        end
    end
end
Nani = isnan([st.stla]);
st(Nani) = [];
save(fout,'st','-v7.3')

function st = intstrt(periods,stall,dataF)

pi = length(periods);
ei = length(dataF);
    for si = 1:length(stall)
        st(si).stn = nan;
        st(si).stla = nan;
        st(si).stlo = nan;
        
        st(si).evla(1:ei) = nan;
        st(si).evlo(1:ei) = nan;
        st(si).dpth(1:ei) = nan;
        
        st(si).v0(1:ei,1:pi) = nan;
        st(si).dv(1:ei,1:pi) = nan;
        st(si).a0(1:ei,1:pi) = nan;
        st(si).da(1:ei,1:pi) = nan;
        st(si).r0(1:ei,1:pi) = nan;
        st(si).dr(1:ei,1:pi) = nan;
        st(si).g0(1:ei,1:pi) = nan;
        st(si).dg(1:ei,1:pi) = nan;
        st(si).azmo(1:ei,1:pi) = nan;
        st(si).periods(1:ei,1:pi) = nan;
    end

