function par=plotwaves(hd,par)
if par.wfGUI==1
    Pband=get(hd.Pband,'string');
    Vgroup=get(hd.VGroup,'string');
    nwfPage=get(hd.NWFPage,'string');
    par.pband=str2num(Pband);
    par.vgroup=str2double(Vgroup);
    par.nf1win=str2double(nwfPage);
    par=plotwavesGUI(hd,par);
else
    Pband=get(hd.Pband,'string');
    Vgroup=get(hd.VGroup,'string');
    nwfPage=get(hd.NWFPage,'string');
    par.pband=str2double(Pband);
    par.vgroup=str2double(Vgroup);
    par.nf1win=str2double(nwfPage);
    par=plotwaves02(hd,par);
end

%
%
function par=plotwavesGUI(hd,par)
tplot=par.timplot;  %---decide the time range to be plotted
nfls=par.nfls;

if nfls<=10
    disp('Less than 10 stations!!');
    return
end

nmax=par.nf1win;
hwin=min(nfls,nmax);
par.nshw=hwin;

subplot(hd.AXSWave);
hold off;
disp([par.istr par.iend])
yscl0=par.yscl0;
for ic=par.istr:par.iend
    ic0=ic-par.istr+1;
    par.kshw(ic)=1;
    if ic>length(par.fls)
          sacf0=par.fls;
    else
        sacf0=cell2mat(par.fls(ic));
    end
     sacf=sacf0;  
        
    
    if exist(sacf,'file')
    par.kchs(ic0)=1;
    par.tim(ic0,1:2)=[NaN NaN];
    [A1,B1,C1,D1]=readsac0(sacf);
    tttt=find(max(abs(D1)));
    dt=A1(1);
    tpk=A1(11:20);
    tpk(isnan(tpk))=[];
    pband=par.pband;
    if ~isempty(par.pband)
        D1=bpfilt(D1',dt,1/pband(2),1/pband(1));
    end
    %env=abs(hilbert(D1));
    dis=A1(51);
    par.dis=dis;
    azm=A1(52);
    tbeg=[A1(6)];
    tend=[A1(7)];
    vgmax=par.vgroup;
    if strcmp(par.knorm,'Constant')        
        if ic==par.istr
            scl=max(abs(D1));
        end
    elseif strcmp(par.knorm,'Maximum')
        scl=max(abs(D1));
    elseif strcmp(par.knorm,'RMS')
        scl=sqrt(sum(D1.*D1)/length(D1));
    end
    scl=scl/yscl0;
    par.dt=dt;
    D1=D1/scl/2;
    %env=env/scl/2;

%     eq=[A1(36) A1(37) A1(38)];
%     st=[A1(32) A1(33) A1(34)];
%     ipnt=round(tplot(1)/dt):round(tplot(2)/dt);
    np=length(D1);
    timp=(1:np)*dt-dt+tbeg;
    timp(tttt);
    plot(timp,ic0+D1,'tag','waves','color','b');
    hold on;
    %plot(timp,ic0+env,'tag','waves','color','g');
    
    ts=dis/par.vgroup;
    plot([ts ts],ic0+[-0.5 0.5],'r');
    
    %%plot orignal times
    npk=length(tpk);
    for ii=1:npk
        plot([tpk(ii) tpk(ii)],ic0+[-0.5 0.5],'k');
    end
    
    %txt=[cell2mat(C1(2)) cell2mat(C1(1))];
    kn=strfind(sacf,'\');
    txt=sacf(kn(length(kn))+1:length(sacf));
    text(tplot(1)+1.0,ic0-0.4,txt, 'color','r','fontsize',8);
    end
end

xlim([tplot(1) tplot(2)]);
ylim([-0.5 hwin+0.5]);
%set(gca,'YDir','reverse');

par.taxs.min=timp(1);
par.taxs.max=timp(length(timp));
par.taxs.dt=dt;
par.nchs=par.nfls;

%
%--------------Plot in seperate window------------------
%
function par=plotwaves02(~,par)
tplot=par.timplot;  %---decide the time range to be plotted
vgrp=par.vgroup;
nfls=par.nfls;

if nfls<=10
    disp('Less than 10 stations!!');
    return
end

nmax=par.nflwin;
if strcmp(par.kdata,'ccn')
    hwin=max(nfls,nmax);
    par.nshw=nfls;
    sum=0;
else
    hwin=min(nfls,nmax);
    par.nshw=hwin;
end

figure(10)
hwf=subplot('position', [0.1 0.25 0.8 0.7]);
hstk=subplot('position',[0.1 0.08 0.8 0.1]);

subplot(hwf);
hold off;
for ic=par.istr:par.iend
    ic0=ic-par.istr+1;
    par.kshw(ic)=1;
    sacf0=cell2mat(par.fls(ic));
    %sacf=[par.datdir sacf0];
    sacf=sacf0;
    if exist(sacf,'file')
    par.kchs(ic0)=1;
    [A1,B1,C1,D1]=readsac0(sacf);
    dt=A1(1);
    thlf=(length(D1)-1)*dt/2;
    ibeg=round((tplot(1)+thlf)/dt)+1;
    iend=round((tplot(2)+thlf)/dt)+1;
    if strcmp(par.kdata,'ccn')
        par.wv(ic0,:)=D1;
    end
    D1=D1(ibeg:iend);
    if par.kflt(1)==1
        D1=bpfilt(D1,dt,par.bpf(1),par.bpf(2));
    end
    dis=A1(51);
    par.dis=dis;
%     azm=A1(52);
    if strcmp(par.kdata,'ccn')
        sum=sum+D1;
    end
    D1=D1/max(abs(D1))/2;
    t0=dis/vgrp;
    par.tsw=t0;

%     eq=[A1(36) A1(37) A1(38)];
%     st=[A1(32) A1(33) A1(34)];
%     ipnt=round(tplot(1)/dt):round(tplot(2)/dt);
    %nn=min(length(ipnt),length(D1));
    nh2=length(D1);
    nh=round(nh2/2-0.5);
    timp=(-nh:1:nh)*dt;
    plot(timp,ic0+D1(1:length(timp)),'color',[0.2 0.2 0.2],'tag','waves');
    hold on;
    clr=[0.7 0.7 0.7];
    plot([t0 t0],[ic0-0.5 ic0+0.5],'color',clr,'tag','TMark1');
    plot([-t0 -t0],[ic0-0.5 ic0+0.5],'color',clr,'tag','TMark1');
    par.tp=interp1(par.rtpa.dst,par.rtpa.tP,dis/111.199);
    %par.tp=A1(12);
    plot([par.tp par.tp],[ic0-0.5 ic0+0.5],'color',clr);
    plot([-par.tp -par.tp],[ic0-0.5 ic0+0.5],'color',clr);
    par.ts=interp1(par.rtpa.dst,par.rtpa.tS,dis/111.199);
    %par.ts=A1(13);
    plot([par.ts par.ts],[ic0-0.5 ic0+0.5],'color',clr);
    plot([-par.ts -par.ts],[ic0-0.5 ic0+0.5],'color',clr);
    dd0=round(dis/111.199);
    stn1=cell2mat(C1(2));
    stn2=cell2mat(C1(1));
    if strcmp(par.kdata,'ccn')
        k1=strfind(sacf0,stn1);
        k2=strfind(sacf0,stn2);
        txt=sacf0(k1:k2+length(stn2)+2+8);
    else
        txt=[cell2mat(C1(2)) cell2mat(C1(1)) ' DIS=' num2str2(dd0,2) sacf0(length(sacf0)-8:length(sacf0))];
    end
    text(tplot(1)+10.0,ic0-0.4,txt, 'color','k','fontsize',8);
    end
end
xlim([tplot(1) tplot(2)]);
ylim([-0.5 hwin+0.5]);
set(gca,'YDir','reverse')
if ~strcmp(par.kdata,'ccn')
    return
end
par.sachd.A=A1;
par.sachd.B=B1;
par.sachd.C=C1;

D1=sum/nfls;

subplot(hstk)
hold off;
hold off;
plot(timp,D1(1:length(timp)),'color',[0.2 0.2 0.2],'HitTest','off','selectionhighlight','off');
hold on;
plot([par.tsw par.tsw],[min(D1) max(D1)],'color',clr,'tag','TMark1');
plot([-par.tsw -par.tsw],[min(D1) max(D1)],'color',clr,'tag','TMark1');
plot([par.ts par.ts],[min(D1) max(D1)],'color',clr,'tag','TMark1');
plot([-par.ts -par.ts],[min(D1) max(D1)],'color',clr,'tag','TMark1');
plot([par.tp par.tp],[min(D1) max(D1)],'color',clr,'tag','TMark1');
plot([-par.tp -par.tp],[min(D1) max(D1)],'color',clr,'tag','TMark1');

xlim([tplot(1) tplot(2)]);
ylim([min(D1) max(D1)]);
par.taxs.min=timp(1);
par.taxs.max=timp(length(timp));
par.taxs.dt=dt;
par.nchs=par.nfls;
par.stk=D1;
xlabel('Time (sec)')
