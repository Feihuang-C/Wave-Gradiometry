%--------------------------------------------
%---irregular array----
%--------------------------------------------
function [dudx,dudy,stc]=WGA_dxdy(st,ev,fband,par)
vavg=par.vavg;
damp=0.01;
fc=mean(fband);
nst=length(st);
    %---filter the waves---
np=length(st(1).dat);
    %Form the equation system to solve for spatial gradient
    %   d=Gx
stc=st(1);    
azpth=azimuth0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))+180;
if azpth>360
    azpth=azpth-360;
end
for is=2:nst
    ui=st(is).dat;
    u0=st(1).dat;
    di(is-1,1:np)=ui-u0;

    dxi=(st(is).st(2)-st(1).st(2));
    dyi=(st(is).st(1)-st(1).st(1));
    if dxi==0;
    dx(is-1)=0;  
    else
    dx(is-1)=distance0(st(1).st(1),st(1).st(2),st(1).st(1),st(is).st(2))*111.199*dxi/abs(dxi);
    end
    if dyi==0
     dy(is-1)=0;
    else
    dy(is-1)=distance0(st(1).st(1),st(1).st(2),st(is).st(1),st(1).st(2))*111.199*dyi/abs(dyi);
    end
    wi(is-1)=compweight(azpth,st(1),st(is),damp,par.kwi);
    er(is-1)=CompErrors(azpth,st(1),st(is),vavg,fc);
end

%remove station pairs with large truncation errors
kchs(1:nst-1)=ones(1,nst-1);
wi=wi/max(wi);
if ~isempty(par.ermax)
    er=er
    idx=find(er<=par.ermax);
    wi=wi(idx);
    di=di(idx,:);
    dx=dx(idx);
    dy=dy(idx);
end

%apply weighting
nr=length(dx);
if par.kwi>0
    for is=1:nr
        %%%%%%%%%%
        wii=wi(is);
        %%%%%%%%%%%
        %wii=1
        if ~isempty(wii)
            di(is,1:np)=di(is,1:np)*wii;
            dx(is)=dx(is)*wii;
            dy(is)=dy(is)*wii;
        end
    end
end

%form the G - matrix
G=[dx' dy'];
dudxy=pinv(G)*di;
dudx=dudxy(1,:);
dudy=dudxy(2,:);

%%
%% Calculate the weights for station pairs.
%%
function wi=compweight(azpth,st1,st2,damp,kwi)
if ~isempty(azpth)
    azst=azimuth0(st1.st(1),st1.st(2),st2.st(1),st2.st(2));
    dst=distance0(st1.st(1),st1.st(2),st2.st(1),st2.st(2))*111.199;
    daz=(azpth-azst)/180*pi;
    wi=dst*cos(daz);
    wi=0.1/(abs(wi)+damp);
    wi=wi^kwi;
    %wi=wi*abs(sin(daz)^(kwi));
else
    wi=[];
end

%%
%% Calculate the errors for truncation of first order taylor series.
%%
function delta=CompErrors(azpth,st1,st2,vel,fc)
if ~isempty(azpth)
    azst=azimuth0(st1.st(1),st1.st(2),st2.st(1),st2.st(2));
    dst=distance0(st1.st(1),st1.st(2),st2.st(1),st2.st(2)); %in degree
    daz=azpth-azst;
    wi=dst*cos(daz/180.0*pi);
    lamda=vel/fc/111.199;
    delta=pi/lamda*abs(wi);
else
    delta=[];
end
