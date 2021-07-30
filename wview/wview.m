%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   - by Liang Chuantao, 2009                      %%%%%%%%%%%
%%%   - liangct@cdut.edu.cn                          %%%%%%%%%%%
%%%   - before running the program, also check input %%%%%%%%%%%
%%%   - parameters in function of "wview_OpeningFcn" %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = wview(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wview_OpeningFcn, ...
                   'gui_OutputFcn',  @wview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) 
                               
    gui_State.gui_Callback = str2func(varargin{1});
end                                                 
if nargout  
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% --- Executes just before wview is made visible.
function wview_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
%Parameters seting in this function
hdl=handles;
fcall='''AxesButtonDown_Callback''';
fcdown=['wview(' fcall ',gcbo,[],guidata(gcbo))'];
set(hdl.wview,'WindowButtonDownFcn',fcdown)

keys='\';
datroot='..\SACdata\'; %root path
fddi=dir([datroot '200*']); %name of Event folders
fptn='T*.BHZ'; %name of sacfiles

for i=1:length(fddi)
fdd(i,:)=fddi(i).name;
end

set(handles.DATDIR,'string',fdd(1,:));
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
timplot=[str2double(Tim1) str2double(Tim2)];
par.timplot=timplot;
par.fptn=fptn;

badD=get(handles.EDTDirTo,'string');
par.badD=badD;

Pband=get(handles.Pband,'string');
par.pband=str2num(Pband);

Vgroup=get(handles.VGroup,'string');
par.vgroup=str2double(Vgroup);

nwfPage=get(handles.NWFPage,'string');
par.nf1win=str2double(nwfPage);

datdir1=get(handles.DATDIR,'string');
par.datdir1=datdir1;

par.ksys=keys;
par.fdd=fdd;
par.datroot=datroot;
par.datdir=[par.datroot par.ksys datdir1  par.ksys];
par.outdir=badD;
par.outdir=[par.datdir par.ksys par.outdir par.ksys];
par.yscl0=1.0;
par.dt=0.5;
par.knorm='Constant';

set(handles.DATDIR,'string',fdd(1,:));
set(handles.DATROOT,'string',datroot);
set(hdl.EDTFPTN,'string',par.fptn);
set(handles.EVlist,'string',fdd);

fptn=get(handles.EDTFPTN,'string');
par.fptn0=fptn;

nfls=0;
dir0=[par.datdir];
datf=[dir0 par.fptn0];
sacfiles0=dir(datf);
nn=length(sacfiles0);

fname={sacfiles0.name};
for jj=1:nn
    par.fls{nfls+jj}=[dir0 cell2mat(fname(jj))];
end
nfls=nfls+nn;
par.nfls=nfls;
par.kshw(1:nfls)=0;
par.kchs(1:nfls)=1;
par.istr=1;
par.iend=par.nf1win;
par.kpmode=0;
par.next1=0;
par.wfGUI=1;
par=plotwaves(hdl,par);
userdata.par=par;
set(gcf,'userdata',userdata);
handles.output = hObject;
guidata(hObject, handles);



function varargout = wview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%
%--Button down function--
%
function AxesButtonDown_Callback(hObject, eventdata, handles)
if gca~=handles.AXSWave
    return
end
if gco==handles.AXSWave
    return
end

userdata=get(gcf,'userdata');
par=userdata.par;
nshw=par.nshw;
if strcmp(get(gcf,'SelectionType'), 'extend')
elseif strcmp(get(gcf,'SelectionType'),'normal')
    strmode=get(handles.BTNMode,'string');
    cpt=get(gca,'currentpoint');
    iw=round(cpt(1,2));
    if strcmp(strmode,'view')
        if iw>nshw
            return
        end
        par.iwf=iw;        
        if par.kpmode==0
            ylim([iw-0.5 iw+0.5]);
            par.kpmode=1;
        elseif par.kpmode==1
            nfls=par.nfls;
                      
            if nfls<=10
                disp('Less than 10 stations!!');
                return
            end
            
            nmax=par.nf1win;
            hwin=min(nfls,nmax);
            par.nshw=hwin;
            ylim([-0.5 hwin+0.5]);
            par.kpmode=0;
        end
    else  %time pick #1
        tim=cpt(1,1);
        if par.kpmode==0
            par.tim(iw,1)=tim;
        elseif par.kpmode==1
            iw=par.iwf;
            par.tim(iw,1)=tim;
        end
        delete(findobj('tag',['pick1.' num2str(iw)]))
        plot([tim tim],[iw-0.25 iw+0.25],'r','tag',['pick1.' num2str(iw)])
        vg=par.dis/tim;
        set(handles.VGroup,'string',num2str(vg));
    end
elseif strcmp(get(gcf,'SelectionType'),'alt')
    cpt=get(gca,'currentpoint');
    iw=round(cpt(1,2));
    strmode=get(handles.BTNMode,'string');
    if strcmp(strmode,'view')
        if iw>nshw
            return
        end
        if par.kchs(iw)==0
            par.kchs(iw)=1;
            set(gco,'color','b')
        else
            par.kchs(iw)=0;
            set(gco,'color','y');
        end
    else %time pick 2
        tim=cpt(1,1);
        if par.kpmode==0
            par.tim(iw,2)=tim;
        elseif par.kpmode==1
            iw=par.iwf;
            par.tim(iw,2)=tim;
        end        
        delete(findobj('tag',['pick2.' num2str(iw)]))
        plot([tim tim],[iw-0.25 iw+0.25],'m','tag',['pick2.' num2str(iw)])
        vg=par.dis/tim;
       set(handles.VGroup,'string',num2str(vg));
    end
end


userdata.par=par;
set(gcf,'userdata',userdata);



function EDTTim1_Callback(hObject, eventdata, handles)
userdata=get(gcf,'userdata');
par=userdata.par;
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'userdata',userdata);




% --- Executes during object creation, after setting all properties.
function EDTTim1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDTTim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDTTim2_Callback(hObject, eventdata, handles)
userdata=get(gcf,'userdata');
par=userdata.par;
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes during object creation, after setting all properties.
function EDTTim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDTTim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EDTFPTN_Callback(hObject, eventdata, handles)
% hObject    handle to EDTFPTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDTFPTN as text
%        str2double(get(hObject,'String')) returns contents of EDTFPTN as a double


% --- Executes during object creation, after setting all properties.
function EDTFPTN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDTFPTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BTNRemove.
function BTNRemove_Callback(hObject, eventdata, handles)
% hObject    handle to BTNRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;

i2=0;
fromdir=par.datdir;
for i1=par.istr:par.iend
    i2=i2+1;
    if par.kchs(i2)==0
        fnm0=cell2mat(par.fls(i1));
        fnm1=[fromdir fnm0];
        if exist(fnm0,'file')
            delete(fnm0);
            msg=[num2str(i1) '-' fnm0 ' is deleted!!']
        end
    end
end


% --- Executes on button press in BTNRePlot.
function BTNRePlot_Callback(hObject, eventdata, handles)
% hObject    handle to BTNRePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes on button press in BTNNext.
function BTNNext_Callback(hObject, eventdata, handles)
% hObject    handle to BTNNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
par.istr=par.istr+par.nf1win;
par.iend=par.iend+par.nf1win;
if par.istr>par.nfls
    msg='The End'
    return
end
if par.iend>par.nfls
    par.iend=par.nfls;
end
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes on button press in BTNBack.
function BTNBack_Callback(hObject, eventdata, handles)
% hObject    handle to BTNBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
par.istr=par.istr-par.nf1win;
par.iend=par.iend-par.nf1win;
if par.iend<1
    return
end
if par.istr<1
    par.istr=1;
end
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes on button press in BTNMoveTo.
function BTNMoveTo_Callback(hObject, eventdata, handles)
% hObject    handle to BTNMoveTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
fromdir=par.datdir;
todir0=get(handles.EDTDirTo,'string');
todir=[par.datdir todir0 par.ksys];

if ~exist(todir,'dir')
    mkdir(todir);
end
i2=0;
for i1=par.istr:par.iend
    i2=i2+1;
    if par.kchs(i2)==0
        fnm0=cell2mat(par.fls(i1));
        fnm1=[fromdir fnm0];
        if exist(fnm0,'file')
            movefile(fnm0,todir);
            msg=[num2str(i1) '-' fnm0 'is Moved to ' todir0]
        end
    end
end


function EDTDirTo_Callback(hObject, eventdata, handles)
% hObject    handle to EDTDirTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
fptn=get(handles.EDTFPTN,'string');
par.fptn0=fptn;
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function EDTDirTo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDTDirTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BTNMode.
function BTNMode_Callback(hObject, eventdata, handles)
% hObject    handle to BTNMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(hObject,'string');
if strcmp(str,'view')
    set(hObject,'string','time pick')
else
    set(hObject,'string','view')
end


% --------------------------------------------------------------------
function MENEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MENEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MENCut_Callback(hObject, eventdata, handles)
% hObject    handle to MENCut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
nfls=par.nfls;

if nfls<=10
    disp('Less than 10 stations!!');
    return
end
for ic=par.istr:par.iend
    ic0=ic-par.istr+1;
    sacf=cell2mat(par.fls(ic));
    if exist(sacf,'file')
        tim=par.tim(ic0,1:2);
        kk=find(~isnan(tim), 1);
        if ~isempty(kk)
            [A1,B1,C1,D1]=readsac0(sacf);
            tbeg=A1(6);
            tim=tim-tbeg;
            dt=A1(1);
            nn=round(tim/dt+1);
            if isnan(nn(1)) || nn(1)<=0
                nn(1)=1;
            else 
                A1(6)=A1(6)+tim(1);
            end
            if isnan(nn(2)) || nn(2)>length(D1)
                nn(2)=length(D1);
            else
                A1(7)=A1(6)+tim(2);
            end
            D1=D1(nn(1):nn(2));
            B1(10)=length(D1);
            writesac0(sacf,A1,B1,C1,D1);
            diap(['Cut & Save:' sacf])
        end
    end
    par.tim(ic0,1:2)=[NaN NaN];
end
subplot(handles.AXSWave);
delete(findobj('color','r'));
delete(findobj('color','m'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --------------------------------------------------------------------
function MenSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MENDeleMark_Callback(hObject, eventdata, handles)
% hObject    handle to MENDeleMark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
nfls=par.nfls;

if nfls<=10
    disp('Less than 10 stations!!');
    return
end

for ic=par.istr:par.iend
    ic0=ic-par.istr+1;
    par.tim(ic0,1:2)=[NaN NaN];
end
subplot(handles.AXSWave);
delete(findobj('color','r'));
delete(findobj('color','m'));
userdata.par=par;
set(gcf,'userdata',userdata);


% --------------------------------------------------------------------
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in POPNorm.
function POPNorm_Callback(hObject, eventdata, handles)
% hObject    handle to POPNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns POPNorm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from POPNorm
userdata=get(gcf,'userdata');
contents = get(hObject,'String');
knorm=contents{get(hObject,'Value')};
userdata.par.knorm=knorm;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function POPNorm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to POPNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function MENPlot_Callback(hObject, eventdata, handles)
% hObject    handle to MENPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MENPlotLim_Callback(hObject, eventdata, handles)
% hObject    handle to MENPlotLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.wview,'userdata');
par=userdata.par;

ttl='Input for Ploting Limit:';
prmpt=[{'tlim [50 2000]:'};{'Vertical Sacling Factor (1):'}];
nln=1;

defAns=[{num2str(par.timplot)};{num2str(par.yscl0)}];
ans = inputdlg(prmpt,ttl,nln,defAns,'on');
if isempty(ans)
    return
else
    par.timplot=str2num(cell2mat(ans(1)));
    par.yscl0=str2num(cell2mat(ans(2)));
    userdata.par=par;
    set(handles.wview,'userdata',userdata);
end



function DATDIR_Callback(hObject, eventdata, handles)
% hObject    handle to DATDIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DATDIR as text
%        str2double(get(hObject,'String')) returns contents of DATDIR as a double


% --- Executes during object creation, after setting all properties.
function DATDIR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DATDIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NWFPage_Callback(hObject, eventdata, handles)
% hObject    handle to NWFPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NWFPage as text
%        str2double(get(hObject,'String')) returns contents of NWFPage as a double


% --- Executes during object creation, after setting all properties.
function NWFPage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NWFPage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DATROOT_Callback(hObject, eventdata, handles)
% hObject    handle to DATROOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DATROOT as text
%        str2double(get(hObject,'String')) returns contents of DATROOT as a double


% --- Executes during object creation, after setting all properties.
function DATROOT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DATROOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pband_Callback(hObject, eventdata, handles)
% hObject    handle to Pband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pband as text
%        str2double(get(hObject,'String')) returns contents of Pband as a double


% --- Executes during object creation, after setting all properties.
function Pband_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VGroup_Callback(hObject, eventdata, handles)
% hObject    handle to VGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VGroup as text
%        str2double(get(hObject,'String')) returns contents of VGroup as a double


% --- Executes during object creation, after setting all properties.
function VGroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NEXT.
function NEXT_Callback(hObject, eventdata, handles)
% hObject    handle to NEXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
fdd=par.fdd;
par.next1=par.next1+1;
set(handles.DATDIR,'string',fdd(par.next1,:));
par.istr=1;
nwfPage=get(handles.NWFPage,'string');
par.iend=str2num(nwfPage);
datdir1=get(handles.DATDIR,'string');
par.datdir1=char(datdir1);
par.datdir=[par.datroot par.ksys par.datdir1  par.ksys];
par.outdir=par.badD;
par.outdir=[par.datdir par.ksys par.outdir par.ksys];
nfls=0;
    dir0=[par.datdir];
    datf=[dir0 par.fptn0];
    sacfiles0=dir(datf);
    nn=length(sacfiles0);
    fname={sacfiles0.name};    
    for jj=1:nn
        par.fls{nfls+jj}=[dir0 cell2mat(fname(jj))];
    end
    nfls=nfls+nn;
% %%%%%%%%%%%%%%%%%%
par.nfls=nfls;
par.kshw(1:nfls)=0;
par.kchs(1:nfls)=1;
par.istr=1;
par.iend=par.nf1win;
par.kpmode=0;
% par.next1=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'Userdata',userdata);


% --- Executes on button press in OUTput_Vg.
function OUTput_Vg_Callback(hObject, eventdata, handles)
% hObject    handle to OUTput_Vg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
abv=num2str(par.vgroup);
pb=num2str(par.pband);
FVG=[par.datdir1,' ',abv,' ',pb];
% abcdef=[par.datroot par.ksys 'evg' num2str(par.pband(1)) '-' num2str(par.pband(2)) '.txt'; ];
abcdef=[par.datroot par.ksys 'evg.txt'; ];

fidd=fopen(abcdef,'at');
fprintf(fidd, '%s\n',FVG);
disp(FVG)


% --- Executes on button press in BACK.
function BACK_Callback(hObject, eventdata, handles)
% hObject    handle to BACK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
fdd=par.fdd;
par.next1=par.next1-1;
if par.next1<1
    disp('this is already the first event')
    return
end
set(handles.DATDIR,'string',fdd(par.next1,:));
par.istr=1;
nwfPage=get(handles.NWFPage,'string');
par.iend=str2num(nwfPage);
datdir1=get(handles.DATDIR,'string');
par.datdir1=char(datdir1);
par.datdir=[par.datroot par.ksys par.datdir1  par.ksys];
par.outdir=par.badD;
par.outdir=[par.datdir par.ksys par.outdir par.ksys];
nfls=0;
    dir0=[par.datdir];
    datf=[dir0 par.fptn0];
    sacfiles0=dir(datf);
    nn=length(sacfiles0);
    fname={sacfiles0.name};    
    for jj=1:nn
        par.fls{nfls+jj}=[dir0 cell2mat(fname(jj))];
    end
    nfls=nfls+nn;
% %%%%%%%%%%%%%%%%%%
par.nfls=nfls;
par.kshw(1:nfls)=0;
par.kchs(1:nfls)=1;
par.istr=1;
par.iend=par.nf1win;
par.kpmode=0;
% par.next1=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'Userdata',userdata);


% --- Executes on selection change in EVlist.
function EVlist_Callback(hObject, eventdata, handles)
% hObject    handle to EVlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(gcf,'userdata');
par=userdata.par;
fdd=par.fdd;
next1=get(handles.EVlist,'value');

par.next1=next1;

set(handles.DATDIR,'string',fdd(par.next1,:));
par.istr=1;
nwfPage=get(handles.NWFPage,'string');
par.iend=str2num(nwfPage);
datdir1=get(handles.DATDIR,'string');
par.datdir1=char(datdir1);
par.datdir=[par.datroot par.ksys par.datdir1  par.ksys];
par.outdir=par.badD;
par.outdir=[par.datdir par.ksys par.outdir par.ksys];
nfls=0;
    dir0=[par.datdir];
    datf=[dir0 par.fptn0];
    sacfiles0=dir(datf);
    nn=length(sacfiles0);
    fname={sacfiles0.name};    
    for jj=1:nn
        par.fls{nfls+jj}=[dir0 cell2mat(fname(jj))];
    end
    nfls=nfls+nn;
% %%%%%%%%%%%%%%%%%%
par.nfls=nfls;
par.kshw(1:nfls)=0;
par.kchs(1:nfls)=1;
par.istr=1;
par.iend=par.nf1win;
par.kpmode=0;
% par.next1=0;
%%%%%%%%%%%%%%%%%%%%%%%%%
Tim1=get(handles.EDTTim1,'string');
Tim2=get(handles.EDTTim2,'string');
par.timplot=[str2double(Tim1) str2double(Tim2)];
par=plotwaves(handles,par);
userdata.par=par;
set(gcf,'Userdata',userdata);


% --- Executes during object creation, after setting all properties.
function EVlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EVlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
