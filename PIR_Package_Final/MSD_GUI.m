function varargout = MSD_GUI(varargin)
% MSD_GUI MATLAB code for MSD_GUI.fig
%      MSD_GUI, by itself, creates a new MSD_GUI or raises the existing
%      singleton*.
%
%      H = MSD_GUI returns the handle to a new MSD_GUI or the handle to
%      the existing singleton*.
%
%      MSD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSD_GUI.M with the given input arguments.
%
%      MSD_GUI('Property','Value',...) creates a new MSD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MSD_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MSD_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MSD_GUI

% Last Modified by GUIDE v2.5 22-Jul-2019 11:36:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
FILENAME="";
PATHNAME="";
handles.data=[];
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MSD_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MSD_GUI_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before MSD_GUI is made visible.
function MSD_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MSD_GUI (see VARARGIN)

% Choose default command line output for MSD_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MSD_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MSD_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

% 
%     [imname,handles.imdirec] = uigetfile('*.*','Select an image file',handles.data.imdirec);
% 
%         i=strfind(imname,'.');
%         handles.data.imext=imname((i(end)+1):end);
%         % This finds everything that is not a number
%         j = regexp(imname(1:i(end)-1),'\D');
%         % Now looking at the first thing that is not a number excluding the
%         % extension to determine the number of zeros.
%         zeros=i(end)-1-j(end);
%         handles.data.imbase=imname(1:(i(end)-1-zeros));
%         handles.data.imzeros=num2str(zeros);
%         fstart=str2double(imname((i(end)-zeros):(i(end)-1)));
%         handles.data.imfstart=num2str(fstart);
%         dirinfo = dir([handles.data.imdirec handles.data.imbase '*.' handles.data.imext]);
%         handles.data.imfend=num2str(str2double(dirinfo(end).name(i(end)-zeros:i(end)-1))-str2double(C));
%     end
%     set(handles.imagedirectory,'string',handles.data.imdirec);
%     set(handles.imagebasename,'string',handles.data.imbase);
%     set(handles.imagezeros,'string',handles.data.imzeros);
%     set(handles.imageextension,'string',handles.data.imext);
%     set(handles.imageframestart,'string',handles.data.imfstart);
%     set(handles.imageframeend,'string',handles.data.imfend);
% %     if strcmp(handles.data.masktype,'dynamic')
% %         load_masklist(handles)
% %     end
%     %load_imlist(handles);
%     guidata(hObject,handles)
    
    



[FILENAME, handles.data.imdirec]  = uigetfile('*.tif','select an image','E:\PID\phantom\');

i=strfind(FILENAME,'.');
handles.data.imext=FILENAME((i(end)+1):end);
% This finds everything that is not a number
j = regexp(FILENAME(1:i(end)-1),'\D')
% Now looking at the first thing that is not a number excluding the
% extension to determine the number of zeros.
zeros=i(end)-1-j(end)
handles.data.imbase=FILENAME(1:(i(end)-1-zeros));
handles.data.imzeros=num2str(zeros);
fstart=str2double(FILENAME((i(end)-zeros):(i(end)-1)));
handles.data.imfstart=num2str(fstart);
dirinfo = dir([handles.data.imdirec handles.data.imbase '*.' handles.data.imext]);

guidata(hObject,handles)





















function Number_Callback(hObject, eventdata, handles)
% hObject    handle to lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lb as text
%        str2double(get(hObject,'String')) returns contents of lb as a double


% --- Executes during object creation, after setting all properties.
function lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MAXCOR_Callback(hObject, eventdata, handles)
% hObject    handle to MAXCOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAXCOR as text
%        str2double(get(hObject,'String')) returns contents of MAXCOR as a double


% --- Executes during object creation, after setting all properties.
function MAXCOR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAXCOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MINCOR_Callback(hObject, eventdata, handles)
% hObject    handle to MINCOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MINCOR as text
%        str2double(get(hObject,'String')) returns contents of MINCOR as a double


% --- Executes during object creation, after setting all properties.
function MINCOR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MINCOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% returns toggle state of checkbox1


% --- Executes on button press in gausb.
function gausb_Callback(hObject, eventdata, handles)
% hObject    handle to gausb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gausb


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N=str2num(get(handles.Number,'string'));
imdir =handles.data.imdirec;%strcat('C:\Users\Adib\Desktop\Projects\Montecarlo_diffiusion\128_onlydiff2\imagefilesMP1_Dt2C0.001_runnumber_1\');
filename=handles.data.imbase;%strcat('image_Motility1concentration_0.001_Dp5_Diff2_')
numdig= strcat('%0',handles.data.imzeros,'i');
filetype=strcat('.',handles.data.imext);
hCheckboxes = [handles.checkbox1 ; handles.gausb;handles.checkbox3;handles.checkbox4;handles.checkbox5;handles.checkbox_strdle]
checkboxValues = get(hCheckboxes, 'Value');
checkboxVal=cell2mat(checkboxValues);
tmax=str2num(get(handles.MAXCOR,'string'));
tmin=str2num(get(handles.MINCOR,'string')); %are min and max time lags
windsize=str2num(get(handles.winsz,'string'));
imagestart=str2num(get(handles.imstart,'string'));
dofilter=checkboxVal(1);    %gaussian Filter on the whole image
gausblur=checkboxVal(2); % Gaussian blur for unresolved particles
mbf=checkboxVal(3); %mean before Gaussian filter
maf=checkboxVal(4); %mean after Gaussianfilter
ensemf=checkboxVal(5);%ensemble in Fourier
stradlemod=checkboxVal(6);
savedir=strcat('dofilter',num2str(dofilter),'_gaussb',num2str(gausblur),'_mbf',num2str(mbf),'maf',num2str(maf),'ensemf',num2str(ensemf))
resave=strcat(imdir,savedir);
mkdir(resave);
handles.data.resave=resave
try handles.data.bkg
    bkga=handles.data.bkg;
catch
    bkga=0;
end


try 
    [diffusionestimate]=MSD_cor(imdir,filename,numdig,filetype,resave,dofilter,N,tmin,tmax,gausblur,mbf,maf,ensemf,windsize,stradlemod,imagestart,bkga);
    saveset = strcat(resave,'\settings');
    save(saveset,'checkboxVal')
    saveset = strcat(resave,'\difestimate_',num2str(N));
    save(saveset,'diffusionestimate')
    set(handles.motility_value,'string',num2str(diffusionestimate.pdfofensemble_Adib(1)))
    conver=str2num(get(handles.conv,'string'));
    frameps=str2num(get(handles.fps,'string'));
    correct_dif=diffusionestimate.pdfofensemble_Adib(1)/frameps*conver^2;
    set(handles.motilityunit,'string',num2str(correct_dif))
    
end
guidata(hObject,handles)
F = findall(0,'type','figure','tag','TMWWaitbar'); 
delete(F);

% --- Executes during object creation, after setting all properties.
function Number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in MEANISUB.
function MEANISUB_Callback(hObject, eventdata, handles)
% hObject    handle to MEANISUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N=str2num(get(handles.Number,'string'));
imdir =handles.data.imdirec;%strcat('C:\Users\Adib\Desktop\Projects\Montecarlo_diffiusion\128_onlydiff2\imagefilesMP1_Dt2C0.001_runnumber_1\');
filename=handles.data.imbase;%strcat('image_Motility1concentration_0.001_Dp5_Diff2_')
numdig= strcat('%0',handles.data.imzeros,'i');
filetype=strcat('.',handles.data.imext);





mCheckboxes = [handles.checkbox6 ; handles.checkbox7;handles.checkbox8;handles.checkbox10]
mcheckboxValues = get(mCheckboxes, 'Value');
mcheckboxVal=cell2mat(mcheckboxValues);

MIN=mcheckboxVal(1)    %MIN Image sub
mean=mcheckboxVal(2)   %Mean Image sub
pcolor=mcheckboxVal(3) %Signal Color
normalalize=mcheckboxVal(4) %Normalize

processed_dir=strcat(imdir,'MIN',num2str(MIN),'mean',num2str(mean),'Signalcolor',num2str(pcolor),'Normalize',num2str(normalalize),'\')
mkdir(processed_dir)
%mean=0;%%if you want mean subtraction set mean to 1
%pcolor=1;%%what color is your signal ? it the particles are white set pcolor=1 else
%%if particle is black pcolor=0
outlier =0; % if in your images there is black or white broken pixel which you
%%want to change to background color
%normalalize=0; % set normalize=1 if you want to normalize the image afterward.
%%It will brighten the signal
dosave=0;

bkga=Meansub(imdir,filename,numdig,filetype,processed_dir,mean,pcolor,outlier,N,normalalize,dosave);
imagesc(bkga)
handles.data.bkg=bkga;
guidata(hObject,handles);



% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in msdgen.
function msdgen_Callback(hObject, eventdata, handles)
% hObject    handle to msdgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

adiboutlier=1;
matout=1;
resave=handles.data.pdfdirec;
filename=handles.data.pdfbase;%strcat('image_Motility1concentration_0.001_Dp5_Diff2_')
%numdig= strcat('%0',handles.data.pdfzeros,'i')
%filetype=strcat('.',handles.data.pdfext)

tmax=str2num(get(handles.MAXCOR,'string'));
tmin=str2num(get(handles.MINCOR,'string'));

[ms2]=PDF2MSD(adiboutlier,matout,resave,filename,tmin,tmax)

handles.pltw=plot(ms2,'r*');
grid on
xlabel('time lag');
ylabel('MSD');
figure(2)
plot(ms2,'r*');
grid on
xlabel('time lag');
ylabel('MSD');
% ylabel(handles.pltw,'ylabel')
% handles.pltw.xlabel='xlabel'







% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FILENAME, handles.data.pdfdirec]  = uigetfile('*.mat','select a PDF file',handles.data.imdirec);

i=strfind(FILENAME,'.');
handles.data.pdfext=FILENAME((i(end)+1):end);
% This finds everything that is not a number
j = regexp(FILENAME(1:i(end)-1),'\D');
% Now looking at the first thing that is not a number excluding the
% extension to determine the number of zeros.
zeros=i(end)-1-j(end);
handles.data.pdfbase=FILENAME(1:(i(end)-1-zeros));
handles.data.pdfzeros=num2str(zeros);
fstart=str2double(FILENAME((i(end)-zeros):(i(end)-1)));
handles.data.pdfstart=num2str(fstart);
dirinfo = dir([handles.data.pdfdirec handles.data.pdfbase '*.' handles.data.pdfext]);

guidata(hObject,handles);





function conv_Callback(hObject, eventdata, handles)
% hObject    handle to conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conv as text
%        str2double(get(hObject,'String')) returns contents of conv as a double


% --- Executes during object creation, after setting all properties.
function conv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double


% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winsz_Callback(hObject, eventdata, handles)
% hObject    handle to winsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winsz as text
%        str2double(get(hObject,'String')) returns contents of winsz as a double


% --- Executes during object creation, after setting all properties.
function winsz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_strdle.
function checkbox_strdle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strdle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strdle



function imstart_Callback(hObject, eventdata, handles)
% hObject    handle to imstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imstart as text
%        str2double(get(hObject,'String')) returns contents of imstart as a double


% --- Executes during object creation, after setting all properties.
function imstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
