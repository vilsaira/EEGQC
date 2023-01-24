%{
EEGQC v5: automated tool for EEG quality analysis.

Copyright (C) 2013, 2015, 2016, Viljami Sairanen (viljami.sairanen@gmail.com)

This software is published under the BSD 2-Clause License
and is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%}

%{
EEG QC v. 5.0
Erroreita karsittu. Tuki Matlab 2014b->. -Teemu Mäkelä
17.7. Luotu oma nansum-funktio, jotta statistic toolboxia ei käytettäisi
31.7. memmapfile-osio korjattu

EEG QC v. 4.0
rivit 1203-1206 kommentoitu (ei poista excelin ylimääräisiä välilehtiä.

EEG QC v. 3.0 on HUS-Kuvantamisen KNF-vastuualueen fyysikoiden
kehittämän laadunvalvontaprotokollan mukaisten mittausten analysoimiseen
kirjoitettu ohjelma. Ohjelman tarkoitus on automatisoida ja siten nopeuttaa
ja tarkentaa protokollan mukaista tulosanalyysia.

Ohjelman käyttäminen ja tulosten analysointi on esitetty Viljami Sairasen
labratyössä "EEG laadunvalvonta" .


---------------------------------------------------------------------------
Kirjaa tähän tekemäsi muutokset, niiden päivämäärä ja mistä löytyy vanha
versio, jos sitä on siirretty.

-- 30.9.2013 -- Viljami Sairanen
Ohjelma on kirjoitettu "EEG Laadunvalvonta"-menetelmäohjeen versiolle 1.0,
jonka voimaantulopäivä oli 6.2.2013. 

---------------------------------------------------------------------------
%% Out-of-Memory ongelmat: kirjoita Matlabin komentoriville 
    >> memory
    Maximum possible array: XXX MB <--- tämä kiinnostaa
Tämän jälkeen tarkista analyysoitavan ASCII-tiedoston koko. Jos koko on
lähes XXX tai suurempi niin datan lukeminen ohjelaan ei vaan onnistu. Yritä
vapauttaa muistia sulkemalla syöppöjä kuten Firefox. Jos tämä ei vielä
riitä, boottaa 32 bittinen tietokoneesi /3GB switchillä (varaudu
komplikaatioihin). Tai vaihtoehtoisesti siirrä datat ja ohjelma
merkuriukselle, 64 bittisellä käyttiksellä kaikki onnistuu. Paitsi
automaattisen raportin generoiminen...

Voit toki ratkaista ongelman lisäämällä ohjelmaan kyvyn lukea ASCII-tiedostoa
lennosta analyysin aikana, jolloin dataa ei tarvitse lukea kerralla
ohjelmaan. Hidastaa, mutta toimii varmasti. (Tämä oli alkuperäinen syy,
miksi dataformaatiksi valittiin ASCII eikä idioottimaisesti rakennettu
EDF, mutta koska rammia riitti ei jaksettu koodata loppuun).
%}

% Initialization code - DO NOT EDIT
function varargout = EEGQC_v5(varargin)
% EEGQC_V5 MATLAB code for EEGQC_v5.fig
%      EEGQC_V5, by itself, creates a new EEGQC_V5 or raises the existing
%      singleton*.
%
%      H = EEGQC_V5 returns the handle to a new EEGQC_V5 or the handle to
%      the existing singleton*.
%
%      EEGQC_V5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGQC_V5.M with the given input arguments.
%
%      EEGQC_V5('Property','Value',...) creates a new EEGQC_V5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EEGQC_v5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EEGQC_v5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EEGQC_v5

% Last Modified by GUIDE v2.5 08-Jun-2015 10:36:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EEGQC_v5_OpeningFcn, ...
                   'gui_OutputFcn',  @EEGQC_v5_OutputFcn, ...
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


% --- Executes just before EEGQC_v5 is made visible.
function EEGQC_v5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EEGQC_v5 (see VARARGIN)
fprintf('Starting...\n');
handles.frequency = 1;
handles.nChannels = 32;
handles.currentChannel = 1;

% set pushbutton images
if exist('bluearrow.PNG', 'file');
    pbIml2r = imread('bluearrow.PNG');
    r = fliplr(pbIml2r(:,:,1));
    g = fliplr(pbIml2r(:,:,2));
    b = fliplr(pbIml2r(:,:,3));
    pbImr2l(:,:,1) = r;
    pbImr2l(:,:,2) = g;
    pbImr2l(:,:,3) = b;
    set(handles.pushbutton_EEG2Bipolar, 'CData', pbIml2r);
    set(handles.pushbutton_Bipolar2EEG, 'Cdata', pbImr2l);
    clear pbIml2r pbImr2l r g b;
else
    set(handles.pushbutton_EEG2Bipolar, 'String', 'E2B');
    set(handles.pushbutton_Bipolar2EEG, 'String', 'B2E');
end

% hide axes and create boxes around graphs
set(handles.axes_Signal, 'Box', 'on');
xlabel(handles.axes_Signal, 'Duration (s)');
ylabel(handles.axes_Signal, 'Current channel (µV)');
set(handles.axes_Frequency, 'Box', 'on',...
    'XTickLabel', [], 'YTickLabel', []);
xlabel(handles.axes_Frequency, 'Frequency (Hz)');
set(handles.axes_Channels, 'Box', 'on',...
    'XTickLabel', [], 'YTickLabel', []);
xlabel(handles.axes_Channels, 'Channel number');
ylabel(handles.axes_Channels, 'Amplitude (µV)');

set(handles.figure_EEGQC, 'Name', ['EEG Quality Control v. 5.0',...
    '                Helppi filua odotellessa: viljami.sairanen@hus.fi'],...
    'HitTest', 'off');

set(handles.listbox_EEG, 'Min', 0);
set(handles.listbox_Bipolar, 'Min', 0);
handles.minmaxZoom = false;

handles.data = [];
handles.time = [];
handles.plot_Signal = plot([],'Parent', handles.axes_Signal);
    
% Choose default command line output for EEGQC_v5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = EEGQC_v5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox_EEG.
function listbox_EEG_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_EEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listbox_EEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_EEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_Bipolar.
function listbox_Bipolar_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Bipolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listbox_Bipolar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Bipolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_EEG2Bipolar.
function pushbutton_EEG2Bipolar_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_EEG2Bipolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'EEGChannels')
    return;
end
if ~isempty(handles.EEGChannels)
    indSelected = get(handles.listbox_EEG, 'Value');
    newBip = [handles.EEGChannels(indSelected, :); handles.BipolarChannels];
    newBip = sortrows(newBip, 2);
    handles.EEGChannels(indSelected, :) = [];
    handles.BipolarChannels = newBip;
    set(handles.listbox_EEG, 'String', handles.EEGChannels(:,1));
    set(handles.listbox_Bipolar, 'String', handles.BipolarChannels(:,1));
    if indSelected(1) > length(handles.EEGChannels(:,1))
        set(handles.listbox_EEG, 'Value', length(handles.EEGChannels(:,1)));
    else
        set(handles.listbox_EEG, 'Value', indSelected(1));
    end
    set(handles.listbox_Bipolar, 'Value', length(handles.BipolarChannels(:,1)));    
end
% Choose default command line output for EEGQC_v5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Bipolar2EEG.
function pushbutton_Bipolar2EEG_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Bipolar2EEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'BipolarChannels')
    return;
end


if ~isempty(handles.BipolarChannels)
    indSelected = get(handles.listbox_Bipolar, 'Value');
    newEEG = [handles.EEGChannels; handles.BipolarChannels(indSelected, :)];
    newEEG = sortrows(newEEG, 2);
    handles.BipolarChannels(indSelected, :) = [];
    handles.EEGChannels = newEEG;
    set(handles.listbox_EEG, 'String', handles.EEGChannels(:,1));
    set(handles.listbox_Bipolar, 'String', handles.BipolarChannels(:,1));
    if indSelected(1) > length(handles.BipolarChannels(:,1))
        set(handles.listbox_Bipolar, 'Value', length(handles.BipolarChannels(:,1)));
    else
        set(handles.listbox_Bipolar, 'Value', indSelected(1));
    end
    set(handles.listbox_EEG, 'Value', length(handles.EEGChannels(:,1)));
end
% Choose default command line output for EEGQC_v5
handles.testi1 = 1;
%%handles.output = hObject;

% Update handles structure
guidata(hObject, handles);    

% --------------------------------------------------------------------
function uipushtool_Import_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.fileName, handles.filePath] = uigetfile('*.txt', 'Select an ASCII file');
if  isequal(handles.fileName, 0)
    return
end

if isfield(handles, 'a')
    % Remove areas from previous study.
    handles = rmfield(handles, 'a');
end

if ~strcmp(handles.filePath(end), filesep)
    handles.filePath(end+1) = filesep;
end
fullPath = strcat(handles.filePath, handles.fileName);

wbh = waitbar(0, 'Importing data...');

% Convert possible commas into points
fid = fopen(fullPath, 'r');
%tmp = fgetl(fid); % header line
%tmp = fgetl(fid); % hopefully the first data line
%if strfind(tmp,',')
    % test if file has comma as decimal separator and fix it if needed
    comma = uint8(',');
    point = uint8('.');
    convfile = memmapfile(fullPath, 'writable', true);
    convfile.Data(transpose(convfile.Data == comma)) = point;
    clear comma point convfile;
%end
fclose(fid);

fid = fopen(fullPath, 'r');
header = fgetl(fid);
header = textscan(header, '%s');
handles.header = header{1};

handles.nChannels = length(handles.header);
for i = 1:handles.nChannels
    handles.header(i,2) = {i};
end
dataformat = repmat('%f\t', 1, handles.nChannels);

% rows in ASCII file
frewind(fid);
handles.nRows = numel(cell2mat(textscan(fid, '%1c%*[^\n]')));

handles.data = zeros(handles.nRows, handles.nChannels, 'single');
rowCount = 0;
headerCount = 0;
frewind(fid);

while rowCount + headerCount < handles.nRows
    clear tmp tmp2
    tmp2 = textscan(fid, dataformat, 'Headerlines', 1);
    lens = unique(cellfun(@length, tmp2));
    if length(lens) > 1
        for i = 1:length(tmp2)
            if length(tmp2{i}) > min(lens)
                tmp2{i} = tmp2{i}(1:min(lens));
            end
        end
    end
    %{
    if length(tmp2{1}) ~= length(tmp2{2})
        tmp2{1} = tmp2{1}(1:end-1);
    end
    
    if length(tmp2) == 65
        len = cellfun(@length, tmp2);
        for i = 1:length(len)
            if len(i) == 0
                tmp2{i} = zeros(1,0);
            end
        end
        if max(len) ~= min(len)
            for i = 1:65
                tmp2{i} = [tmp2{i}; zeros(max(len)-length(tmp2{i}),1)];
            end
        end
    end
    %}
    %fprintf('tmp2\n')
    %assignin('base','tmp2',tmp2)
    %handles.data = handles.data(:,1:handles.nChannels);
    tmp = cell2mat(tmp2);
    
    headerCount = headerCount + 1;
    %if isempty(tmp)
    %    tmp = zeros(1,length(tmp2));
    %end
    %size(handles.data(rowCount+1:size(tmp,1)+rowCount, :))
    %size(tmp)
    handles.data(rowCount+1:size(tmp,1)+rowCount, :) = tmp;
    
    rowCount = rowCount + size(tmp,1);
    
    waitbar((rowCount+headerCount)/handles.nRows);
end
fclose(fid);
close(wbh);

handles.data(end-headerCount:end,:) = [];
assignin('base','handles',handles)
[handles.frequency, handles.time] = inputFrequency(length(handles.data(:,1)));
set(handles.textHz,'String',[num2str(handles.frequency), ' Hz']);
handles.EEGChannels = handles.header(1:floor(3*end/4-1),:);
handles.BipolarChannels = handles.header(floor(3*end/4):end,:);

set(handles.listbox_EEG, 'String', handles.EEGChannels(:,1),...
    'Max', handles.nChannels, 'Value', (1:length(handles.EEGChannels(:,1))));
set(handles.listbox_Bipolar, 'String', handles.BipolarChannels(:,1),...
    'Max', handles.nChannels, 'Value', (1:length(handles.BipolarChannels(:,1))));

handles.plot_minmaxPositions = zeros(handles.nChannels, 4);

% Plot
handles.plot_minVPPSignal = ...
    plot(handles.axes_Signal, -5, -5, 'rx', 'MarkerSize' ,10);
hold(handles.axes_Signal, 'on');
handles.plot_maxVPPSignal = ...
    plot(handles.axes_Signal, -5, -5, 'rx', 'MarkerSize' ,10);
handles.plot_minVPPSignalDurationLine = ...
    line([-5 -5], [-5 5], 'Parent', handles.axes_Signal,...
    'LineWidth', 1, 'LineStyle', ':', 'Color', 'red');
handles.plot_maxVPPSignalDurationLine = ...
    line([-5 -5], [-5 5], 'Parent', handles.axes_Signal,...
    'LineWidth', 1, 'LineStyle', ':', 'Color', 'red');
handles.plot_minVPPSignalWindowLine = ...
    line([-5 -5], [-5 5], 'Parent', handles.axes_Signal,...
    'LineWidth', 1, 'Marker', '+', 'Color', 'red');
handles.plot_maxVPPSignalWindowLine = ...
    line([-5 -5], [-5 5], 'Parent', handles.axes_Signal,...
    'LineWidth', 1, 'Marker', '+',  'Color', 'red');
handles.plot_Signal = plot(handles.time,...
    handles.data(:, handles.currentChannel),...
    'Parent', handles.axes_Signal);

%hold(handles.axes_Signal, 'on');

yMax = max(max(handles.data));
yMin = min(min(handles.data));

handles.origXlim = [handles.time(1), handles.time(end)];

for i = 1:99
    handles.plot_Signal_Areas(i) = ...
        patch([0 0 0 0], [yMax+[100, 100] yMin-[100, 100]],'red','FaceAlpha', 0.1,'Parent', handles.axes_Signal);
    %rectangle('Position',[-30 20 0 yMax-yMin],'Parent', handles.axes_Signal,'FaceAlpha', 0.1);
    %area([0, 0], yMax+[100, 100],'Parent', handles.axes_Signal, 'BaseValue', yMin-100);
    %tmp = get(handles.plot_Signal_Areas(i), 'Children');
    %set(tmp, 'FaceAlpha', 0.1);
end


hold(handles.axes_Signal, 'off');

xlabel(handles.axes_Signal, 'Duration (s)');
plotUpdate(hObject, handles, [handles.time(1), handles.time(end)]);


% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);     

% ------------------------------------------------------- plot update
function plotUpdate(hObject, handles, xLim,forceYLim)

if isempty(handles.data)
    return;
end

if nargin == 2
    xLim = get(handles.axes_Signal, 'XLim');
end

set(handles.plot_Signal, 'XData', handles.time,...
    'YData', handles.data(:, handles.currentChannel));

if (get(handles.checkbox_yAxisLock, 'Value') == 1)
    yLim = handles.YLim;
else
    yLim = [min(handles.data(:,handles.currentChannel)),...
        max(handles.data(:,handles.currentChannel))];
end

hold(handles.axes_Signal, 'on')

if isfield(handles, 'a')

    n = length(handles.a);
    if n > 99
        n = 99;
    end
    for i = 1:n
        set(handles.plot_Signal_Areas(i), 'XData',...
            [handles.a(i).startTime+[0, handles.a(i).duration] handles.a(i).startTime+[handles.a(i).duration,0]]);
        if handles.a(i).active
            % Active area has green color and it's specs are shown in right
            % top corner of gui.
            set(handles.plot_Signal_Areas(i), 'FaceColor', 'Green');            
            set(handles.text_StartTime, 'String', ['Start time (s): ',...
                num2str(handles.a(i).startTime)]);
            set(handles.text_Duration, 'String', ['Duration (s): ',...
                num2str(handles.a(i).duration)]);
            set(handles.text_Window, 'String', ['Window (s): ',...
                num2str(handles.a(i).window)]);                        
        else
            % Not active areas have blue color.
            set(handles.plot_Signal_Areas(i), 'FaceColor', 'Blue');
        end    
    end
    
    
end
hold(handles.axes_Signal, 'off')

ylabel(handles.axes_Signal, ['CH ',...
    '\bf' handles.header{handles.currentChannel}, '\rm (µV)']);

if nargin == 4
    set(handles.axes_Signal, 'XLim', xLim, 'YLim', forceYLim);
    %handles.YLim = yLim;
    %handles.XLim = xLim;
else
    set(handles.axes_Signal, 'XLim', xLim, 'YLim', yLim);
    handles.YLim = yLim;
    handles.XLim = xLim;
end

% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles)

% --------------------------------------------------------------------
function uipushtool_frequency_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.data)
    return;
end

[handles.frequency, handles.time] = inputFrequency(length(handles.data(:,1)));
set(handles.textHz,'String',[num2str(handles.frequency), ' Hz']);
% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
plotUpdate(hObject, handles, [handles.time(1), handles.time(end)]);

% --------------------------------------------------------------------
function [freq, t] = inputFrequency(N)

% ask for user to give sampling frequency
freq = str2double(inputdlg('Enter frequency (Hz).',...
    'EEGQC', 1, {'250'}));
while isnan(freq) || ~(freq > 0)
    freq = str2double(inputdlg('Frequency input must be numeric and larger than zero.',...
    'EEGQC', 1, {'250'}));
end
t = (1/freq)*(1:N);


function axes_Signal_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_EEGQC_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure_EEGQC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get current point and axes limits
if isempty(handles.data)
    return
end
cp = get(handles.axes_Signal, 'CurrentPoint');
xLim = get(handles.axes_Signal, 'XLim');
yLim = get(handles.axes_Signal, 'YLim');
handles.buttonDown = true;
switch get(hObject, 'SelectionType')
    case 'normal' % left click
        % check that it lies on handles.axes_Signal
        if (cp(1,1) > xLim(1)) && (cp(1,1) < xLim(2)) &&...
                (yLim(1) < cp(1,2)) && (cp(1,2) < yLim(2))
            if ~isfield(handles, 'a')
                % if there is zero areas, make one
                if cp(1,1) + 20 > xLim(2)
                    cp(1,1) = xLim(2)-20;
                end
                handles.a(1) = analyzeArea(floor(cp(1,1)));
                activeAreaNumber = 1;
            else
                % else check if an area was clicked
                n = length(handles.a);
                hit = zeros(n,1);
                for i = 1:n
                    handles.a(i).active = false;
                    if (handles.a(i).startTime < cp(1,1)) &&...
                            (cp(1,1) < handles.a(i).startTime + handles.a(i).duration)
                        handles.a(i).active = true;
                        activeAreaNumber = i;
                        hit(i) = true;
                    else
                        hit(i) = false;
                    end                   
                end
                if all(~hit)
                    % area was not clicked, add new area
                    if cp(1,1) + 20 > xLim(2)
                         cp(1,1) = xLim(2)-20;
                    end
                    handles.a(end+1) = analyzeArea(floor(cp(1,1)));
                    activeAreaNumber = n+1;
                end
            end
        handles.activeAreaNumber = activeAreaNumber;
        handles.pickPoint = cp;
        handles.diffPickStart = handles.a(handles.activeAreaNumber).startTime - cp(1,1);
        guidata(hObject, handles);
        plotUpdate(hObject, handles,xLim);

        end            
    case 'alt' %right click 
        % nothing        
    case 'open' % double click
    otherwise
        % nothing    
end

% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_EEGQC_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure_EEGQC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.buttonDown = false;
% show results for the active area
    if isfield(handles, 'a');
        %n = length(handles.a);
        %for i = 1:n
        i = handles.activeAreaNumber;
        iArea = (0:handles.a(i).duration*handles.frequency)+...
            handles.a(i).startTime*handles.frequency;
        %assignin('base','data',handles.data);
        %assignin('base','iArea',iArea);
        if iArea(1) == 0
            iArea=iArea(2:end);
        end
        iData = handles.data(iArea, :);
        iTime = handles.time(iArea);

        % Rolling window with 1s step size
        windowsInArea = handles.a(i).duration-handles.a(i).window+1;

        VPP = zeros(windowsInArea, size(handles.data,2));

        % The first window
        k = 1;
        iWindow = 1:handles.a(i).window*handles.frequency;

        N = length(iData);
        if iWindow(end) > N
            % Example case: duration is 20s and window is 11s. The
            % first window (the first 11s) is ok, but second window
            % would be form 12s to 21s. Since this goes over the duration
            % of area it might pick some unwanted errors from the
            % data. Hence the second window is fixed to 20s-11s =
            % 9s...20s.
            iWindow = N - handles.a(i).window*handles.frequency:N;
        end
        % Find the local max and min of every channel in current
        % window and save it in VPP array.
  
        %[vppMinVal, vppMinInd] = min(iData(iWindow,:));
        %[vppMaxVal, vppMaxInd] = max(iData(iWindow,:));
        %minTime = iTime(iWindow(vppMinInd));
        %maxTime = iTime(iWindow(vppMaxInd));
        %VPP(k, :) = vppMaxVal - vppMinVal;
        VPP(k, :) = max(iData(iWindow,:)) - min(iData(iWindow,:));
        % The rest of windows
        for k = 2:windowsInArea
            iWindow = iWindow + 1*handles.frequency;
            %[newMinVal, newMinInd] = min(iData(iWindow,:));
            %[newMaxVal, newMaxInd] = max(iData(iWindow,:));
            %VPP(k, :) = newMaxVal - newMinVal;
            VPP(k, :) = max(iData(iWindow,:)) - min(iData(iWindow,:));
            
            %for ind = 1:length(VPP(k, :))
            %    if VPP(k, ind) < any(VPP(VPP(:,ind) ~= 0))
            %        vppMinVal(ind) = newMinVal(ind);
            %        newTime = iTime(iWindow(newMinInd));
            %        minTime(ind) = newTime(ind);
            %    end
            %    if VPP(k,ind) > any(VPP(VPP(:,ind) ~= 0))
            %        vppMaxVal(ind) = newMaxVal(ind);
            %        newTime = iTime(iWindow(newMaxInd));
            %        maxTime(ind) = newTime(ind);
            %    end
            %end               

        end
        
        handles.VPP = VPP;
        % VPP kolumnit vastaavat kanavia
        % VPP rivit vastaavat Windoweja
        %[maxval, maxind] = max(VPP); % maxind = window mistä maksimi löytyi
        %[minval, minind] = min(VPP); % minind = window mistä minimi löytyi
        % kanava kerrallaan kerätään minimi VPP max ja min positio sekä
        % maksimi VPP max ja min positio.
          
   
        % Get selected EEG and Bipolar channels from listboxes
        selectedEEGChannels = get(handles.listbox_EEG, 'Value');
        selectedBipolarChannels = get(handles.listbox_Bipolar, 'Value');
        if selectedEEGChannels == 0
            iEEGChannels = [];
        else
            iEEGChannels = cell2mat(handles.EEGChannels(...
                selectedEEGChannels, 2));
        end
        if selectedBipolarChannels == 0
            iBipolarChannels = [];
        else
            iBipolarChannels = cell2mat(handles.BipolarChannels(...
                selectedBipolarChannels, 2));
        end

        % FFT of the selected area
        
        iFT = fft(nansum(iData(:, [iEEGChannels; iBipolarChannels]), 2));
        iFT = iFT(1:floor(N/2));
        iFT(1:2) = 0; % 1:2 are error prone elements and unnecessary in this case.
        [maxValue, maxIndex] = max(abs(iFT));
        iFT = abs(iFT)./maxValue;
        frequencyAxis = handles.frequency*(0:N/2-1)./N;

        % Plot FFT data
        axes(handles.axes_Frequency);
        
        plot(frequencyAxis, iFT);
        text(frequencyAxis(maxIndex), 0.9, ['\leftarrow ',...
            num2str(frequencyAxis(maxIndex)), ' Hz'],...
            'Parent', handles.axes_Frequency);
        xlabel(handles.axes_Frequency, 'Frequency (Hz)');
        %freqLim = frequencyAxis(maxIndex)+[-15, +15];
        %if freqLim(1) < 0
        %    freqLim(1) = 0;
        %end
        %if freqLim(2) > frequencyAxis(end)
        %    freqLim(2) = frequencyAxis(end);
        %end
        %set(handles.axes_Frequency, 'XLim', freqLim);           

        % Plot VPP data  
        plot(iEEGChannels, max(VPP(:, iEEGChannels)), 'b^',...
            'MarkerSize', 7, 'MarkerFaceColor', 'b',...
            'Parent', handles.axes_Channels);
        hold(handles.axes_Channels, 'on');
        plot(iEEGChannels, min(VPP(:, iEEGChannels)), 'bv',...
            'MarkerSize', 7, 'MarkerFaceColor', 'none',...
            'Parent', handles.axes_Channels);
        plot(iBipolarChannels, max(VPP(:, iBipolarChannels)), 'r^',...
            'MarkerSize', 7, 'MarkerFaceColor', 'r',...
            'Parent', handles.axes_Channels);
        plot(iBipolarChannels, min(VPP(:, iBipolarChannels)), 'rv',...
            'MarkerSize', 7, 'MarkerFaceColor', 'none',...
            'Parent', handles.axes_Channels);

        handles.plot_channelLocations =...
            [iEEGChannels', iBipolarChannels';...
            max(VPP(:, iEEGChannels)), max(VPP(:, iBipolarChannels));...
            min(VPP(:, iEEGChannels)), min(VPP(:, iBipolarChannels))];

        line([0, handles.nChannels+1],...
            mean(max(VPP(:, iEEGChannels))).*[1, 1],...
            'Color', 'Blue', 'Parent', handles.axes_Channels);
        line([0, handles.nChannels+1],...
            mean(min(VPP(:, iEEGChannels))).*[1, 1],...
            'Color', 'Blue', 'Parent', handles.axes_Channels,...
            'LineStyle', '--');
        line([0, handles.nChannels+1],...
            mean(max(VPP(:, iBipolarChannels))).*[1, 1],...
            'Color', 'Red', 'Parent', handles.axes_Channels);
        line([0, handles.nChannels+1],...
            mean(min(VPP(:, iBipolarChannels))).*[1, 1],...
            'Color', 'Red', 'Parent', handles.axes_Channels,...
            'LineStyle', '--');
        handles.plot_activeChannel =...
            plot(handles.axes_Channels, -5, -5, 'go',...
            'MarkerSize', 10, 'LineWidth', 2);

        hold(handles.axes_Channels, 'off');
        set(handles.axes_Channels, 'XLim', [0, handles.nChannels+1],...
            'XMinorTick', 'on', 'XGrid', 'on');
        ylabel(handles.axes_Channels, 'Amplitude (µV)');
        xlabel(handles.axes_Channels, 'Channel number');
        legend(handles.axes_Channels,...
            'EEG max', 'EEG min', 'Bipolar max', 'Bipolar min',...
            'Mean EEG max', 'Mean EEG min', 'Mean Bipolar Max',...
            'Mean Bipolar Min', 'Location', 'NorthEastOutside');     
        %end
    end
% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%eventdata.Key = 'return';
%figure_EEGQC_WindowKeyPressFcn(hObject, eventdata, handles);

% --- Executes on mouse motion over figure - except title and menu.
function figure_EEGQC_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure_EEGQC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 if ~isempty(get(handles.axes_Signal,'UserData'))
     return;
 end

% Move analyze areas
if isfield(handles, 'buttonDown') && isfield(handles, 'a')
    if handles.buttonDown
        cp = get(handles.axes_Signal, 'CurrentPoint');
        xLim = get(handles.axes_Signal, 'XLim');

        if cp(1,1) + handles.diffPickStart < xLim(1)
            cp(1,1) = xLim(1) - handles.diffPickStart;
        end
        if cp(1,1) + handles.a(handles.activeAreaNumber).duration +...
                handles.diffPickStart > xLim(2)
            cp(1,1) = xLim(2)-handles.a(handles.activeAreaNumber).duration-...
                handles.diffPickStart;
        end
        
        handles.a(handles.activeAreaNumber).startTime = round(cp(1,1)+...
            handles.diffPickStart);
        
        
        % Choose default command line output for EEGQC_v5
        %handles.output = hObject;
        % Update handles structure
        guidata(hObject, handles);
        plotUpdate(hObject, handles, xLim);
        
    end
end
% Zoom into channel min/max points

if isfield(handles, 'plot_channelLocations') && isfield(handles,'a')    
    cp = get(handles.axes_Channels, 'CurrentPoint');
    channelAxesXLim = get(handles.axes_Channels, 'XLim');
    channelAxesYLim = get(handles.axes_Channels, 'YLim');
    if cp(1,1) > channelAxesXLim(1) && cp(1,1) < channelAxesXLim(2) &&...
            cp(1,2) > channelAxesYLim(1) && cp(1,2) < channelAxesYLim(2)
        % check if nearest point is closer than 0.5
        [xdist, xmini] = min(abs(cp(1,1)-handles.plot_channelLocations(1,:)));
        [ymaxdist, ymaxind] = min(abs(cp(1,2) - handles.plot_channelLocations(2,xmini)));
        [ymindist, yminind] = min(abs(cp(1,2) - handles.plot_channelLocations(3,xmini)));
        
        if xdist < 0.5 && ymaxdist < 0.5 || ymindist < 0.5
            tmpdif = abs(handles.plot_channelLocations(2,xmini)-...
                handles.plot_channelLocations(3,xmini));
            if tmpdif < 1
                tmpdif = 1;
            end
            if ymaxdist < tmpdif/2               
                [maxval, maxind] = max(handles.VPP);
                iWindow = (1:handles.a(handles.activeAreaNumber).window*handles.frequency)+...
                    (maxind(xmini)-1)*handles.frequency+...
                    handles.a(handles.activeAreaNumber).startTime*handles.frequency;                
                [maxval, maxind] = max(handles.data(iWindow, xmini));
                [minval, minind] = min(handles.data(iWindow, xmini));
                maxind = maxind + iWindow(1)-1;
                minind = minind + iWindow(1)-1;

                % Show green circle around selected area min/max
                set(handles.plot_activeChannel,...
                    'XData', handles.plot_channelLocations(1,xmini),...
                    'YData', handles.plot_channelLocations(2,xmini));
                % plot minimum and maximum of maximum VPP of selected
                % channel in selected area.
                set(handles.plot_minVPPSignalDurationLine,...
                    'XData', handles.a(handles.activeAreaNumber).startTime+...
                    [0, handles.a(handles.activeAreaNumber).duration],...
                    'YData', minval.*[1 1]);
                set(handles.plot_maxVPPSignalDurationLine,...
                    'XData', handles.a(handles.activeAreaNumber).startTime+...
                    [0, handles.a(handles.activeAreaNumber).duration],...
                    'YData', maxval.*[1 1]);
                set(handles.plot_minVPPSignalWindowLine,...
                    'XData', [iWindow(1), iWindow(end)]/handles.frequency,...
                    'YData', minval.*[1 1]);
                set(handles.plot_maxVPPSignalWindowLine,...
                    'XData', [iWindow(1), iWindow(end)]/handles.frequency,...
                    'YData', maxval.*[1 1]);
                set(handles.plot_minVPPSignal,...
                    'XData', handles.time(minind),...
                    'YData', minval);
                set(handles.plot_maxVPPSignal,...
                    'XData', handles.time(maxind),...
                    'YData', maxval);
            end
            if ymindist < tmpdif/2
                [minval, minind] = min(handles.VPP);
                iWindow = (1:handles.a(handles.activeAreaNumber).window*handles.frequency)+...
                    (minind(xmini)-1)*handles.frequency+...
                    handles.a(handles.activeAreaNumber).startTime*handles.frequency;
                [maxval, maxind] = max(handles.data(iWindow, xmini));
                [minval, minind] = min(handles.data(iWindow, xmini));
                maxind = maxind + iWindow(1)-1;
                minind = minind + iWindow(1)-1;
                % Show green circle around selected area min/max
                set(handles.plot_activeChannel,...
                    'XData', handles.plot_channelLocations(1,xmini),...
                    'YData', handles.plot_channelLocations(3,xmini));       
                % plot minimum and maximum of minimum VPP of selected
                % channel in selected area.
                set(handles.plot_minVPPSignalDurationLine,...
                    'XData', handles.a(handles.activeAreaNumber).startTime+...
                    [0, handles.a(handles.activeAreaNumber).duration],...
                    'YData', minval.*[1 1]);
                set(handles.plot_maxVPPSignalDurationLine,...
                    'XData', handles.a(handles.activeAreaNumber).startTime+...
                    [0, handles.a(handles.activeAreaNumber).duration],...
                    'YData', maxval.*[1 1]);
                set(handles.plot_minVPPSignalWindowLine,...
                    'XData', [iWindow(1), iWindow(end)]/handles.frequency,...
                    'YData', minval.*[1 1]);
                set(handles.plot_maxVPPSignalWindowLine,...
                    'XData', [iWindow(1), iWindow(end)]/handles.frequency,...
                    'YData', maxval.*[1 1]);
                set(handles.plot_minVPPSignal,...
                    'XData', handles.time(minind),...
                    'YData', minval);
                set(handles.plot_maxVPPSignal,...
                    'XData', handles.time(maxind),...
                    'YData', maxval);
            end
            set(handles.plot_Signal, 'YData', handles.data(:, xmini));
            xLim = handles.a(handles.activeAreaNumber).startTime+ [0,...
                handles.a(handles.activeAreaNumber).duration];
            yLim = [min(handles.data(xLim(1)*handles.frequency:xLim(2)*handles.frequency,xmini))*1.1,...
                max(handles.data(xLim(1)*handles.frequency:xLim(2)*handles.frequency,xmini))*1.1];
            set(handles.axes_Signal, 'XLim', xLim, 'YLim', yLim);
            ylabel(handles.axes_Signal, ['Ch ', ...
                '\bf' handles.header{xmini}, '\rm (µV)']);

            handles.minmaxZoom = true;
        end
    else
        if isfield(handles, 'plot_activeChannel')
            % Hide lines on signal axes and green circle on channel axes
            set(handles.plot_activeChannel,...
                'XData', -5,...
                'YData', -5);
                set(handles.plot_minVPPSignalDurationLine,...
                    'XData', -10.*[1 1]);
                set(handles.plot_maxVPPSignalDurationLine,...
                    'XData', -10.*[1 1]);
                set(handles.plot_minVPPSignalWindowLine,...
                    'XData', -10.*[1 1]);
                set(handles.plot_maxVPPSignalWindowLine,...
                    'XData', -10.*[1 1]);
                set(handles.plot_minVPPSignal,...
                    'XData', -5,...
                    'YData', -5);
                set(handles.plot_maxVPPSignal,...
                    'XData', -5,...
                    'YData', -5);
        end
        if handles.minmaxZoom
            set(handles.plot_Signal, 'YData', handles.data(:, handles.currentChannel));
            set(handles.axes_Signal, 'XLim', handles.XLim);
            set(handles.axes_Signal, 'YLim', handles.YLim);   
            ylabel(handles.axes_Signal, ['CH ',...
                '\bf' handles.header{handles.currentChannel}, '\rm (µV)']);
            handles.minmaxZoom = false;
        end
    end
end

% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes on scroll wheel click while the figure is in focus.
function figure_EEGQC_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure_EEGQC (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'activeAreaNumber')
    return;
end

active = handles.activeAreaNumber;
cp = get(handles.axes_Signal, 'CurrentPoint');
xLim = get(handles.axes_Signal, 'XLim');
yLim = get(handles.axes_Signal, 'YLim');

if (cp(1,1) > xLim(1)) && (cp(1,1) < xLim(2)) &&...
    (yLim(1) < cp(1,2)) && (cp(1,2) < yLim(2))
    switch get(hObject, 'SelectionType')
        case 'normal'
            %% Hiiren vasemman napin painon jälkeen scroll muuttaa Durationia
            if  eventdata.VerticalScrollCount < 0
                if handles.a(active).duration < 60
                    handles.a(active) = handles.a(active).changeDuration(+1);
                    if handles.a(active).startTime + handles.a(active).duration > handles.time(end)+1
                        handles.a(active) = handles.a(active).changeStartTime(-1);
                    end
                end
            elseif  eventdata.VerticalScrollCount > 0
                if handles.a(active).duration > 1
                    handles.a(active) = handles.a(active).changeDuration(-1);
                end
                if handles.a(active).duration < handles.a(active).window
                    handles.a(active).window = handles.a(active).duration;
                end
            end

        case 'alt'
            %% Hiiren oikean napin painon jälkeen scoll muuttaa Windowin pituutta
            if  eventdata.VerticalScrollCount < 0
                if handles.a(active).window < floor(handles.a(active).duration-1)
                    handles.a(active) = handles.a(active).changeWindow(+1);
                end
            elseif  eventdata.VerticalScrollCount > 0
                if handles.a(active).window > 1
                    handles.a(active) = handles.a(active).changeWindow(-1);
                end
            end
        otherwise
    end
    %guidata(hObject,handles);
    %plotUpdate(hObject, handles, get(handles.axes_Signal, 'XLim'));
    guidata(hObject, handles);
plotUpdate(hObject, handles, get(handles.axes_Signal, 'XLim'));

end

%set(handles.plot_Signal_Areas(active), 'XData',...
%    handles.a(active).startTime+[0, handles.a(active).duration]);

% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
%plotUpdate(hObject, handles, handles.origXlim);

%figure_EEGQC_WindowButtonUpFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on figure_EEGQC or any of its controls.
function figure_EEGQC_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_EEGQC (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.data)
    return;
end

xLim = get(handles.axes_Signal, 'XLim');
xDif = xLim(2)-xLim(1);
%yLim = get(handles.axes_Signal, 'YLim');

switch eventdata.Key
    case 'return'
        %{
        % show results for the active area
        if isfield(handles, 'a');
            %n = length(handles.a);
            %for i = 1:n
            i = handles.activeAreaNumber;
            iArea = (0:handles.a(i).duration*handles.frequency)+...
                handles.a(i).startTime*handles.frequency;
            iData = handles.data(iArea, :);
            
            % Rolling window with 1s step size
            windowsInArea = handles.a(i).duration-handles.a(i).window+1;

            VPP = zeros(windowsInArea, size(handles.data,2));
            
            % The first window
            k = 1;
            iWindow = 1:handles.a(i).window*handles.frequency;

            N = length(iData);
            if iWindow(end) > N
                % Example case: duration is 20s and window is 11s. The
                % first window (the first 11s) is ok, but second window
                % would be form 12s to 21s. Since this goes over the duration
                % of area it might pick some unwanted errors from the
                % data. Hence the second window is fixed to 20s-11s =
                % 9s...20s.
                iWindow = N - handles.a(i).window*handles.frequency:N;
            end
            % Find the local max and min of every channel in current
            % window and save it in VPP array.
            
            % --------------------
            % uncomment block below to see Rolling window in action!
            
            %{
            nF = figure;
            nA = axes('Parent', nF);            
            plot(nA, handles.time(iArea(iWindow)), iData(iWindow,1));
            hold(nA, 'on');
            [maxv, maxi] = max(iData(iWindow,1));
            [minv, mini] = min(iData(iWindow,1)); 
            plot(handles.time(iArea(iWindow(maxi))),...
                iData(iWindow(maxi),1),'r^');
            plot(handles.time(iArea(iWindow(mini))),...
                iData(iWindow(mini),1),'rv');
            title(nA, ['VPP: ' num2str(maxv-minv)]);
            hold(nA, 'off');
            pause
            %}
           
            
            VPP(k,:) = max(iData(iWindow,:)) - min(iData(iWindow,:));
            % The rest of windows
            for k = 2:windowsInArea
                iWindow = iWindow + 1*handles.frequency;
                VPP(k,:) = max(iData(iWindow,:)) - min(iData(iWindow,:));
                
                % --------------------
                % uncomment block below to see Rolling window in action!                
                %{
                plot(nA, handles.time(iArea(iWindow)), iData(iWindow,1));
                hold(nA, 'on');
                [maxv, maxi] = max(iData(iWindow,1));
                [minv, mini] = min(iData(iWindow,1));                      
                plot(handles.time(iArea(iWindow(maxi))),...
                    iData(iWindow(maxi),1),'r^');
                plot(handles.time(iArea(iWindow(mini))),...
                    iData(iWindow(mini),1),'rv');
                title(nA, ['VPP: ' num2str(maxv-minv)]);
                hold(nA, 'off');
                pause
                %}
                
            end

            % Get selected EEG and Bipolar channels from listboxes
            selectedEEGChannels = get(handles.listbox_EEG, 'Value');
            selectedBipolarChannels = get(handles.listbox_Bipolar, 'Value');
            if selectedEEGChannels == 0
                iEEGChannels = [];
            else
                iEEGChannels = cell2mat(handles.EEGChannels(...
                    selectedEEGChannels, 2));
            end
            if selectedBipolarChannels == 0
                iBipolarChannels = [];
            else
                iBipolarChannels = cell2mat(handles.BipolarChannels(...
                    selectedBipolarChannels, 2));
            end
            
            % FFT of the selected area
            iFT = fft(sum(iData(:, [iEEGChannels; iBipolarChannels]), 2));
            iFT = iFT(1:floor(N/2));
            iFT(1:2) = 0; % 1:2 are error prone elements and unnecessary in this case.
            [maxValue, maxIndex] = max(abs(iFT));
            iFT = abs(iFT)./maxValue;
            frequencyAxis = handles.frequency*(0:N/2-1)./N;
            
            % Plot FFT data
            plot(frequencyAxis, iFT, 'Parent', handles.axes_Frequency);
            text(frequencyAxis(maxIndex), 0.9, ['\leftarrow ',...
                num2str(frequencyAxis(maxIndex)), ' Hz'],...
                'Parent', handles.axes_Frequency);
            xlabel(handles.axes_Frequency, 'Frequency (Hz)');
            %freqLim = frequencyAxis(maxIndex)+[-15, +15];
            %if freqLim(1) < 0
            %    freqLim(1) = 0;
            %end
            %if freqLim(2) > frequencyAxis(end)
            %    freqLim(2) = frequencyAxis(end);
            %end
            %set(handles.axes_Frequency, 'XLim', freqLim);               
            
            % Plot VPP data  
                plot(iEEGChannels, max(VPP(:, iEEGChannels)), 'b^',...
                'MarkerSize', 7, 'MarkerFaceColor', 'b',...
                'Parent', handles.axes_Channels);
            hold(handles.axes_Channels, 'on');
            plot(iEEGChannels, min(VPP(:, iEEGChannels)), 'bv',...
                'MarkerSize', 7, 'MarkerFaceColor', 'none',...
                'Parent', handles.axes_Channels);
            plot(iBipolarChannels, max(VPP(:, iBipolarChannels)), 'r^',...
                'MarkerSize', 7, 'MarkerFaceColor', 'r',...
                'Parent', handles.axes_Channels);
            plot(iBipolarChannels, min(VPP(:, iBipolarChannels)), 'rv',...
                'MarkerSize', 7, 'MarkerFaceColor', 'none',...
                'Parent', handles.axes_Channels);
            
            handles.plot_channelLocations =...
                [iEEGChannels', iBipolarChannels';...
                max(VPP(:, iEEGChannels)), max(VPP(:, iBipolarChannels));...
                min(VPP(:, iEEGChannels)), min(VPP(:, iBipolarChannels))];
            
            line([0, handles.nChannels+1],...
                mean(max(VPP(:, iEEGChannels))).*[1, 1],...
                'Color', 'Blue', 'Parent', handles.axes_Channels);
            line([0, handles.nChannels+1],...
                mean(min(VPP(:, iEEGChannels))).*[1, 1],...
                'Color', 'Blue', 'Parent', handles.axes_Channels,...
                'LineStyle', '--');
            line([0, handles.nChannels+1],...
                mean(max(VPP(:, iBipolarChannels))).*[1, 1],...
                'Color', 'Red', 'Parent', handles.axes_Channels);
            line([0, handles.nChannels+1],...
                mean(min(VPP(:, iBipolarChannels))).*[1, 1],...
                'Color', 'Red', 'Parent', handles.axes_Channels,...
                'LineStyle', '--');
            handles.plot_activeChannel =...
                plot(handles.axes_Channels, -5, -5, 'go',...
                'MarkerSize', 10, 'LineWidth', 2);       
                       
            hold(handles.axes_Channels, 'off');
            set(handles.axes_Channels, 'XLim', [0, handles.nChannels+1],...
                'XMinorTick', 'on', 'XGrid', 'on');
            ylabel(handles.axes_Channels, 'Amplitude (µV)');
            xlabel(handles.axes_Channels, 'Channel number');
            legend(handles.axes_Channels,...
                'EEG max', 'EEG min', 'Bipolar max', 'Bipolar min',...
                'Mean EEG max', 'Mean EEG min', 'Mean Bipolar Max',...
                'Mean Bipolar Min', 'Location', 'NorthEastOutside');     
            %end
        end
        %}
    case 'leftarrow'
        if xLim(1) > handles.time(1)
            xLim(1) = xLim(1)-0.3*xDif;
        end
        if xLim(1) < handles.time(1)
            xLim(1) = handles.time(1);
        end
        xLim(2) = xLim(1)+xDif;
        
        set(handles.axes_Signal, 'XLim', xLim);
        
    case 'rightarrow'
        if xLim(2) < handles.time(end)
            xLim(2) = xLim(2)+0.3*xDif;
        end
        if xLim(2) > handles.time(end)
            xLim(2) = handles.time(end);
        end
        xLim(1) = xLim(2)-xDif;
        set(handles.axes_Signal, 'XLim', xLim);       
    case 'uparrow'
        if handles.currentChannel > 1
            handles.currentChannel = handles.currentChannel - 1;
        end
        guidata(hObject, handles);
        plotUpdate(hObject, handles,get(handles.axes_Signal,'XLim'),get(handles.axes_Signal,'YLim'));
    case 'downarrow'
        if handles.currentChannel < size(handles.data,2)
            handles.currentChannel = handles.currentChannel + 1;
        end
        guidata(hObject, handles);
        plotUpdate(hObject, handles,get(handles.axes_Signal,'XLim'),get(handles.axes_Signal,'YLim'));
    case 'delete'
        if isfield(handles, 'a')
            for i = 1:length(handles.a)
                if handles.a(i).active
                    handles.a(i) = [];
                    set(handles.plot_Signal_Areas(i), 'XData', [0,0,0,0]);
                    tmp = handles.plot_Signal_Areas(i);
                    handles.plot_Signal_Areas(i:end-1) =...
                        handles.plot_Signal_Areas(i+1:end);
                    handles.plot_Signal_Areas(end) = tmp;
                    
                    % remove area specs texts
                    set(handles.text_StartTime, 'String', 'Start time (s):');
                    set(handles.text_Duration, 'String', 'Duration (s):');                        
                    set(handles.text_Window, 'String', 'Window (s):');
                    if isempty(handles.a)
                        handles=rmfield(handles,'a');                        
                    end
                    break;
                end
            end
            guidata(hObject, handles);
            plotUpdate(hObject, handles, xLim);
        end
    otherwise
end

% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure




% --- Executes on button press in checkbox_yAxisLock.
function checkbox_yAxisLock_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_yAxisLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_yAxisLock

if ~isfield(handles,'header')
    return;
end

handles.YLim = get(handles.axes_Signal, 'YLim');

if get(hObject, 'Value') == 1
    str = strcat('<html> Locked on channel <b>',...
        handles.header{handles.currentChannel},...
        '</b> y-axis</html>');
    set(handles.checkbox_yAxisLock, 'String', str);
else
    set(handles.checkbox_yAxisLock, 'String', 'Lock current y-axis');
end
% Choose default command line output for EEGQC_v5
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool_Report_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_Report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Choose default command line output for EEGQC_v5
if ~isfield(handles, 'a');
    return;
end

% open a word file
[outPutFile, outPutPath] = uiputfile('*.xls', 'Excel report name');
if outPutFile==0
    return;
end

prog = waitbar(0,'Building report');
set(handles.axes_Signal,'UserData',1)

if ~strcmp(outPutPath(end), filesep)
    outPutPath(end+1) = filesep;
end

outPutFull = strcat(outPutPath, outPutFile);

% Run Excel
excelApplication = actxserver('Excel.Application');
%set(excelApplication, 'Visible', 1);

if ~exist(outPutFull, 'file')
    % Open new Excel file
    excelWorkbook = invoke(excelApplication.Workbooks, 'Add');
    excelSheets = excelApplication.ActiveWorkBook.Sheets;
    % Remove initial 2 extra sheets
    %excelSheet = get(excelSheets, 'Item', 3);
    %invoke(excelSheet, 'Delete');
    %excelSheet = get(excelSheets, 'Item', 2);
    %invoke(excelSheet, 'Delete');
else
    % Open existing Excel file
    excelWorkbook = invoke(excelApplication.Workbooks, 'Open', outPutFull);
    excelSheets = excelApplication.ActiveWorkBook.Sheets;
end

% Get last sheet and test if it's first element is clear
excelSheet = get(excelSheets, 'Item', 1);
excelSheetRange = get(excelSheet, 'Range', 'A1');
while ~isnan(get(excelSheetRange, 'Value'))
    invoke(excelSheets, 'Add');
    excelSheet = get(excelSheets, 'Item', 1);
    excelSheetRange = get(excelSheet, 'Range', 'A1');
end

excelSheet.Name = strcat(date, '(', num2str(excelSheets.Count), ')');

waitbar(0.2,prog)
excelSheetRange = get(excelSheet, 'Range', 'F2');
excelSheetRange.ColumnWidth = 4;
excelSheetRange = get(excelSheet, 'Range', 'G2');
excelSheetRange.ColumnWidth = 120;
excelSheetRange = get(excelSheet, 'Range', 'A2:E2');
excelSheetRange.ColumnWidth = 15;
set(excelSheetRange, 'Value', {'~Frequency (Hz)', 'EEG max (µV)',...
    'EEG min (µV)', 'Bipolar max (µV)' , 'Bipolar min (µV)'});
set(excelSheetRange.Font, 'Bold', true);

tmpFigure = figure('Position', [0 0 1000 450], 'visible', 'off',...
    'Color', [38/255 155/255 158/255], 'InvertHardcopy', 'off');

n = length(handles.a);
for i = 1:n
    
    excelColumns = strcat('A', num2str(i+2), ':E', num2str(i+2));
    excelSheetRange = get(excelSheet, 'Range', excelColumns);
    excelSheetRange.NumberFormat = '0,00';
    if mod(i,2) == 0
        excelSheetRange.Interior.ColorIndex = 24;
    end
    iArea = (0:handles.a(i).duration*handles.frequency)+...
        handles.a(i).startTime*handles.frequency;
    iData = handles.data(iArea, :);
    
    % Rolling window with 1s step size
    windowsInArea = handles.a(i).duration-handles.a(i).window+1;
    
    VPP = zeros(windowsInArea, size(handles.data,2));
    
    % The first window
    k = 1;
    iWindow = 1:handles.a(i).window*handles.frequency;
    
    N = length(iData);
    if iWindow(end) > N
        % Example case: duration is 20s and window is 11s. The
        % first window (the first 11s) is ok, but second window
        % would be form 12s to 21s. Since this goes over the duration
        % of area it might pick some unwanted errors from the
        % data. Hence the second window is fixed to 20s-11s =
        % 9s...20s.
        iWindow = N - handles.a(i).window*handles.frequency:N;
    end
    % Find the local max and min of every channel in current
    % window and save it in VPP array.
    
    VPP(k,:) = max(iData(iWindow,:)) - min(iData(iWindow,:));
    % The rest of windows
    for k = 2:windowsInArea
        iWindow = iWindow + 1*handles.frequency;
        VPP(k,:) = max(iData(iWindow,:)) - min(iData(iWindow,:));
    end
    
    % Get selected EEG and Bipolar channels from listboxes
    selectedEEGChannels = get(handles.listbox_EEG, 'Value');
    selectedBipolarChannels = get(handles.listbox_Bipolar, 'Value');
    if selectedEEGChannels == 0
        iEEGChannels = [];
    else
        iEEGChannels = cell2mat(handles.EEGChannels(...
            selectedEEGChannels, 2));
    end
    if selectedBipolarChannels == 0
        iBipolarChannels = [];
    else
        iBipolarChannels = cell2mat(handles.BipolarChannels(...
            selectedBipolarChannels, 2));
    end
    
    % Select active channel
    set(handles.plot_Signal_Areas, 'FaceColor', 'Blue');
    set(handles.plot_Signal_Areas(i), 'FaceColor', 'Green');
    set(handles.text_StartTime, 'String', ['Start time (s): ',...
        num2str(handles.a(i).startTime)]);
    set(handles.text_Duration, 'String', ['Duration (s): ',...
        num2str(handles.a(i).duration)]);
    set(handles.text_Window, 'String', ['Window (s): ',...
        num2str(handles.a(i).window)]);
    
    % FFT of the selected area
    iFT = fft(sum(iData(:, [iEEGChannels; iBipolarChannels]), 2));
    iFT = iFT(1:floor(N/2));
    iFT(1:2) = 0; % 1:2 are error prone elements and unnecessary in this case.
    [maxValue, maxIndex] = max(abs(iFT));
    iFT = abs(iFT)./maxValue;
    frequencyAxis = handles.frequency*(0:N/2-1)./N;
    
    % Plot FFT data
    plot(frequencyAxis, iFT, 'Parent', handles.axes_Frequency);
    text(frequencyAxis(maxIndex), 0.9, ['\leftarrow ',...
        num2str(frequencyAxis(maxIndex)), ' Hz'],...
        'Parent', handles.axes_Frequency);
    xlabel(handles.axes_Frequency, 'Frequency (Hz)');
    %freqLim = frequencyAxis(maxIndex)+[-15, +15];
    %if freqLim(1) < 0
    %    freqLim(1) = 0;
    %end
    %if freqLim(2) > frequencyAxis(end)
    %    freqLim(2) = frequencyAxis(end);
    %end
    %set(handles.axes_Frequency, 'XLim', freqLim);
    
    % Export VPP data to excel
    
    set(excelSheetRange, 'Value', [frequencyAxis(maxIndex),...
        max(max(VPP(:, iEEGChannels))),...
        min(min(VPP(:, iEEGChannels))),...
        max(max(VPP(:, iBipolarChannels))),...
        min(min(VPP(:, iBipolarChannels)))]);
    
    % Plot VPP data
    plot(iEEGChannels, max(VPP(:, iEEGChannels)), 'b^',...
        'MarkerSize', 7, 'MarkerFaceColor', 'b',...
        'Parent', handles.axes_Channels);
    hold(handles.axes_Channels, 'on');
    plot(iEEGChannels, min(VPP(:, iEEGChannels)), 'bv',...
        'MarkerSize', 7, 'MarkerFaceColor', 'none',...
        'Parent', handles.axes_Channels);
    plot(iBipolarChannels, max(VPP(:, iBipolarChannels)), 'r^',...
        'MarkerSize', 7, 'MarkerFaceColor', 'r',...
        'Parent', handles.axes_Channels);
    plot(iBipolarChannels, min(VPP(:, iBipolarChannels)), 'rv',...
        'MarkerSize', 7, 'MarkerFaceColor', 'none',...
        'Parent', handles.axes_Channels);
    
    line([0, handles.nChannels+1],...
        mean(max(VPP(:, iEEGChannels))).*[1, 1],...
        'Color', 'Blue', 'Parent', handles.axes_Channels);
    line([0, handles.nChannels+1],...
        mean(min(VPP(:, iEEGChannels))).*[1, 1],...
        'Color', 'Blue', 'Parent', handles.axes_Channels,...
        'LineStyle', '--');
    line([0, handles.nChannels+1],...
        mean(max(VPP(:, iBipolarChannels))).*[1, 1],...
        'Color', 'Red', 'Parent', handles.axes_Channels);
    line([0, handles.nChannels+1],...
        mean(min(VPP(:, iBipolarChannels))).*[1, 1],...
        'Color', 'Red', 'Parent', handles.axes_Channels,...
        'LineStyle', '--');
    handles.plot_activeChannel =...
        plot(handles.axes_Channels, -5, -5, 'go',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    hold(handles.axes_Channels, 'off');
    set(handles.axes_Channels, 'XLim', [0, handles.nChannels+1],...
        'XMinorTick', 'on', 'XGrid', 'on');
    ylabel(handles.axes_Channels, 'Amplitude (µV)');
    xlabel(handles.axes_Channels, 'Channel number');
    legend(handles.axes_Channels,...
        'EEG max', 'EEG min', 'Bipolar max', 'Bipolar min',...
        'Mean EEG max', 'Mean EEG min', 'Mean Bipolar Max',...
        'Mean Bipolar Min', 'Location', 'NorthEastOutside');
    
    % Paste Channels axes to Excel 300 DPI
    
    excelColumns = strcat('G', num2str((i*23-21)));
    excelSheetRange = get(excelSheet, 'Range', excelColumns);
    excelSheetRange.Select;
    
    tmpAxes = copyobj(handles.axes_Channels, tmpFigure);
    legend(tmpAxes,...
        'EEG max', 'EEG min', 'Bipolar max', 'Bipolar min',...
        'Mean EEG max', 'Mean EEG min', 'Mean Bipolar Max',...
        'Mean Bipolar Min', 'Location', 'NorthEastOutside');
    set(tmpAxes, 'units', 'normalized', 'position',...
        [0.05 0.11 0.77 0.78]);
    title(tmpAxes, {['Alueen: ', num2str(i),...
        ' kanavakohtainen VPP (µV) ~taajuudella ', ...
        num2str(frequencyAxis(maxIndex)), ' Hz'],...
        ['Start time: ', num2str(handles.a(i).startTime), ' s',...
        ', Duration: ', num2str(handles.a(i).duration), ' s',...
        ', Window: ', num2str(handles.a(i).window), ' s']},...
        'FontSize', 12, 'FontWeight', 'bold');
    
    print('-dmeta', tmpFigure, '-r300');
    invoke(excelApplication.Selection, 'PasteSpecial');
    
    delete(tmpAxes);
    
end
waitbar(0.8,prog)

excelSheetRange = get(excelSheet, 'Range', 'A1');
excelSheetRange.Select;
set(excelSheetRange, 'Value', ['Raportti tiedostosta: ',...
    handles.filePath, handles.fileName]);

if ~exist(outPutFull, 'file')
    invoke(excelWorkbook, 'SaveAs', outPutFull, 1);
else
    invoke(excelWorkbook, 'Save');
end
invoke(excelWorkbook, 'Close');
invoke(excelApplication, 'Quit');
delete(excelApplication);
delete(tmpFigure);
waitbar(1,prog)

close(prog);
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);


set(handles.axes_Signal,'UserData',[])


% --- Executes on button press in pushReset.
function pushReset_Callback(hObject, eventdata, handles)
plotUpdate(hObject, handles, [handles.time(1), handles.time(end)]);

function out = nansum(in,dim)

if isempty(in)
    out = [];
    return;
end

if nargin == 1
    tmp = size(dim);
    [~,dim] = min(tmp(tmp>1));
    if isempty(dim)
        dim = 1;
    end
end

in(isnan(in)) = 0;

out = sum(in,dim);

