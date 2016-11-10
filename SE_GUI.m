function varargout = SE_GUI(varargin)
%SE_GUI M-file for SE_GUI.fig
%      SE_GUI, by itself, creates a new SE_GUI or raises the existing
%      singleton*.
%
%      H = SE_GUI returns the handle to a new SE_GUI or the handle to
%      the existing singleton*.
%
%      SE_GUI('Property','Value',...) creates a new SE_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SE_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SE_GUI('CALLBACK') and SE_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SE_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SE_GUI

% Last Modified by GUIDE v2.5 03-Sep-2014 11:06:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SE_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SE_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before SE_GUI is made visible.
function SE_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SE_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% %% ==========Image button initialization================
% % Set the icons of the buttons. These functions are defined later in the
% % M-code. The 'playbutton' and the 'stopbutton' functions read in image
% % files and then extracts enough data from them to form a small icons. The
% % 'fasterbutton' and 'slowerbutton' function do not do anything at this
% % point.
% set(handles.play,'cdata',playbutton);
% set(handles.stop,'cdata',stopbutton);
% set(handles.faster,'cdata',fasterbutton);
% set(handles.slower,'cdata',slowerbutton);
% 
% % Set graphics parameters.
% set(handles.figure1,'color','k')


% UIWAIT makes SE_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%% ==================Begin Initialization===============
% addpath('src/OD_NMF');
addpath('src');
addpath('settings');
addpath('src/func');
addpath('src/MS');
mkdir('wav');

global prmUI;

prmUI.mode_MS = 0;
prmUI.mode_BNMF = 0;
prmUI.mode_BNMFNA = 1;


%%On-the-fly recording related initialization
prmUI.record_flag = 0; %flag == 0; stop_state, flag == 1; recording_state
prmUI.PTT = 1;

%%Speech and Noise Basis Plotting
%a-priori basis load
% load(['src/BNMF/basis/noise/','3s_3rank_1024fft','/noise_','h','_KL','50_recent','.mat']);
% load(['src/BNMF/basis/speech/','0.3s_6rank_1024fft','/speech_(21-25,31-35)_KL','50','.mat']);

% axes(handles.axes_SpeechBasis);
% contour(B_s);
% xlabel('Basis Rank');
% ylabel('Frequency Bin');
% 
% axes(handles.axes_NoiseBasis);
% contour(B_n);
% xlabel('Basis Rank');
% ylabel('Frequency Bin');

%tgrabaudio('stop',16000, 0.1);

% prmUI.B_n_update = B_n; %Initial B_n

%Initial status in text
set(handles.txt_status, 'String','Mic Idle');
set(handles.txt_status, 'ForegroundColor',[1.0, 0, 0]);

%Full plot open
%file I/O initialization
path_in = 'wav/in.RAW';
path_denoised = 'wav/out.RAW';
path_noise = 'wav/e_noise.RAW';
path_speech = 'wav/e_speech.RAW';

fin = fopen(path_in,'rb');
fout = fopen(path_denoised,'rb');
fnoise = fopen(path_noise,'rb');
fspeech = fopen(path_speech,'rb');

if fin > -1 && fout > -1
    s_f_in = fread(fin, inf, 'int16');
    s_f_out = fread(fout, inf, 'int16');
%     s_f_en = fread(fnoise, inf, 'int16');
%     s_f_es = fread(fspeech, inf, 'int16');
else
    s_f_in = 0; s_f_out = 0; s_f_en = 0; s_f_es = 0;
end

prmUI.s_f_in = s_f_in;
prmUI.s_f_out = s_f_out;
% prmUI.s_f_en = s_f_en;
% prmUI.s_f_es = s_f_es;

%% Plot all file into spectrogram
axes(handles.axes_InputSignal);
myspectrogram(s_f_in, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
xlabel('Frame Index');
ylabel('Frequency (Hz)');

% axes(handles.axes_EstimSpeech);
% myspectrogram(s_f_es, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
% axes(handles.axes_EstimNoise);
% myspectrogram(s_f_en, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
axes(handles.axes_EnhanceSpeech);
myspectrogram(s_f_out, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');

% --- Outputs from this function are returned to the command line.
function varargout = SE_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function SE_GUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SE_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in bttn_Play2.
function bttn_Play2_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_Play2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prmUI;
persistent p;

state=get(handles.bttn_Play2,'Value'); 
tmp = prmUI.s_f_es ./32767;
p = audioplayer(tmp,16000);
play(p);
if state == 1
    set(handles.bttn_Play2,'String','||'); 
    resume(p);
else
    set(handles.bttn_Play2,'String','Play'); 
    pause(p);
end

% --- Executes on button press in bttn_Play1.
function bttn_Play1_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_Play1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prmUI;
persistent p;

state=get(handles.bttn_Play1,'Value'); 
tmp = prmUI.s_f_in ./32767;
p = audioplayer(tmp,16000);
play(p);
if state == 1
    set(handles.bttn_Play1,'String','||'); 
    resume(p);
else
    set(handles.bttn_Play1,'String','Play'); 
    pause(p);
end





% --- Executes on button press in bttn_RecordSpeech.
function bttn_RecordSpeech_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_RecordSpeech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bttn_RecordSpeech
global prmUI;

if prmUI.record_flag == 1
    if prmUI.PTT == 0;
        prmUI.PTT = 1;
        set(handles.bttn_RecordSpeech,'String','Push to Finish');
        set(handles.bttn_RecordSpeech, 'ForegroundColor', [1, 0, 0]);
    else
        prmUI.PTT = 0;
        set(handles.bttn_RecordSpeech,'String','Push to Speak');
        set(handles.bttn_RecordSpeech, 'ForegroundColor', [0, 0, 1]);
    end
end

% --- Executes on button press in bttn_Play3.
function bttn_Play3_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_Play3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prmUI;
persistent p;

state=get(handles.bttn_Play3,'Value'); 
tmp = prmUI.s_f_en ./32767;
p = audioplayer(tmp,16000);
play(p);
if state == 1
    set(handles.bttn_Play3,'String','||'); 
    resume(p);
else
    set(handles.bttn_Play3,'String','Play'); 
    pause(p);
end

% --- Executes on button press in bttn_Play4.
function bttn_Play4_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_Play4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prmUI;
persistent p;

state=get(handles.bttn_Play4,'Value'); 
tmp = prmUI.s_f_out ./32767;
p = audioplayer(tmp,16000);
play(p);
if state == 1
    set(handles.bttn_Play4,'String','||'); 
    resume(p);
else
    set(handles.bttn_Play4,'String','Play'); 
    pause(p);
end


% --- Executes on button press in bttn_OnOff.
function bttn_OnOff_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_OnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bttn_OnOff
global prmUI;
persistent obj;

%a-priori basis load
% if prmUI.mode_BNMFNA
%     load(['src/BNMF/basis/noise/','3s_3rank_1024fft','/noise_','h','_KL','50_recent','.mat']);
% else
%     load(['src/BNMF/basis/noise/','3s_3rank_1024fft','/noise_','h','_KL','50_fix','.mat']);
% end

%MS initialization
init_MS;

%BNMF initialization
initial_setting_Proposed_Techwin_201603_RT;

%Load Event basis
load(['basis/Clean_train_car_alarm/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = B_DFT_sub; B_Mel_x = B_Mel_sub;
load(['basis/Clean_train_scream/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];
load(['basis/Clean_train_TIMIT_test/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];
B_Mel_x = B_DFT_x;

%Load Noise basis
load(['basis/CHiME3_bgn_ch6/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_50.mat']);
B_DFT_d = B_DFT_sub;  B_Mel_d = B_Mel_sub;
B_Mel_d = B_DFT_d;

frame_len = p.framelength;
frame_shift = p.frameshift;
shift_ratio = frame_len / frame_shift; %should be integer!!

%Buffer initialization
g=init_buff(B_Mel_x, B_Mel_d, B_DFT_x, B_DFT_d, p);
s_in = zeros(frame_len, 1);
s_out = zeros(frame_len, 1);
% s_noise = zeros(frame_len, 1);
% s_speech = zeros(frame_len, 1);


rec_cycle = 1;
len = frame_shift * rec_cycle; ch = 1; fs=16000; buff_frm=10; interval=10; interval_basis = 100; 
RT_scale=get(handles.edit_RTscale, 'String');
RT_scale = str2double(RT_scale);
pRT_scale = RT_scale; mRT_scale = -1 * RT_scale;
if prmUI.record_flag == 0
% %     record(recObj);
% %     set(handles.text_RecStatus, 'String', 'Recording...');
% %     set(handles.text_RecStatus, 'ForegroundColor',[1,0,0]);
    
    prmUI.record_flag = 1;
    set(handles.bttn_OnOff,'String','Record Off'); 
    set(handles.bttn_OnOff, 'ForegroundColor',[0, 0, 0]);
    
    %file I/O initialization
    path_in = 'wav/in.RAW';
    path_denoised = 'wav/out.RAW';
%     path_noise = 'wav/e_noise.RAW';
%     path_speech = 'wav/e_speech.RAW';

    fin = fopen(path_in,'wb');
    fout = fopen(path_denoised,'wb');
%     fnoise = fopen(path_noise,'wb');
%     fspeech = fopen(path_speech,'wb');    

    %Initialize Audio recording
    [obj]=dsp_record(0, fs, ch, len, 'start');
    
    frame_Input_buff = zeros(buff_frm*frame_shift,1);
    update_cnt_frm = 1;
    update_cnt_basis = 1;
    init_cnt = 1;
    total_cnt = 1;
    cnt_initMS = 0; 
    
    %Plot Initialization
    axes(handles.axes_InputSignal);
    axes(handles.axes_NoiseBasis);
    
    while prmUI.record_flag == 1        
        %s_in_sub=tgrabaudio(frame_shift);
        [obj, s_in_rec]=dsp_record(obj, fs, ch, len, 'record');
        
        s_in_rec = s_in_rec(:,1) .* 32767;
        rec_cycle_cnt = 1;
        while rec_cycle_cnt <= rec_cycle
                %%do realtime processing here
            s_in_sub = s_in_rec(1+(rec_cycle_cnt-1)*frame_shift : rec_cycle_cnt*frame_shift);

            %Frame_wise queing
            s_in(1:frame_len-frame_shift) = s_in(frame_shift+1:frame_len');
            s_in(frame_len-frame_shift+1:frame_len) = s_in_sub;

            if total_cnt >= shift_ratio
                %Put frame-wise algorithm here
                if prmUI.PTT == 0 && prmUI.mode_BNMFNA == 1
                    %update status text
    %                 set(handles.txt_status, 'String','Adapting Noise...');
    %                 set(handles.txt_status, 'ForegroundColor',[0, 0, 1.0]);
    %                 [g] = bnmf_adapt_RT(s_in, p.blk_len_adapt, g, p);
                elseif prmUI.PTT == 1 && prmUI.mode_MS == 0 %BNMF enhancement
                    if prmUI.mode_BNMF == 1 %Conventional NMF (basis fixed)
                        initial_setting_SNMF_Techwin_201603_RT;
                    else
                        initial_setting_Proposed_Techwin_201603_RT; %Proposed
                    end

                    fwrite(fin, s_in_sub, 'int16');
                    [~, ~, s_d_frame, g] = bnmf_sep_event_RT_IS16(s_in', total_cnt, g, p);

                    s_out(1:frame_len-frame_shift) = s_out(frame_shift+1:frame_len);
                    s_out(frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
                    s_out = s_out + s_d_frame';
                    s_out_sub = s_out(1:frame_shift);
                    fwrite(fout, s_out_sub, 'int16');

    %                 s_speech(1:frame_len-frame_shift) = s_speech(frame_shift+1:frame_len);
    %                 s_speech(frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
    %                 s_speech = s_speech + s_est_frame;
    %                 s_speech_sub = s_speech(1:frame_shift);
    %                 fwrite(fspeech, s_speech_sub, 'int16');

    %                 s_noise(1:frame_len-frame_shift) = s_noise(frame_shift+1:frame_len);
    %                 s_noise(frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
    %                 s_noise = s_noise + n_est_frame;
    %                 s_noise_sub = s_noise(1:frame_shift);
    %                 fwrite(fnoise, s_noise_sub, 'int16');
                elseif prmUI.PTT == 1 && prmUI.mode_MS == 1 %MS enhancement
                    if cnt_initMS == 0 %initial call of MS enhancement
                        [s_out_sub, z]=ssubmmse(s_in_sub,fs);
                        cnt_initMS = 1;
                    else
                        [s_out_sub, z]=ssubmmse(s_in_sub,z);
                    end

                    fwrite(fin, s_in_sub, 'int16');
                    fwrite(fout, s_out_sub, 'int16');
                    %blank for sepech, noise estimation buffer that doesnt
                    %exist in MS mode
                    blank_sub = zeros(size(s_in_sub));
%                     fwrite(fspeech, blank_sub, 'int16');
%                     fwrite(fnoise, blank_sub, 'int16');
                end

               %% Buffer for real-time plotting
                if init_cnt <= buff_frm
                    frame_Input_buff((init_cnt-1)*frame_shift+1:(init_cnt)*frame_shift, :)=s_in_sub;
                    init_cnt = init_cnt + 1;
                else
                    frame_Input_buff(1:(buff_frm-1)*frame_shift, :) = frame_Input_buff(frame_shift+1:buff_frm*frame_shift, :);
                    frame_Input_buff((buff_frm-1)*frame_shift+1:(buff_frm)*frame_shift, :)=s_in_sub;
                end

                if update_cnt_frm == interval
                    plot(handles.axes_InputSignal, frame_Input_buff(:,1));
                    axis(handles.axes_InputSignal, [1 frame_shift*buff_frm mRT_scale pRT_scale]);
                    xlabel(handles.axes_InputSignal, 'Time Sample');
                    ylabel(handles.axes_InputSignal, 'Sample Value (16bit)');
                    drawnow;

                    %bug prevention codes
                    if prmUI.record_flag  == 0 && prmUI.PTT == 0
                        path_in = 'wav/in.RAW';
                        fin = fopen(path_in,'rb');
                        s_f_in = fread(fin, inf, 'int16');
                        axes(handles.axes_InputSignal);
                        myspectrogram(s_f_in, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
                        xlabel(handles.axes_InputSignal, 'Frame Index');
                        ylabel(handles.axes_InputSignal, 'Frequency (Hz)');
                    end
                    update_cnt_frm = 1;
                end

                if update_cnt_basis == interval_basis
                    %%Noise Basis Update Plotting
                    if prmUI.mode_BNMFNA
    %                     prmUI.B_n_update = g.B_Mel_d;
                        B_n = g.B_DFT_d;
    %                     contour(handles.axes_NoiseBasis, B_n);
                        xlabel(handles.axes_NoiseBasis, 'Basis Rank');
                        ylabel(handles.axes_NoiseBasis, 'Frequency Bin');
                        drawnow;

    %                     save('src/BNMF/basis/noise/3s_3rank_1024fft/noise_h_KL50_recent.mat', 'B_n');
                    end
                    update_cnt_basis =  1;
                end

                update_cnt_frm = update_cnt_frm + 1;
                update_cnt_basis = update_cnt_basis + 1;
            end
            total_cnt = total_cnt + 1;

            %Status Update
            if prmUI.record_flag == 1
                if prmUI.PTT == 1;
                    if prmUI.mode_MS
                        %update status text
                        set(handles.txt_status, 'String','(MS) Enhancing Speech...');
                        set(handles.txt_status, 'ForegroundColor',[0, 0, 1]);
                    elseif prmUI.mode_BNMF
                        set(handles.txt_status, 'String','(BNMF) Enhancing Speech...');
                        set(handles.txt_status, 'ForegroundColor',[0, 0, 1]);
                    elseif prmUI. mode_BNMFNA
                        set(handles.txt_status, 'String','(Proposed) Enhancing Speech...');
                        set(handles.txt_status, 'ForegroundColor',[0, 0, 1]);
                    end
                else
                    if prmUI.mode_MS
                        %update status text
                        set(handles.txt_status, 'String','(MS) Adpating Noise...');
                        set(handles.txt_status, 'ForegroundColor',[1, 0, 0]);
                    elseif prmUI.mode_BNMF
                        set(handles.txt_status, 'String','(BNMF) Idle...');
                        set(handles.txt_status, 'ForegroundColor',[1, 0, 0]);
                    elseif prmUI. mode_BNMFNA
                        set(handles.txt_status, 'String','(Proposed) Adpating Noise...');
                        set(handles.txt_status, 'ForegroundColor',[1, 0, 0]);
                    end  
                end  
            end
            rec_cycle_cnt = rec_cycle_cnt + 1;
        end
    end
else
    prmUI.record_flag  = 0;
    prmUI.PTT = 1;
    
    dsp_record(obj, fs, ch, len, 'start');
    
    set(handles.bttn_OnOff,'String','Record On'); 
    set(handles.bttn_OnOff, 'ForegroundColor',[1, 0, 0]);    

    %update status text
    set(handles.txt_status, 'String','Mic Idle');
    set(handles.txt_status, 'ForegroundColor',[1.0, 0, 0]);    
    
    %dsp_record(obj, fs, ch, len, 'stop');
    fclose('all');
    
    %Full plot open
    %file I/O initialization
    path_in = 'wav/in.RAW';
    path_denoised = 'wav/out.RAW';
    path_noise = 'wav/e_noise.RAW';
    path_speech = 'wav/e_speech.RAW';

    fin = fopen(path_in,'rb');
    fout = fopen(path_denoised,'rb');
    fnoise = fopen(path_noise,'rb');
    fspeech = fopen(path_speech,'rb');    
    
    s_f_in = fread(fin, inf, 'int16');
    s_f_out = fread(fout, inf, 'int16');
    s_f_en = fread(fnoise, inf, 'int16');
    s_f_es = fread(fspeech, inf, 'int16');
    
    prmUI.s_f_in = s_f_in;
    prmUI.s_f_out = s_f_out;
    prmUI.s_f_en = s_f_en;
    prmUI.s_f_es = s_f_es;
    
    %% Plot all file into spectrogram
    axes(handles.axes_InputSignal);
    myspectrogram(s_f_in, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
    xlabel('Frame Index');
    ylabel('Frequency (Hz)');
    
%     axes(handles.axes_EstimSpeech);
%     myspectrogram(s_f_es, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
%     axes(handles.axes_EstimNoise);
%     myspectrogram(s_f_en, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
    axes(handles.axes_EnhanceSpeech);
    myspectrogram(s_f_out, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');    
end


% --- Executes on button press in radio_MS.
function radio_MS_Callback(hObject, eventdata, handles)
% hObject    handle to radio_MS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_MS
global prmUI;
prmUI.mode_MS = 1;
prmUI.mode_BNMF = 0;
prmUI.mode_BNMFNA = 0;

axes(handles.axes_SpeechBasis);
plot(0);
axis off;

axes(handles.axes_NoiseBasis);
plot(0);
axis off;

% --- Executes on button press in radio_BNMF.
function radio_BNMF_Callback(hObject, eventdata, handles)
% hObject    handle to radio_BNMF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_BNMF
global prmUI;
prmUI.mode_MS = 0;
prmUI.mode_BNMF = 1;
prmUI.mode_BNMFNA = 0;

%Load Event basis
load(['basis/Clean_train_car_alarm/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = B_DFT_sub; B_Mel_x = B_Mel_sub;
load(['basis/Clean_train_scream/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];
load(['basis/Clean_train_TIMIT_test/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];

%Load Noise basis
load(['basis/CHiME3_bgn_ch6/','IS16_Splice0-SNMF-preemph_p2_DD1','/R_50.mat']);
B_DFT_d = B_DFT_sub; B_Mel_d = B_Mel_sub;

axes(handles.axes_SpeechBasis);
contour(B_DFT_x);
xlabel('Basis Rank');
ylabel('Frequency Bin');

% load(['src/BNMF/basis/noise/','3s_3rank_1024fft','/noise_','h','_KL','50_fix','.mat']);
axes(handles.axes_NoiseBasis);
contour(B_DFT_d);
xlabel('Basis Rank');
ylabel('Frequency Bin');


% --- Executes on button press in radio_BNMFNA.
function radio_BNMFNA_Callback(hObject, eventdata, handles)
% hObject    handle to radio_BNMFNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_BNMFNA
global prmUI;
prmUI.mode_MS = 0;
prmUI.mode_BNMF = 0;
prmUI.mode_BNMFNA = 1;

%Load Event basis
load(['basis/Clean_train_car_alarm/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = B_DFT_sub; B_Mel_x = B_Mel_sub;
load(['basis/Clean_train_scream/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];
load(['basis/Clean_train_TIMIT_test/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);
B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];

%Load Noise basis
load(['basis/CHiME3_bgn_ch6/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_50.mat']);
B_DFT_d = B_DFT_sub; B_Mel_d = B_Mel_sub;

% load(['src/BNMF/basis/speech/','0.3s_6rank_1024fft','/speech_(21-25,31-35)_KL','50','.mat']);
axes(handles.axes_SpeechBasis);
contour(B_DFT_x);
xlabel('Basis Rank');
ylabel('Frequency Bin');

% load(['src/BNMF/basis/noise/','3s_3rank_1024fft','/noise_','h','_KL','50_recent','.mat']);
axes(handles.axes_NoiseBasis);
contour(B_DFT_d);
xlabel('Basis Rank');
ylabel('Frequency Bin');


% --- Executes on button press in bttn_SaveB_n.
function bttn_SaveB_n_Callback(hObject, eventdata, handles)
% hObject    handle to bttn_SaveB_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global prmUI;
% B_n = prmUI.B_n_update;
% save('src/BNMF/basis/noise/3s_3rank_1024fft/noise_h_KL50_recent.mat', 'B_n');


% --- Executes on button press in radio_Spec.
function radio_Spec_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Spec
global prmUI;
if prmUI.record_flag == 0
    %% Plot all file into spectrogram
    axes(handles.axes_InputSignal);
    myspectrogram(prmUI.s_f_in, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
    xlabel('Frame Index');
    ylabel('Frequency (Hz)');
    
    axes(handles.axes_EstimSpeech);
    myspectrogram(prmUI.s_f_es, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
    axes(handles.axes_EstimNoise);
    myspectrogram(prmUI.s_f_en, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');
    axes(handles.axes_EnhanceSpeech);
    myspectrogram(prmUI.s_f_out, 16000, [32, 8], @hanning, 512, [-99 0], false, 'default');    
end


% --- Executes on button press in radio_Wave.
function radio_Wave_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Wave
global prmUI;
if prmUI.record_flag == 0
    %% Plot all file into spectrogram
    len = length(prmUI.s_f_in);
    axes(handles.axes_InputSignal);
    plot(prmUI.s_f_in);
    maxval = 1.1*max(abs(prmUI.s_f_in));
    axis([1 len -maxval maxval]);
    xlabel('Time Sample');
    ylabel('Sample Value (16bit)');
    
    axes(handles.axes_EstimSpeech);
    plot(prmUI.s_f_es);
    axis([1 len -maxval maxval]);
    
    axes(handles.axes_EstimNoise);
    plot(prmUI.s_f_en);
    axis([1 len -maxval maxval]);
    
    axes(handles.axes_EnhanceSpeech);
    plot(prmUI.s_f_out);
    axis([1 len -maxval maxval]);
end



function edit_RTscale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RTscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RTscale as text
%        str2double(get(hObject,'String')) returns contents of edit_RTscale as a double


% --- Executes during object creation, after setting all properties.
function edit_RTscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RTscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
