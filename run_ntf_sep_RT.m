function run_ntf_sep_RT(DB_input, DB_output, B1_x, B1_d, B2_x, B2_d, p)

addpath('src')
mkdir([DB_output,'/event/']);
mkdir([DB_output,'/noise/']);

%% -----------------BNTF separation------------------
disp('------Blockwise NMF(NTF) separation------');
InFileList = dir(DB_input);
for j=3:p.filegap:length(InFileList)
    %Load multi-channel I/O
    path_in = 0;
    path_denoise = 0;
    path_event = 0;
    path_noise = 0;
    for i = 1:p.ch
        path_in_tmp = [DB_input,'/',InFileList(j+i-1).name];
        file_name = InFileList(j+i-1).name;
        path_denoise_tmp = [DB_output,'/',file_name];
        path_event_tmp = [DB_output,'/event/',file_name];
        path_noise_tmp = [DB_output,'/noise/',file_name];
        if i == 1
            path_in = path_in_tmp;
            path_denoise = path_denoise_tmp;
            path_event = path_event_tmp;
            path_noise = path_noise_tmp;
        else
            path_in = [path_in; path_in_tmp];
            path_denoise = [path_denoise; path_denoise_tmp];
            path_event = [path_event;path_event_tmp];
            path_noise = [path_noise;path_noise_tmp];
        end
    end
    disp(['Enhancing: ',file_name]);
    f_tmp = fopen(path_denoise(p.ch,:));
    if f_tmp == -1 || p.ForceRewrite
        NTF_sep_event_RT(path_in, path_event, path_noise, path_denoise, B1_x, B1_d, B2_x, B2_d, p);
    else
        fclose(f_tmp);
    end
end
%--------------------------------------------------------------------------

fclose('all');

