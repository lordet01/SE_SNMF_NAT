function run_IMCRA(DB_input, DB_output, p)

addpath('src')


disp('------IMCRA-NMF enhancement------');
InFileList = dir(DB_input);
for j=3:p.filegap:length(InFileList)
    %Load multi-channel I/O
    path_in = 0;
    path_denoise = 0;
    for i = 1:p.ch
        path_in_tmp = [DB_input,'/',InFileList(j+i-1).name];
        file_name = InFileList(j+i-1).name;
        path_denoise_tmp = [DB_output,'/',file_name];
        if i == 1
            path_in = path_in_tmp;
            path_denoise = path_denoise_tmp;
        else
            path_in = [path_in; path_in_tmp];
            path_denoise = [path_denoise; path_denoise_tmp];
        end
    end
    disp(['Enhancing: ',file_name]);
    f_tmp = fopen(path_denoise(p.ch,:));
    if f_tmp == -1 || p.ForceRewrite
        proc_IMCRA(path_in, path_denoise);
    else
        fclose(f_tmp);
    end
end
%--------------------------------------------------------------------------

fclose('all');

