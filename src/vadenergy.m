function [vad,v_start_final,v_end_final]=vadenergy(x_in, p)

fs = p.fs;
bg_len = p.bg_len;
min_voiced_len = p.min_voiced_len;
min_unvoiced_len = p.min_unvoiced_len;
thr = p.thr;

%% BG Energy calculation
% if bg_mean == 0
    bg_mean = mean(abs(x_in(1:bg_len)));
% end

%% frame by frame VAD decision

x_in_abs = abs(x_in); 
vad = zeros(size(x_in));

frame_len = 0.02 *fs; % 20ms
frame_shift = frame_len / 2;
frame_num = floor(length(x_in) / frame_shift);

i = 1;
for k = 1:frame_num-1
    frame_crnt = x_in_abs(i:i+frame_len-1);
    
   if ( ((mean(frame_crnt) - bg_mean) / mean(frame_crnt)) > thr )
       vad(i:i+frame_len-1) = ones(size(vad(i:i+frame_len-1)));
   end
   
   i = i + frame_shift;
end


%% VAD smoothing 
v_state = 0;
v_cnt = 0;
v_start = 1;
v_end = 1;

uv_cnt = 0;
uv_start = 1;
uv_end = 1;

v_start_final = 1;
v_end_final = length(x_in);

update_flag = 1;
for i = 2:length(vad)
    %% (impulsive unvoiced compensation)
    if vad(i-1) == 1 && vad(i) == 0
        v_state = 0;
        
        uv_start = i;
    end
    
    if v_state == 0;
        uv_cnt = uv_cnt + 1;
    end
    
    if vad(i-1) == 0 && vad(i) == 1
        v_state = 1;
        uv_end = i-1;
        if uv_cnt < min_unvoiced_len
            vad(uv_start:uv_end) = ones(size(vad(uv_start:uv_end)));
            
            %v_cnt = v_cnt_prior + length(vad(uv_start:uv_end));
        end
        
        %uv_cnt_prior = uv_cnt;
        uv_cnt = 0;
    end    
end

for i = 2:length(vad)
    %% (impulsive voiced compensation)
    if vad(i-1) == 0 && vad(i) == 1
        v_state = 1;
        
        v_start = i;
        
        %v_start_final update
        if update_flag == 1
            v_start_final = v_start;
        end
    end
    
    if v_state == 1;
        v_cnt = v_cnt + 1;
    end
    
    if vad(i-1) == 1 && vad(i) == 0
        v_state = 0;
        v_end = i-1;
        
        %v_end_final update
        if update_flag == 1
            v_end_final = v_end;
            update_flag = 0;
        end
        
        if v_cnt < min_voiced_len
            vad(v_start:v_end) = zeros(size(vad(v_start:v_end)));
            
            update_flag = 1;
            
            %uv_cnt = uv_cnt_prior + length(vad(v_start:v_end));
        end
        
        %v_cnt_prior = v_cnt;
        v_cnt = 0;
    end
end


end


