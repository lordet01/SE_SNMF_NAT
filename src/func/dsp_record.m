function [obj, frame]=dsp_record(obj, fs, ch, len, stat)

if strcmp(stat, 'start')
    %Initialization of recording obj.
    obj = dsp.AudioRecorder('SampleRate', fs, ...
        'NumChannels', ch, ...
        'SamplesPerFrame', len, ...
        'OutputDataType', 'double', ...
        'QueueDuration', 1);
    disp('Microphone recording initialized...');
    
    frame = 0;
end

if strcmp(stat, 'stop')
    release(obj);
    disp('Microphone recording stopped...');
    frame = 0;
end

if strcmp(stat, 'record')
    frame = step(obj);
end