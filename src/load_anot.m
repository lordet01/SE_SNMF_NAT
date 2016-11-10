function [v_start, v_end] = load_anot(filename, len, p)

ANNO_PATH = ['training_anno/',filename,'_sid.txt'];

if 1 + fopen(ANNO_PATH,'rb') %#ok<BDLOG> %Check whether annotation file exists
    ANNO_SAMPLE = load(ANNO_PATH);
    ANNO_SAMPLE = ceil(ANNO_SAMPLE.*p.fs);

    if ANNO_SAMPLE(1) == 0
        ANNO_SAMPLE(1) = 1;
    end

    if ANNO_SAMPLE(2) > len
        ANNO_SAMPLE(2) = len;
    end
    
    v_start = ANNO_SAMPLE(1); v_end = ANNO_SAMPLE(2);
    fclose('all');
end





end