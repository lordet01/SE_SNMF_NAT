function x=wavread_rawM(path)
fid = fopen(path, 'rb');
x = fread(fid,inf,'int16');
fclose(fid);
end