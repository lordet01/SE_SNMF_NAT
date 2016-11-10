function wavwrite_rawM(path, x)
x = cast(x, 'int16');
fid = fopen(path, 'wb');
fwrite(fid,x,'int16');
fclose(fid);
end