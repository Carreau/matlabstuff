function ret = pp(s)
    format = '';
    [pathstr, name, ext] = fileparts(s.rawdata.name);
    [pathstr, sub , ~  ] = fileparts(pathstr);
    A=[sprintf('%s/%s%s',sub,name,ext)];
    format = strcat(format,'%s\t');

    A2 = s.rawdata.datevec(1:5)
    format = strcat(format,'%04i-%02i-%02i %02i:%02i\t');
    
    A3 = s.approachAngleDeg;
    format = strcat(format,'%4.0f\t');

    A4 = s.param.speed.value;
    format = strcat(format,'%3i\t');
     
    A4 = s.param.resting_time.value;
    format = strcat(format,'%3i\t');
    
    ret = sprintf(format,A,A2,A3,A4);

end
