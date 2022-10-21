clear, clc;
fid = fopen('text.txt');
tline = fgetl(fid);
enc = 0;
CA =[];M = [];
lc = 0;
while ischar(tline)
    
    % get the second quantity
    newStr = split(tline);
    
   
   CA =  [CA,  str2double(string(newStr(3)))];
    

    tline = fgetl(fid);
    enc = enc+1;
    
    
    lc = lc+1;
end

CA = [CA(1), CA(6:10), CA(2:5)];

part1 = "%f +%f dx_1 + %f dx_2 + %f dx_3 + %f dx_4  \\ + %f dx_5 +%f  dx_6  +%f  dx_7  + %f dx_8  + %f dx_9  \n ";

fprintf(part1, CA);


    
     

    




