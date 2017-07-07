function [K,R,t,im_size]=readCalibrationOrientation_Strecha(image_path,image_name)

filename=strcat(image_path,image_name,'.camera');
calibration_file = fopen(filename,'r');

K=[str2num(fgetl(calibration_file));...
    str2num(fgetl(calibration_file));...
    str2num(fgetl(calibration_file))];

fgetl(calibration_file);

R=[str2num(fgetl(calibration_file));...
    str2num(fgetl(calibration_file));...
    str2num(fgetl(calibration_file))].';

t=-R*[str2num(fgetl(calibration_file))].';

im_size=[str2num(fgetl(calibration_file))].';

fclose(calibration_file);

end


