clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

for beam = 1:128
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5 + depth;
    time_to_point(ii) = d(ii)/speed;
end

delay_matrix(jj,:,beam) = time_to_point; %delays

end

for aa = 1:128
    delayed_channel(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
end


end


for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_channel)));
max_data = max(max(max(delayed_channel)));
imagesc(delayed_channel(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (pointTargetData.mat), Problem 4');

summed_channels = sum(delayed_channel,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (pointTargetData.mat), Problem 4');

load('anecoicCystData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

for beam = 1:128
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5 + depth;
    time_to_point(ii) = d(ii)/speed;
end

delay_matrix(jj,:,beam) = time_to_point; %delays

end

for aa = 1:128
    delayed_channel(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
end


end


for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_channel)));
max_data = max(max(max(delayed_channel)));
imagesc(delayed_channel(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (anecoicCystData.mat), Problem 4');

summed_channels = sum(delayed_channel,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (anecoicCystData.mat), Problem 4');

