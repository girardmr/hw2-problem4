clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs; 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

channel = [[-63.5:1:63.5]];

for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m
data_matrix = data([1:jj+150],:,:);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end
time_from_zero = time_to_point(64);
time_from_zero_v = ones(1,length(time_to_point))*time_from_zero;
time_delay = time_to_point - time_from_zero_v;

time_array = time_array_all(1:max(size(data_matrix)));

for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed = time_array+delay;
    delayed_channel([1:jj+150],aa,:) = interp1(time_array,data_matrix([1:jj+150],aa,:),time_array_delayed,'linear');
    delayed_channel_row(jj,aa,:) = delayed_channel(jj,aa,:);
end


end

for ll = 1:numel(delayed_channel_row)
    if isnan(delayed_channel_row(ll))==1
        delayed_channel_row(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_channel_row)));
max_data = max(max(max(delayed_channel_row)));
imagesc(delayed_channel_row(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (pointTargetData.mat)');

summed_channels = sum(delayed_channel_row,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (pointTargetData.mat)');

load('anecoicCystData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs; 

for ii = 1:max(size(data))
    time_array(ii) = ii/fs;
end

channel = [-63.5:1:63.5];

for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m
data_matrix = data([1:jj+150],:,:);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end
time_from_zero = time_to_point(64);
time_from_zero_v = ones(1,length(time_to_point))*time_from_zero;
time_delay = time_to_point - time_from_zero_v;

time_array = time_array_all(1:max(size(data_matrix)));

for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed = time_array+delay;
    delayed_channel([1:jj+150],aa,:) = interp1(time_array,data_matrix([1:jj+150],aa,:),time_array_delayed,'linear');
    delayed_channel_row(jj,aa,:) = delayed_channel(jj,aa,:);
end

end

for ll = 1:numel(delayed_channel_row)
    if isnan(delayed_channel_row(ll))==1
        delayed_channel_row(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_channel_row)));
max_data = max(max(max(delayed_channel_row)));
imagesc(delayed_channel_row(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (anecoicCystData.mat)');

summed_channels = sum(delayed_channel_row,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (anecoicCystData.mat)');

