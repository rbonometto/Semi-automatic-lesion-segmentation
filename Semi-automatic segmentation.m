% Local threshold algorithm
clc, clear all, close all
load MRIdata.mat

map = gray;
axial = vol;

% Montage axial slices
figure()
montage(axial, map)
title('Axial Slices')

% Montage sagittal slices
for i = 1:size(vol,1)
    % sagittal(:,:,i) = vol(i,:,:);
    slice = vol(i,:,:);
    slice = reshape(slice,[256 112]);
%     T0 = maketform('affine', [0 -2.5; 1.7 0; 0 0]);
%     R2 = makeresampler({'cubic','nearest'},'fill');
%     sagittal(:,:,i) = imtransform(slice,T0,R2);
    sagittal(:,:,i) = imrotate(slice,90);
end

figure()
montage(sagittal,map)
title('Sagittal Slices')


%% Extracting sagittal slice 135 to select the ROI
disp('SAGITTAL SLICES PROCESSING')
M = sagittal(:,:,135);
figure()
imshow(M,map); 
disp('Drag and drop to select the ROI; then, right click + "CropImage"')
[M_crop, coord] = imcrop(M); 
% [xmin ymin width height]
coord(1:2) = floor(coord(1:2));
coord(3:4) = ceil(coord(3:4));
x_min_sag = coord(1); y_min_sag = coord(2); width_sag = coord(3); height_sag =coord(4); 

close

%% Skull stripping and crop of the ROI for every sagittal slice
disp('Performing skull stripping on sagittal slices...')
skull_threshold = 107; % Threshold for skull stripping
for j = 1:size(sagittal,3)
     % Skull stripping
     binaryImage = sagittal(:,:,j) > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = sagittal(:,:,j);
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min_sag : y_min_sag+height_sag, x_min_sag : x_min_sag+width_sag);
     ROI_vol_sag(:,:,j) = skullFreeImage;        
end

%% Workflow for lesion segmentation --> sagittal slices
disp('Performing the processing of sagittal slices... ')
k = 7; % Constant for histogram smoothing
T = [];
masks_sag = {};

for j = 107:144 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_sag(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 245) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    masks_sag{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');

    figure('visible', 'off')
    subplot(221), imshow(variable), title('Original ROI')
    subplot(222), imshow(masks_sag{j}), title('Binarized ROI')
    subplot(223), bar(new_hist), xline([loc_peak1 loc_peak2 T(j)]), title('Histogram')
    subplot(224), imshow(BW), title('Edge detection with Sobel filter')
    
end

%% SELECTING AND POSTPROCESSING OF THE MASKS
masks_sag_corr = {};
for i = 107: length(masks_sag)
    masks_sag_corr{i} = bwareafilt(masks_sag{i},1);
    masks_sag_corr{i} = imfill(masks_sag_corr{i}, 'holes');  
    
    figure('visible', 'off')
    subplot(131),imshow(ROI_vol_sag(:,:,i)), title('Original ROI')
    subplot(132), imshow(masks_sag{i}), title('Mask selected by thresholding')
    subplot(133), imshow(masks_sag_corr{i}), title('Filled mask')
end 

%% FREEHAND SEGMENTATION FOR SOME SLICES
disp('Please, manually select the tumor on the following slices: ')
% Adjusting the contrast to facilitate manual segmentation
contrasto = (ROI_vol_sag(:,:,107));

figure()
imshow(contrasto), title('Manually segment the tumor:')
m = drawfreehand;
bw = createMask(m);
masks_sag_corr{107}=bw;

figure('visible', 'off')
subplot(131), imshow(ROI_vol_sag(:,:,107)), title('Original ROI')
subplot(132), imshow(masks_sag{107}), title('Mask selected by thresholding')
subplot(133), imshow(bw), title('Mask selected by hand')

contrasto2 = (ROI_vol_sag(:,:,144));
figure()
imshow(contrasto2), title('Manually segment the tumor:')
m2 = drawfreehand;
bw2 = createMask(m2);
masks_sag_corr{144}=bw2;

figure('visible','off')
subplot(131), imshow(ROI_vol_sag(:,:,144)), title('Original ROI')
subplot(132), imshow(masks_sag{144}), title('Mask selected by thresholding')
subplot(133), imshow(bw2), title('Mask selected by hand')

for i = 128:131
    contrast = (ROI_vol_sag(:,:,i));
    figure()
    imshow(contrast), title('Manually segment the tumor:')
    m3 = drawfreehand;
    bw3 = createMask(m3);
    masks_sag_corr{i} = bw3;
    figure('visible', 'off')
    subplot(131)
    imshow(ROI_vol_sag(:,:,i)), title('Original ROI')
    subplot(132)
    imshow(masks_sag{i}), title('Mask selected by thresholding')
    subplot(133)
    imshow(bw3), title('Mask selected by hand')

end

close all

%% PLOTTING THE FINAL MASKS SUBSEQUENTLY USED FOR THE COMPUTATION OF AREA AND VOLUME

for i = 107:144
    figure('visible', 'off')
    imshow(masks_sag_corr{i})
end

%% COMPUTATION OF THE AREA FOR EACH SLICE 
disp('Computing the area of the tumor on each sagittal slice...')
area_sag = [];
for k = 107:length(masks_sag_corr)
    hh = regionprops(masks_sag_corr{k}, 'Area');
    area_sag = [area_sag hh.Area];
end

area_mm = area_sag.*(pixdim(1)).*(pixdim(2));

%% COMPUTATION OF THE VOLUME
disp('Computing the volume of the tumor from the sagittal slices...')
volume_mm = 0;
for k = 1:length(area_mm)
    volume_slice = area_mm(k)*pixdim(3);
    volume_mm = volume_mm + volume_slice;
end

fprintf('Volume computed from sagittal slices: %f mm^3\n', volume_mm);

%% AXIAL SLICES ANALYSIS
disp('AXIAL SLICES PROCESSING')
% Selecting the ROI on slice 77
M = axial(:,:,77);
figure;
disp('Drag and drop to select the ROI; then, right click + "CropImage"')
[M_crop, coord] = imcrop(M);
coord(1:2) = floor(coord(1:2));
coord(3:4) = ceil(coord(3:4));
x_min = coord(1);
y_min = coord(2);
width = coord(3);
height=coord(4);

close all

%% Skull stripping and crop of the ROI for every axial slice
disp('Performing skull stripping on axial slices...')
skull_threshold = 107; % Threshold for skull stripping
for j = 1:size(axial,3)
     % Skull stripping
     binaryImage = axial(:,:,j) > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = axial(:,:,j);
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min : y_min+height, x_min : x_min+width);
     ROI_vol_ax(:,:,j) = skullFreeImage;        
end

%% Workflow for lesion segmentation --> axial slices
disp('Performing the processing of axial slices... ')
k = 7; 
T = [];
masks_ax = {};

for j = 65:89 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_ax(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 250) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    masks_ax{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');

    figure('visible', 'off')
    subplot(221), imshow(variable), title('Original ROI')
    subplot(222), imshow(masks_ax{j}), title('Binarized ROI')
    subplot(223), bar(new_hist); xline([loc_peak1 loc_peak2 T(j)]), title('Histogram')
    subplot(224), imshow(BW), title('Edge detection with Sobel filter')

end

%% SELECTING AND POSTPROCESSING OF THE MASKS
for i = 65: length(masks_ax)
    masks_ax_corr{i} = bwareafilt(masks_ax{i},1);
    masks_ax_corr{i} = imfill(masks_ax_corr{i}, 'holes'); 

    figure('visible', 'off')
    subplot(131), imshow(ROI_vol_ax(:,:,i)), title('Original ROI')
    subplot(132), imshow(masks_ax{i}), title('Mask selected by thresholding')
    subplot(133), imshow(masks_ax_corr{i}), title('Filled mask')
end 


%% PLOTTING THE FINAL MASKS SUBSEQUENTLY USED FOR THE COMPUTATION OF AREA AND VOLUME

for i = 65:length(masks_ax_corr)
    figure('visible', 'off')
    imshow(masks_ax_corr{i})
end

%% COMPUTATION OF THE AREA FOR EACH SLICE 
disp('Computing the area of the tumor on each sagittal slice...')
area_ax = [];
for k = 65:length(masks_ax_corr)
    hhax = regionprops(masks_ax_corr{k}, 'Area');
    area_ax = [area_ax hhax.Area];
end

area_mm_ax = area_ax.*(pixdim(1)).*(pixdim(2));

%% COMPUTATION OF THE VOLUME
disp('Computing the volume of the tumor from the sagittal slices...')
volume_mm_ax = 0;
for k = 1:length(area_mm_ax)
    volume_sliceax = area_mm_ax(k)*pixdim(3);
    volume_mm_ax = volume_mm_ax + volume_sliceax;
end

fprintf('Volume computed from axial slices: %f mm^3\n', volume_mm_ax);


%% EVALUATING THE PERFORMANCES OF THE IMPLEMENTED WORKFLOW WITH RESPECT TO DIFFERENT LEVELS OF NOISE
GN1={};
GN2={};
SP1={};
SP2={};
for i=1:112
    % adding noise to axial slices
    GN1 = [GN1 imnoise(axial(:,:,i), 'gaussian', 0, 1e-3)];
    GN2 = [GN2 imnoise(axial(:,:,i), 'gaussian', 0, 0.1)];
    SP1 = [SP1 imnoise(axial(:,:,i), 'salt & pepper', 0.05)];
    SP2 = [SP2 imnoise(axial(:,:,i), 'salt & pepper', 0.2)];
end

GN1s={};
GN2s={};
SP1s={};
SP2s={};
for i = 1:256
    % adding noise to sagittal slices
    GN1s = [GN1s imnoise(sagittal(:,:,i), 'gaussian', 0, 1e-3)];
    GN2s = [GN2s imnoise(sagittal(:,:,i), 'gaussian', 0, 0.1)];
    SP1s = [SP1s imnoise(sagittal(:,:,i), 'salt & pepper', 0.05)];
    SP2s = [SP2s imnoise(sagittal(:,:,i), 'salt & pepper', 0.2)];
end 

% example of different types of noise superimposed to axial slices
figure()
subplot(231), imshow(vol(:,:,77), []), title('original')
subplot(232), imshow(GN1{77}, []), title('GN1')
subplot(233), imshow(GN2{77}, []), title('GN2')
subplot(234), imshow(SP1{77}, []), title('SP1')
subplot(235), imshow(SP2{77}, []), title('SP2')

% example of different types of noise superimposed to sagittal slices
figure()
subplot(231), imshow(sagittal(:,:,135), []), title('original')
subplot(232), imshow(GN1s{135}, []), title('GN1')
subplot(233), imshow(GN2s{135}, []), title('GN2')
subplot(234), imshow(SP1s{135}, []), title('SP1')
subplot(235), imshow(SP2s{135}, []), title('SP2')


%% CASE 1: GAUSSIAN NOISE, low, axial slices
for j = 1:length(GN1)
     % Skull stripping
     binaryImage = GN1{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = GN1{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min : y_min+height, x_min : x_min+width);
     ROI_vol_GN1(:,:,j) = skullFreeImage;        
end

for j = 65:89 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_GN1(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 250) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_GN1{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
   
end

for i = 65: length(masks_GN1)
    masks_GN1_corr{i} = bwareafilt(masks_GN1{i},1);
    masks_GN1_corr{i} = imfill(masks_GN1_corr{i}, 'holes');

end 

area_GN1 = [];
for k = 65:length(masks_GN1_corr)
    hhax = regionprops(masks_GN1_corr{k}, 'Area');
    area_GN1 = [area_ax hhax.Area];
end

area_mm_GN1 = area_GN1.*(pixdim(1)).*(pixdim(2));

volume_mm_GN1 = 0;
for k = 1:length(area_mm_GN1)
    volume_sliceGN1 = area_mm_GN1(k)*pixdim(3);
    volume_mm_GN1 = volume_mm_GN1 + volume_sliceGN1;
end

%% GAUSSIAN NOISE, high, axial slices
for j = 1:length(GN2)
     % Skull stripping
     binaryImage = GN2{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = GN2{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min : y_min+height, x_min : x_min+width);
     ROI_vol_GN2(:,:,j) = skullFreeImage;        
end

for j = 65:89 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_GN2(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 250) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_GN2{j} = segmented_image;
    

    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');

end

for i = 65: length(masks_GN2)
    masks_GN2_corr{i} = bwareafilt(masks_GN2{i},1);
    masks_GN2_corr{i} = imfill(masks_GN2_corr{i},'holes');
end 

area_GN2 = [];
for k = 65:length(masks_GN2_corr)
    hhax = regionprops(masks_GN2_corr{k}, 'Area');
    area_GN2 = [area_ax hhax.Area];
end

area_mm_GN2 = area_GN2.*(pixdim(1)).*(pixdim(2));

volume_mm_GN2 = 0;
for k = 1:length(area_mm_GN2)
    volume_sliceGN2 = area_mm_GN2(k)*pixdim(3);
    volume_mm_GN2 = volume_mm_GN2 + volume_sliceGN2;
end

%% SALT & PEPPER NOISE, low, axial slices
for j = 1:length(SP1)
     % Skull stripping
     binaryImage = SP1{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = SP1{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min : y_min+height, x_min : x_min+width);
     ROI_vol_SP1(:,:,j) = skullFreeImage;        
end

for j = 65:89 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_SP1(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 250) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_SP1{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
    
end

for i = 65: length(masks_SP1)
    masks_SP1_corr{i} = bwareafilt(masks_SP1{i},1);
    masks_SP1_corr{i} = imfill(masks_SP1{i},'holes');
end 

area_SP1 = [];
for k = 65:length(masks_SP1_corr)
    hhax = regionprops(masks_SP1_corr{k}, 'Area');
    area_SP1 = [area_ax hhax.Area];
end

area_mm_SP1 = area_SP1.*(pixdim(1)).*(pixdim(2));

volume_mm_SP1 = 0;
for k = 1:length(area_mm_SP1)
    volume_sliceSP1 = area_mm_SP1(k)*pixdim(3);
    volume_mm_SP1 = volume_mm_SP1 + volume_sliceSP1;
end


%% SALT AND PEPPER NOISE, high, axial slices
for j = 1:length(SP2)
     % Skull stripping
     binaryImage = SP2{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = SP2{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min : y_min+height, x_min : x_min+width);
     ROI_vol_SP2(:,:,j) = skullFreeImage;        
end

for j = 65:89 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_SP2(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 250) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_SP2{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
    
end

for i = 65: length(masks_SP2)
    masks_SP2_corr{i} = bwareafilt(masks_SP2{i},1);
    masks_SP2_corr{i} = imfill(masks_SP2{i},'holes');
end 

area_SP2 = [];
for k = 65:length(masks_SP2_corr)
    hhax = regionprops(masks_SP2_corr{k}, 'Area');
    area_SP2 = [area_ax hhax.Area];
end

area_mm_SP2 = area_SP2.*(pixdim(1)).*(pixdim(2));

volume_mm_SP2 = 0;
for k = 1:length(area_mm_SP2)
    volume_sliceSP2 = area_mm_SP2(k)*pixdim(3);
    volume_mm_SP2 = volume_mm_SP2 + volume_sliceSP2;
end

%% GAUSSIAN NOISE, low, sagittal slices
for j = 1:length(GN1s)
     % Skull stripping
     binaryImage = GN1s{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = GN1s{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min_sag : y_min_sag+height_sag, x_min_sag : x_min_sag+width_sag);
     ROI_vol_GN1s(:,:,j) = skullFreeImage;        
end

for j = 107:144 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_GN1s(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 245) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_GN1s{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
   
end

%% GAUSSIAN NOISE, high, sagittal slices
for j = 1:length(GN2s)
     % Skull stripping
     binaryImage = GN2s{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = GN2s{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min_sag : y_min_sag+height_sag, x_min_sag : x_min_sag+width_sag);
     ROI_vol_GN2s(:,:,j) = skullFreeImage;        
end

for j = 107:144 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_GN2s(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 245) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_GN2s{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
    
end

%% SALT & PEPPER NOISE, low, sagittal slice
for j = 1:length(SP1s)
     % Skull stripping
     binaryImage = SP1s{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = SP1s{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min_sag : y_min_sag+height_sag, x_min_sag : x_min_sag+width_sag);
     ROI_vol_SP1s(:,:,j) = skullFreeImage;        
end

for j = 107:144 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_SP1s(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 245) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_SP1s{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');
    
end

%% SALT & PEPPER NOISE, high, sagittal slices
for j = 1:length(SP2s)
     % Skull stripping
     binaryImage = SP2s{j} > skull_threshold;
     labeledImage = bwlabel(binaryImage);
     binaryImage = ismember(labeledImage,1);
     binaryImage = imdilate(binaryImage, true(6));
     skullFreeImage = SP2s{j};
     skullFreeImage(binaryImage) = 0;
        
     % ROI extraction
     skullFreeImage = skullFreeImage(y_min_sag : y_min_sag+height_sag, x_min_sag : x_min_sag+width_sag);
     ROI_vol_SP2s(:,:,j) = skullFreeImage;        
end

for j = 107:144 % Manually selecting the slices where the tumor can be seen
    variable = ROI_vol_SP2s(:,:,j);
    variable = medfilt2(variable);
    variable(variable >= 245) = 0;

    % Histogram to find the threshold
    [counts,edges] = imhist(variable,256);
    
    % Histogram smoothing to make finding the maxima easier
    for i = k+1 :length(edges)-k
        new_hist(i) = sum(counts(i-k:i+k))/(2*k+1); 
    end 

    cropped_hist = new_hist(60:end); % To avoid a threshold that's too low
    local_maxima = islocalmax(cropped_hist, 'MinSeparation',50);
    loc_peaks = find(local_maxima);
    peak_values = cropped_hist(loc_peaks);
    peak_matrix = [loc_peaks; peak_values];
    sorted_matrix = sortrows(peak_matrix', -2); % Sorting the peaks in descending order of value
    if length(loc_peaks) > 1  
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = sorted_matrix(2,1) + 60;
    else 
        loc_peak1 = sorted_matrix(1,1) + 60;
        loc_peak2 = 200; % Last few slices only show one peak in the histogram
    end
    
    % Selecting the threshold between the two maxima
    T(j) = (loc_peak1 + loc_peak2)/2;
    segmented_image = imbinarize(variable, T(j)/255);
    % segmented_image = imfill(segmented_image, 'holes');
    masks_SP2s{j} = segmented_image;
    
    % Edge detection with the Sobel filter
    [BW, threshOut] = edge(variable,'sobel');

end

%% EXAMPLE OF THE REMOVAL OF DIFFERENT TYPES OF NOISE with a median filter
figure()
subplot(251)
imshow(ROI_vol_ax(:,:,77)), title('original ROI')
subplot(252)
imshow(ROI_vol_GN1(:,:,77)), title('low gaussian noise')
subplot(253)
imshow(ROI_vol_GN2(:,:,77)), title('high gaussian noise')
subplot(254)
imshow(ROI_vol_SP1(:,:,77)), title('low salt and pepper noise')
subplot(255)
imshow(ROI_vol_SP2(:,:,77)), title('high salt and pepper noise')
subplot(256)
imshow(medfilt2(ROI_vol_ax(:,:,77))), title('original ROI filtered with a median filter')
subplot(257)
imshow(medfilt2(ROI_vol_GN1(:,:,77))), title('low GN filtered with median filter')
subplot(258)
imshow(medfilt2(ROI_vol_GN2(:,:,77))), title('high GN filtered with median filter')
subplot(259)
imshow(medfilt2(ROI_vol_SP1(:,:,77))), title('low S&P noise filtered with median filter')
subplot(2,5,10)
imshow(medfilt2(ROI_vol_SP2(:,:,77))), title('high S&P noise filtered with median filter')


figure()
subplot(251)
imshow(ROI_vol_sag(:,:,135)), title('original ROI')
subplot(252)
imshow(ROI_vol_GN1s(:,:,135)), title('low gaussian noise')
subplot(253)
imshow(ROI_vol_GN2s(:,:,135)), title('high gaussian noise')
subplot(254)
imshow(ROI_vol_SP1s(:,:,135)), title('low salt and pepper noise')
subplot(255)
imshow(ROI_vol_SP2s(:,:,135)), title('high salt and pepper noise')
subplot(256)
imshow(medfilt2(ROI_vol_sag(:,:,135))), title('original ROI filtered with median filter')
subplot(257)
imshow(medfilt2(ROI_vol_GN1s(:,:,135))), title('low GN filtered with median filter')
subplot(258)
imshow(medfilt2(ROI_vol_GN2s(:,:,135))), title('high GN filtered with median filter')
subplot(259)
imshow(medfilt2(ROI_vol_SP1s(:,:,135))), title('low S&P noise filtered with median filter')
subplot(2,5,10)
imshow(medfilt2(ROI_vol_SP2s(:,:,135))), title('high S&P noise filtered with median filter')






