%% SCRIPT_NAME.m
%% Author: Erin Gilbert
%% Created: May 24 2018
%% Modified: Sept 16 2018

%% Usage: matlab /r "file=RAWFILE;file_white=RAWWHITE;file_dark=DARKFILE;matlab_soybean_pixel_extraction.m"
%% Purpose: Using thresholds (Mainly NDVI) to create a mask of an image with black being background (non-plant) and white being plant. These masks are then used to populate a JSON file with the pixel coordinates and cell values of hyperspectral images of plants.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Define Environment  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%function files required:
%parseHdrInfo.m, norm01.m, SNV_function.m


function soybean_extraction_test_kmeans(file, file_white, file_dark)
clearvars

%%%%%%%%%%%%%%%%%%%%%%%
%% Load All The Data %%
%%%%%%%%%%%%%%%%%%%%%%%

%%Assign hyperspectral .raw file to variable 'file'.
file=RAWFILE;

%parsing the title for information
filename = strsplit(file, '/');
filedirectory = join(filename(1:end-1),"/");
filedirectory = strcat(filedirectory, "/");
filename = filename{length(filename)};
fileinfo = strsplit(filename, ".");
fileinfo = fileinfo{1};
fileinfo_split = strsplit(fileinfo, "_");
potlabels = strsplit(fileinfo_split{length(fileinfo_split)}, "-");
potlabels = potlabels(2:end);


%%Get info from the header file
header= strcat(fileinfo, '.hdr');
[wavelengths, spatial, frames, spectral, tint, settings] =  parseHdrInfo(char(filedirectory), header);%MAKE DYNAMIC
wavelengths = round(wavelengths);

%%Read file with multibandread function, and assign to variable 'data'.
data = multibandread(file, [frames, spatial, spectral], 'uint16', 0, 'bil', 'ieee-le');


%%Assign white reference hyperspectral .raw file to variable 'file_white'.
file_white=RAWWHITE;

%%Assign dark hyperspectral .raw file to variable 'file_dark'.
file_dark=DARKFILE;

%%Read white reference file with multibandread function, and assign to variable 'data_w'.
data_w = multibandread(file_white, [10, spatial, spectral], 'uint16', 0, 'bil', 'ieee-le');

%%Read dark reference file with multibandread function, and assign to
data_d = multibandread(file_dark, [10, spatial, spectral], 'uint16', 0, 'bil', 'ieee-le');



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Average dark reference into 2 dimensions (width and wavelengths).
data_davg = mean(data_d, 1);

%%Average white reference into 2 dimensions (width and wavelengths).
data_wavg = mean(data_w, 1);

%%Normalize data to white and dark references.
data_norm_num = bsxfun(@minus, data, data_davg);
data_norm_denom = data_wavg - data_davg;  
data_norm = bsxfun(@rdivide, data_norm_num, data_norm_denom);

%Transpose data
data_norm_t = permute(data_norm, [2, 1, 3]);

%Calculate NDVI
c = norm01(data_norm_t(:, :, find(wavelengths == 682))); %682nm
d = norm01(data_norm_t(:, :, find(wavelengths == 801))); %801nm
NDVI = ((d-c)./(d+c));

% Create color image
R = norm01(data_norm_t(:, :, find(wavelengths == 640)));
G = norm01(data_norm_t(:, :, find(wavelengths == 473)));
B = norm01(data_norm_t(:, :, find(wavelengths == 532)));

gamma = .7;
RGBall = cat(3,R.^gamma,B.^gamma,G.^gamma);
%figure, imshow(RGBall), title('RGB'); %generate RGB image


%% K-means clustering and visualization

lab_he=rgb2lab(RGBall);


ab = lab_he(:,:,2:3);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

% 5 clusters RGB
nColors = 5;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(NDVI(:),nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
                                  
pixel_labels = reshape(cluster_idx,nrows,ncols);
%imshow(pixel_labels,[]);

segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = RGBall;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end

figure,imshow(segmented_images{1}), title('objects in cluster 1');
figure,imshow(segmented_images{2}), title('objects in cluster 2');
figure, imshow(segmented_images{3}), title('objects in cluster 3');
figure,imshow(segmented_images{4}), title('objects in cluster 4');
figure,imshow(segmented_images{5}), title('objects in cluster 5');
%figure, imshow(segmented_images{2} +segmented_images{3}), title('combined objects');


whiteImage = 255 * ones(80, frames, 3, 'uint8');

whiteImage = insertText(whiteImage, [10 10],fileinfo,'FontSize', 40,'BoxColor','w');
imshow(whiteImage);

qc_image_name = strcat('/Users/gilbe952/Desktop/',fileinfo, '_kmeans.png');

imwrite([im2double(whiteImage); RGBall;segmented_images{1};segmented_images{2};segmented_images{3}; segmented_images{4}; segmented_images{5}], qc_image_name);


%% Object identification through creating a mask


BW = im2bw(segmented_images{3}, 0);
%figure, imshow(BW);
BW2 = bwareafilt(BW, [1000 inf]);
%figure, imshow(BW2);

e = data_norm_t(:, :, 1);
h = bwconncomp(BW);

%Label the objects
L = labelmatrix(h);
colRGB = label2rgb(L);
figure, imshow(colRGB);

%%Calculate parameters for objects and display on image for the 1st wavelength.
more = regionprops(h, e, 'all');

%%%If there are more than 6 objects, combine objects that are close
while  h.NumObjects > length(potlabels)
    x = [];
    y = [];
    for i = 1:h.NumObjects
        x = [x, more(i).Centroid(1)];
        y = [y, more(i).Centroid(2)];
    end
    for xIndex = 1 : length(x)
        for yIndex = xIndex : length(y)
            distances(xIndex, yIndex) = sqrt((x(xIndex)-x(yIndex))^2 + (y(xIndex)-y(yIndex))^2);
        end
    end
    distances(distances==0)= NaN; %remove zeros from distance matrix
    [M,I] = min (distances(:));
    [I_row, I_col] = ind2sub(size(distances),I);
    
    %combine the objects
    h.PixelIdxList{I_row}=sort([h.PixelIdxList{I_row}; h.PixelIdxList{I_col}]);
    h.PixelIdxList(I_col)=[];
    h.NumObjects = h.NumObjects -1;
    more = regionprops(h, e, 'all');
end

%Label the objects
L = labelmatrix(h);
colRGB = label2rgb(L);
figure, imshow(colRGB);

whiteImage = 255 * ones(80, frames, 3, 'uint8');

whiteImage = insertText(whiteImage, [10 10],fileinfo,'FontSize', 40,'BoxColor','w');
%imshow(whiteImage);


qc_image_name = strcat('/home/hirschc3/gilbe952/sds_bsr/objects/matlab/',fileinfo, '.png');
imwrite([im2double(whiteImage);RGBall;cat(3, NDVI, NDVI, NDVI);im2double(colRGB)], qc_image_name);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Auxilary PCA Code  %%
%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Getting the data
% 
% %%Pixel values for all
% for i2 = 1:size(data_norm_t, 3)
%     all_pix(:, :, i2) = regionprops(h, data_norm_t(:, :, i2), 'PixelValues');
% end
% pixeltable = permute(all_pix, [1 3 2]);%data in struct (37x290) with 1 field (PixelValues, each with nx1 double).
% pv2 = struct2cell(pixeltable);
% table2 = permute(pv2, [2 3 1]);
% 
% %%PCA on subset of wavelengths (500-900 nm (so 56 to 227))
% pv = cell2mat(table2(:, find(wavelengths == 500):find(wavelengths == 900)));
% pv = cell2mat(table2(:, find(wavelengths == 500):find(wavelengths == 900)));
% 
% %% Preprocessing methods
% 
% %standard normal variate (scatter correction)
% 
% [pv_snv]=SNV_function(pv);
% 
% [coeff, scores, latent, tsquared, explained, mu] = pca(pv_snv);
% 
% 
% % 1st PC
% t = struct2table(more);
% idlist = t.PixelIdxList;
% IDlist = cell2mat(idlist);
% img = NDVI;
% for j = 1:length(IDlist)
%     img(IDlist(j)) = scores(j, 1);
% end
% PC1 = img.*BW2;
% %figure, imshow(PC1, []);
% 
% % 2nd PC
% img = NDVI;
% for j = 1:length(IDlist)
%     img(IDlist(j)) = scores(j, 2);
% end
% PC2 = img.*BW2;
% %figure, imshow(PC2, []);
% 
% % 3rd PC
% img = NDVI;
% for j = 1:length(IDlist)
%     img(IDlist(j)) = scores(j, 3);
% end
% PC3 = img.*BW2;
%figure, imshow(PC3, []);


%figure, imshow(cat(3, PC1.^gamma, PC2.^gamma,PC3.^gamma));
%figure, imshow(cat(3, R.^gamma, G.^gamma, d.^gamma));

% %% LOADINGS
% loadings = array2table(coeff(:, 1:3));
% 
% % plot of loadings for PC 1-3.
% figure, plot(coeff(:,1:3),'DisplayName','coeff(:,1:3)')
% 
% csvwrite('PCA_scores_SNV.csv', scores);
% csvwrite('PCA_loadings_SNV.csv', coeff);
% csvwrite('PCA_explained_SNV.csv', explained);


