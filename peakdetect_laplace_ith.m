function [region_properties, L2, bw2, num2, mF, sF, area, cent] = peakdetect_laplace_ith(img,ith, show_control_figures)

minarea = 10;           % --> Minimale Regionsgröße in Pixeln
maxarea = 40;          % --> Maximale Regionsgröße in Pixeln
pth = 20;              % --> Percentile für die Bestimmung des thresholds
laplacian_shape = 0.2;

% Bild doublen zum Rechnen:
img = double(img);

if show_control_figures==1
    figure
    subplot(1,3,1)
    imagesc(img)
    title('Raw image')
end
% Vorfiltern des Bildes:
h = fspecial('average',3);    % Moving average - Filter
img = imfilter(img,h);
if show_control_figures==1
    subplot(1,3,2)
    imagesc(img)
    title('Image after moving average')
end


% Anwendung Laplace-Filter:
h = fspecial('laplacian', laplacian_shape);
filtered_img = imfilter(img, h);
filtered_img( filtered_img > 0 ) = 0; % Werte größer Null entfernen



% Normalize image:
filtered_img = double(filtered_img);
filtered_img = (filtered_img-min(filtered_img(:)))/(max(filtered_img(:))-min(filtered_img(:)));
filtered_img = 1-filtered_img;

if show_control_figures==1
    subplot(1,3,3)
    imagesc(filtered_img)
    title('Image after Laplace filter')
end


% Determining upper pth-th percentile of intensity values
pth = pth / 100;
[counts, bins] = imhist(filtered_img);
k = 0;
while sum(counts(end-k:end)) < pth*sum(counts)
    k = k + 1;
end
level = bins(end-k+1);

% Schwarz-weiß Bild mit 'level' als Schwelle bestimmen:
bw = im2bw(filtered_img, level);



% Remove regions, that are in contact with border:
bw = imclearborder(bw);

if show_control_figures ==1
    figure
    imagesc(bw)
    title('Binary image after laplacian filtering')
end

% Label all regions with integer numbers and measure 'Area' of all regions:
[L, num] = bwlabel(bw, 4);
region_area = regionprops(L, 'Area');
region_mean = regionprops(L,img, 'MeanIntensity');
% Take only regions which lies in between the specified values:
idx = find([region_area.Area] > minarea & [region_area.Area] < maxarea & [region_mean.MeanIntensity] > ith);
bw2 = ismember(L,idx);
clear region_area L num bw

if show_control_figures ==1
    figure
    imagesc(bw2)
    title('Binary image after filtering for size and intensity')
end

% Fill regions, if there is a hole in it:
bw2 = imfill(bw2, 'holes');
if show_control_figures ==1
    figure
    imagesc(bw2)
    title('Binary image after filtering for size and intensity with holes filled')
end
% Thicken regions by 1 pixel:
% bw2 = bwmorph(bw2,'thicken',1);
% if show_control_figures ==1
%     figure
%     imagesc(bw2)
%     title('Binary image after filtering for size and intensity with holes filled and thickened')
% end 
% Label all regions with integer numbers:
[L2, num2] = bwlabel(bw2, 4);

% Measure defined properties of the regions:
region_properties = regionprops(L2, filtered_img, 'Area', 'Centroid', 'EquivDiameter', 'WeightedCentroid','PixelIdxList','PixelList');

for k=1:length(region_properties)
    region_properties(k).summen_helligkeit=sum(img(region_properties(k).PixelIdxList));
    region_properties(k).mittlere_helligkeit=mean(img(region_properties(k).PixelIdxList));
end

%mF1,sF1,area1,cent1,
for k=1:length(region_properties)
    mF(k)=region_properties(k).mittlere_helligkeit;
    sF(k)=region_properties(k).summen_helligkeit;
    area(k)=region_properties(k).Area;
    cent(k,:)=region_properties(k).Centroid;
end
    
%testbild
if show_control_figures==1
figure
subplot(1,2,2)
imagesc(L2>0); 
title(sprintf('Found regions: %i', size(region_properties,1)))
subplot(1,2,1)
imagesc(img); 
title('Search image')
end