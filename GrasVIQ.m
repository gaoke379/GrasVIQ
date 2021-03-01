% % %*********************************************************************************
% * Filename: GrasVIQ.m
% *
% * Description:
% * Matlab script for Grass Vein Image Quantification (GrasVIQ)
% *
% * Authors: Ke Gao, Michael Boeding, Dr. Filiz Bunyak
% *
% * Copyright (C) 2021. Ke Gao, Michael Boeding, Dr. Filiz Bunyak
% * and Curators of the University of Missouri, a public corporation.
% * All Rights Reserved.
% *
% * Created by
% * Ke Gao
% * Department of Electrical Engineering and Computer Science
% * University of Missouri-Columbia
% * Columbia, MO 65211
% * kegao@mail.missouri.edu
% *
% *
% * For more information, please contact;
% * Dr. Filiz Bunyak
% * Department of Electrical Engineering and Computer Science
% * 219 Naka Hall
% * University of Missouri-Columbia
% * Columbia, MO 65211-2060
% * bunyak@missouri.edu
% *
% % %*********************************************************************************


clc; clear;
close all;
% clf;

% %-- Parameters
IS_JUVENILE = 1; % set to 1 if processing juvenile leaf images or 0 for adult leaf images (1: juvenile; 0: adult)
DEBUG = 0; % visualize intermediate results for debugging
EXPORT_CSV = 1; % export extracted measurements to csv files
EXPORT_IMG = 1; % save segmentation and classification images to png files
IMG_FILE_PATH = './data/'; % path for input images
IMG_FILE_TYPE = 'png'; % file extension (e.g. png, jpg, jpeg)
SCALE_FACTOR = 0.6167; % pixels per micrometer (px/um); assume all images in the given input path share the same scale

% %-- Image file path
if ~exist(IMG_FILE_PATH, 'dir')  % check if image path is valid
   msg = 'ERROR: invalid image path!';
   error(msg);
end
imgDir = dir(fullfile(IMG_FILE_PATH, strcat('*.', IMG_FILE_TYPE)));
nImages = length(imgDir);
disp(IMG_FILE_PATH);
if (nImages == 0)  % throw an error if no images with specified extension are found
   msg = ['ERROR: no ', IMG_FILE_TYPE, ' images found!'];
   error(msg);
end

% %-- Output file path
savePath = fullfile(IMG_FILE_PATH, 'GrasVIQ_Results/');
if ~exist(savePath, 'dir')
    mkdir(IMG_FILE_PATH, 'GrasVIQ_Results');
end
csvNameVeinWidth = strcat(savePath, 'vein_thickness.csv');
csvNameAllStats = strcat(savePath, 'measurements.csv');

for fr = 1:nImages
    
    % %-- Load the input image
    imgRawName = imgDir(fr).name;
    imgPath = fullfile(IMG_FILE_PATH, imgRawName);
    X = ['Processing:  ', num2str(fr), '/', num2str(nImages), ' - ', imgRawName];
    disp(X);
    img = imread(imgPath);
    if (DEBUG)
        figure(1); subplot(2,2,1); imshow(img); title('Input Image');
    end
    imgGray = rgb2gray(img);
    imgGray = double(imgGray);
    [rowsOrig, colsOrig] = size(imgGray);
    
    % %-- Pre-process the grayscale image
    % %-- Correct nonuniform illumination
    se0 = strel('disk', round(min(rowsOrig, colsOrig)/12));
    imgGrayBG = imopen(imgGray,se0); % remove all the foreground using morphological opening
    imgGray = imgGray - imgGrayBG; % subtract the BG approximation image from the original image to correct nonuniform illumination
    
    % %-- Add padding to the grayscale image
    patchSize = round(rowsOrig/40); % height of each patch
    imgGrayPad = padarray(imgGray, [patchSize 0], 'replicate', 'both');
    [rowsPadImg, ~] = size(imgGrayPad);
    
    % %-- Process each patch
    imgEdgeBwPad = zeros(size(imgGrayPad));
    patchRowMin = 1:1:rowsPadImg-patchSize+1; % consecutive patches
    for i = 1:length(patchRowMin)
        % %-- patch crop and normalization
        patch = imgGrayPad(patchRowMin(i):patchRowMin(i)+patchSize-1, :);
        patchNorm = 255 - (255 * (patch - min(patch(:))) / (max(patch(:)) - min(patch(:))));
        % %-- vertical projection and normalization
        patchSumV = sum(patchNorm, 1);
        patchSumV = patchSumV / max(patchSumV);
        % %-- binarization using adaptive thresholding
        patchProjV = AdaptiveThreshold(patchSumV, round(min(rowsOrig, colsOrig)/17), 0.1, 0);
        % %-- find edge
        patchEdge = abs(diff(patchProjV));
        patchEdgeCol = find(patchEdge);
        patchEdgeRow = round(patchRowMin(i)+patchSize/2) * ones(length(patchEdgeCol),1);
        imgEdgeBwPad(patchEdgeRow, patchEdgeCol') = 255;
    end
    imgEdgeBw = imgEdgeBwPad(patchSize+1:patchSize+rowsOrig, :); % remove the padding
    
    % %-- Refine the detected vein edges
    se1 = strel('disk', 1);
    se2 = strel('disk', 5);
    imgEdgeBw = imdilate(imgEdgeBw, se1);
    if (min(rowsOrig, colsOrig) > 500) % apply morphologically close only if the image is large enough
        imgEdgeBw = imclose(imgEdgeBw, se2);
    end
    % %-- scatter plot on the original image
    if (DEBUG)
        figure(1); subplot(2,2,2); imshow(img); hold on; title('Edge Detection');
        [veinEdgeRows, veinEdgeCols] = find(imgEdgeBw);
        scatter(veinEdgeCols, veinEdgeRows, 'r', '.');
        hold off;
    end
    
    % %-- Connected component labeling
    imgVertVeinMask = zeros(size(imgGray));
    imgNonEdgeBw = 255 - imgEdgeBw;
    imgNonEdgeBw = bwareaopen(imgNonEdgeBw, (round(min(rowsOrig, colsOrig)/50))^2); % remove small objects
    if (IS_JUVENILE)
        se3 = strel('rectangle', [25,1]); % juvenile leaf
    else
        se3 = strel('rectangle', [round(min(rowsOrig, colsOrig)/15),1]); % adult leaf
    end
    imgNonEdgeBw = imopen(imgNonEdgeBw, se3); % remove the bridges connecting two vertical veins
    imgLabel = bwlabel(imgNonEdgeBw);
    imgGrayNorm = (255 * (imgGray - min(imgGray(:))) / (max(imgGray(:)) - min(imgGray(:)))); % normalize grayscale image
    labelStats = regionprops('table', logical(imgNonEdgeBw), 'Area');
    labelAreas = labelStats.Area;
    if (IS_JUVENILE)
        areaThresh = median(labelAreas) * 0.5; % juvenile leaf
    else
        areaThresh = median(labelAreas) * 0.8; % adult leaf
    end
    
    % %-- Compute foreground/background mask
    [histCounts,~] = imhist(uint8(imgGrayNorm),16); % calculate histogram
    fgThresh = otsuthresh(histCounts) * 255; % compute a global threshold for normalized grayscale image
    fgAreaList = [];
    for i = 1:length(labelAreas)
        if (labelAreas(i) < areaThresh)
            imgLabel(imgLabel==i) = 0; % remove small blobs
        else
            % %-- check intensity
            labelMean = mean(mean(imgGrayNorm(imgLabel==i)));
            if (labelMean < fgThresh)
                imgVertVeinMask(imgLabel==i) = 1;
                fgAreaList = [fgAreaList; labelAreas(i)];
            end
        end
    end
    imgVertVeinMaskOverlay = MarkFG(img, imgVertVeinMask, [255 0 0], [0.5 0.5]);
    if (DEBUG)
        figure(1); subplot(2,2,3); imshow(imgVertVeinMaskOverlay); title('Vertical Veins');
    end
    
    if (IS_JUVENILE)
        % %------------------- Juvenile leaf image ------------------% %
        % %-- Detect over-segmented secondary veins and quaternary veins
        % %-- Refine the FG mask
        imgBgMask = 1 - imgVertVeinMask;
        imgBgLabel = bwlabel(imgBgMask); % connected component labeling for BG regions
        bgLabelStats = regionprops('table', logical(imgBgMask), 'Area'); % compute area for BG regions
        bgLabelAreas = bgLabelStats.Area;
        bgAreaThresh = mean(bgLabelAreas) * 3; % threshold for BG region area (skip if larger)
        imgSecVeinMask = zeros(size(imgGrayNorm)); % mask for secondary veins
        imgQuatVeinMask = zeros(size(imgGrayNorm)); % mask for quaternary veins
        quatVeinDistList = []; % list of vertical distance between quaternary veins
        quatVeinWidthThresh = 2; % threshold for quaternary vein thickness
        if (min(rowsOrig, colsOrig) > 500)
            quatVeinWidthThresh = 5; % increase the threshold if image size is big
        end
        % %-- Process each BG region
        for i = 1:length(bgLabelAreas)
            if (bgLabelAreas(i) < bgAreaThresh) % skip large BG regions (under-segmented veins may exist)
                bgLabelMean = mean(mean(imgGrayNorm(imgBgLabel==i)));
                bgLabelStd = std(imgGrayNorm(imgBgLabel==i));
                if ((bgLabelMean < fgThresh*1.2) && (bgLabelStd < 13)) % secondary veins (low mean & low std)
                    imgVertVeinMask(imgBgLabel==i) = 1;
                    imgSecVeinMask(imgBgLabel==i) = 1;
                else
                    % %-- detect quaternary veins
                    imgGrayNormBgOneLabel = imgGrayNorm;
                    imgGrayNormBgOneLabel(imgBgLabel~=i) = 0;
                    rowSumBgOneLabel = sum(imgGrayNormBgOneLabel,2);
                    rowSumBgOneLabel = rowSumBgOneLabel / max(rowSumBgOneLabel);
                    rowSumBgOneLabel = smooth(rowSumBgOneLabel, round(min(rowsOrig, colsOrig)/34)); % smooth the 1-D signal
                    bgOneLabelProjH = AdaptiveThreshold(rowSumBgOneLabel, round(min(rowsOrig, colsOrig)/9), 0.2, 0); % binarization using adaptive thresholding
                    bgOneLabelProjH = 1 - bgOneLabelProjH;
                    bgOneLabelProjLabel = bwlabel(bgOneLabelProjH);
                    bgOneLabelProjLabelStats = regionprops('table', logical(bgOneLabelProjH), 'Area', 'Centroid');
                    bgOneLabelProjLabelAreas = bgOneLabelProjLabelStats.Area;
                    bgOneLabelProjLabelCentroid = bgOneLabelProjLabelStats.Centroid;
                    maskTemp = zeros(size(imgGrayNorm));
                    for i1 = 1:length(bgOneLabelProjLabelAreas)
                        if (bgOneLabelProjLabelAreas(i1) > quatVeinWidthThresh) % only keep wide local peaks (potentially quaternary veins)
                            rowIdxTemp = round(bgOneLabelProjLabelCentroid(i1,2));
                            maskTemp(max(1,rowIdxTemp-floor(bgOneLabelProjLabelAreas(i1)/2)):min(size(imgGrayNorm,1),rowIdxTemp+floor(bgOneLabelProjLabelAreas(i1)/2)),:) = 1;
                        end
                    end
                    imgBgOneLabelQuatVeinMask = zeros(size(imgGrayNorm));
                    imgBgOneLabelQuatVeinMask(imgBgLabel==i) = 1;
                    imgBgOneLabelQuatVeinMask = and(imgBgOneLabelQuatVeinMask, maskTemp); % detected quaternary vein in the current BG area
                    se4 = strel('rectangle', [round(min(rowsOrig, colsOrig)/55),1]);
                    imgBgOneLabelQuatVeinMask = imclose(imgBgOneLabelQuatVeinMask, se4); % merge over-segmented quaternary veins
                    imgQuatVeinMask = or(imgQuatVeinMask, imgBgOneLabelQuatVeinMask);
                    % %-- compute quaternary vein interval distance (vertically)
                    quatVeinOneLabel = bwlabel(imgBgOneLabelQuatVeinMask); % connected component labeling for quaternary veins in the current BG area
                    if (max(unique(quatVeinOneLabel)) > 1) % only compute distance if two or more quaternary veins exist in the current BG area
                        quatVeinOneLabelStats = regionprops('table', logical(quatVeinOneLabel), 'Centroid');
                        quatVeinOneLabelCentroidY = sort(quatVeinOneLabelStats.Centroid(:,2)); % sort the Y coordinate of quaternary veins
                        for i1 = 2:size(quatVeinOneLabelCentroidY,1)
                            quatVeinDist = quatVeinOneLabelCentroidY(i1) - quatVeinOneLabelCentroidY(i1-1); % compute distance between two adjacent quaternary veins
                            quatVeinDistList = [quatVeinDistList; quatVeinDist];
                        end
                    end
                end
            end
        end       
        % %-- Refine secondary veins mask
        imgVertVeinLabel = bwlabel(imgVertVeinMask); % connected component labeling for all vertical veins
        secVeinIDs = unique(imgSecVeinMask.*imgVertVeinLabel); % find ID for possible secondary veins
        imgSecVeinMask = zeros(size(imgGrayNorm)); % clear secondary vein mask
        areaRatioCHThresh = 0.9; % threshold for ratio of vein area to convex hull area
        for i = 1:length(secVeinIDs)
            secVeinID = secVeinIDs(i);
            if (secVeinID == 0)
                continue; % skip ID = 0 which is background
            end
            % %-- check if the current secondary vein is valid
            currSecVeinMask = zeros(size(imgGrayNorm));
            currSecVeinMask(imgVertVeinLabel==secVeinID) = 1; % mask for the current secondary vein
            areaCurrSecVein = sum(currSecVeinMask(:));
            currSecVeinMask(1,:) = 1;
            currSecVeinMask(end,:) = 1;
            currSecVeinCH = imfill(currSecVeinMask, 'holes'); % simulate convex hull by filling the holes inside the vein boundary
            currSecVeinCHShort = currSecVeinCH(2:end-1,:);
            currSecVeinCH = padarray(currSecVeinCHShort, [1 0], 'replicate', 'both');
            areaRatioCH = areaCurrSecVein / sum(currSecVeinCH(:)); % compute the ratio of the vein area to convex hull area
            if (areaRatioCH > areaRatioCHThresh) % valid secondary vein if the ratio is above a threshold
                imgSecVeinMask(imgVertVeinLabel==secVeinID) = 1; % update secondary vein mask
            end
        end
        imgSecVeinLabel = bwlabel(imgSecVeinMask); % connected component labeling for secondary veins
        numSecVein = max(unique(imgSecVeinLabel)); % number of secondary veins
        %-- compute secondary vein thickness
        vertVeinWidthList = []; % width/thickness for each vertical vein (secondary + tertiary veins including abnormal veins)
        secVeinLabelStats = regionprops('table', logical(imgSecVeinLabel), 'Area', 'Extrema');
        secVeinLabelAreas = secVeinLabelStats.Area; % area for each secondary vein
        secVeinLabelExtrema = secVeinLabelStats.Extrema; % extrema points for each secondary vein
        for i = 1:max(unique(imgSecVeinLabel)) % process each secondary vein
            veinExtrema = secVeinLabelExtrema{i,1};
            veinHeight = max(veinExtrema(:,2)) - min(veinExtrema(:,2)); % vein height
            veinWidth = round(secVeinLabelAreas(i) / veinHeight); % vein width/thickness
            vertVeinWidthList = [vertVeinWidthList veinWidth];
        end
    else
        
        % %------------------- Adult leaf image ------------------% %
        % %-- Detect quaternary veins
        % %-- Refine the FG mask
        imgBgMask = 1 - imgVertVeinMask;
        imgBgLabel = bwlabel(imgBgMask); % connected component labeling for BG regions
        bgLabelStats = regionprops('table', logical(imgBgMask), 'Area'); % compute area for BG regions
        bgLabelAreas = bgLabelStats.Area;
        bgAreaThresh = mean(fgAreaList) * 2; % threshold for BG region area (skip if larger)
        imgQuatVeinMask = zeros(size(imgGrayNorm)); % mask for quaternary veins
        quatVeinDistList = []; % list of vertical distance between quaternary veins
        % %-- Process each BG region
        for i = 1:length(bgLabelAreas)
            if (bgLabelAreas(i) < bgAreaThresh) % skip large BG regions (under-segmented veins may exist)
                % %-- detect quaternary veins
                imgGrayNormBgOneLabel = imgGrayNorm;
                imgGrayNormBgOneLabel(imgBgLabel~=i) = 0;
                rowSumBgOneLabel = sum(imgGrayNormBgOneLabel,2);
                rowSumBgOneLabel = rowSumBgOneLabel / max(rowSumBgOneLabel);
                rowSumBgOneLabel = smooth(rowSumBgOneLabel, 30); % smooth the 1-D signal
                bgOneLabelProjH = AdaptiveThreshold(rowSumBgOneLabel, 120, 0.2, 0); % binarization using adaptive thresholding
                bgOneLabelProjH = 1 - bgOneLabelProjH;
                bgOneLabelProjLabel = bwlabel(bgOneLabelProjH);
                bgOneLabelProjLabelStats = regionprops('table', logical(bgOneLabelProjH), 'Area', 'Centroid');
                bgOneLabelProjLabelAreas = bgOneLabelProjLabelStats.Area;
                bgOneLabelProjLabelCentroid = bgOneLabelProjLabelStats.Centroid;
                maskTemp = zeros(size(imgGrayNorm));
                for i1 = 1:length(bgOneLabelProjLabelAreas)
                    if (bgOneLabelProjLabelAreas(i1) > 5) % only keep wide local peaks (potentially quaternary veins)
                        rowIdxTemp = round(bgOneLabelProjLabelCentroid(i1,2));
                        maskTemp(max(1,rowIdxTemp-floor(bgOneLabelProjLabelAreas(i1)/2)):min(size(imgGrayNorm,1),rowIdxTemp+floor(bgOneLabelProjLabelAreas(i1)/2)),:) = 1;
                    end
                end
                imgBgOneLabelQuatVeinMask = zeros(size(imgGrayNorm));
                imgBgOneLabelQuatVeinMask(imgBgLabel==i) = 1;
                imgBgOneLabelQuatVeinMask = and(imgBgOneLabelQuatVeinMask, maskTemp); % detected quaternary vein in the current BG area
                se4 = strel('rectangle', [19,1]);
                imgBgOneLabelQuatVeinMask = imclose(imgBgOneLabelQuatVeinMask, se4); % merge over-segmented quaternary veins
                imgQuatVeinMask = or(imgQuatVeinMask, imgBgOneLabelQuatVeinMask);
                % %-- compute quaternary vein interval distance (vertically)
                quatVeinOneLabel = bwlabel(imgBgOneLabelQuatVeinMask); % connected component labeling for quaternary veins in the current BG area
                if (max(unique(quatVeinOneLabel)) > 1) % only compute distance if two or more quaternary veins exist in the current BG area
                    quatVeinOneLabelStats = regionprops('table', logical(quatVeinOneLabel), 'Centroid');
                    quatVeinOneLabelCentroidY = sort(quatVeinOneLabelStats.Centroid(:,2)); % sort the Y coordinate of quaternary veins
                    for i1 = 2:size(quatVeinOneLabelCentroidY,1)
                        quatVeinDist = quatVeinOneLabelCentroidY(i1) - quatVeinOneLabelCentroidY(i1-1); % compute distance between two adjacent quaternary veins
                        quatVeinDistList = [quatVeinDistList; quatVeinDist];
                    end
                end
            end
        end       
        % %-- Detect secondary veins
        imgSecVeinMask = zeros(size(imgGrayNorm)); % mask for secondary veins
        imgVertVeinLabel = bwlabel(imgVertVeinMask); % connected component labeling for all vertical veins
        vertVeinLabelStats = regionprops('table', logical(imgVertVeinLabel), 'Area', 'Extrema');
        vertVeinLabelAreas = vertVeinLabelStats.Area; % area for each vertical vein
        vertVeinAreaThresh = median(vertVeinLabelAreas) * 1.5; % area threshold to classify a secondary vein
        areaRatioCHThresh = 0.9; % threshold for ratio of vein area to convex hull area
        for i = 1:length(vertVeinLabelAreas) % check the area of each vertical vein
            if (vertVeinLabelAreas(i) > vertVeinAreaThresh) % secondary vein if area is larger than a threshold
                currSecVeinMask = zeros(size(imgGrayNorm));
                currSecVeinMask(imgVertVeinLabel==i) = 1;
                areaCurrSecVein = sum(currSecVeinMask(:));
                currSecVeinMask(1,:) = 1;
                currSecVeinMask(end,:) = 1;
                currSecVeinCH = imfill(currSecVeinMask, 'holes'); % simulate convex hull by filling the holes inside the vein boundary
                currSecVeinCHShort = currSecVeinCH(2:end-1,:);
                currSecVeinCH = padarray(currSecVeinCHShort, [1 0], 'replicate', 'both');
                areaRatioCH = areaCurrSecVein / sum(currSecVeinCH(:)); % compute the ratio of the vein area to convex hull area
                if (areaRatioCH > areaRatioCHThresh) % valid secondary vein if the ratio is above a threshold
                    imgSecVeinMask(imgVertVeinLabel==i) = 1; % update secondary vein mask
                end
            end
        end
        imgSecVeinLabel = bwlabel(imgSecVeinMask); % connected component labeling for secondary veins
        numSecVein = max(unique(imgSecVeinLabel)); % number of secondary veins
        %-- compute secondary vein thickness
        vertVeinWidthList = []; % width/thickness for each vertical vein (secondary + tertiary veins including abnormal veins)
        secVeinLabelStats = regionprops('table', logical(imgSecVeinLabel), 'Area', 'Extrema');
        secVeinLabelAreas = secVeinLabelStats.Area; % area for each secondary vein
        secVeinLabelExtrema = secVeinLabelStats.Extrema; % extrema points for each secondary vein
        for i = 1:max(unique(imgSecVeinLabel)) % process each secondary vein
            veinExtrema = secVeinLabelExtrema{i,1};
            veinHeight = max(veinExtrema(:,2)) - min(veinExtrema(:,2)); % vein height
            veinWidth = round(secVeinLabelAreas(i) / veinHeight); % vein width/thickness
            vertVeinWidthList = [vertVeinWidthList veinWidth];
        end
    end
    
    % %-- Check each tertiary vein and label abnormal veins (special cases)
    imgAbnormalVeinMask = zeros(size(imgGrayNorm)); % abnormal vein mask
    imgTerVeinMask = double(xor(imgSecVeinMask, imgVertVeinMask)); % tertiary vein mask
    imgTerVeinLabel = bwlabel(imgTerVeinMask); % connected component labeling for tertiary veins
    for i = 1:max(unique(imgTerVeinLabel)) % check each tertiary vein
        currTerVeinMask = zeros(size(imgGrayNorm));
        currTerVeinMask(imgTerVeinLabel==i) = 1; % mask for the current tertiary vein
        currVeinMask = currTerVeinMask; % a copy of current tertiary vein mask
        areaCurrTerVein = sum(currTerVeinMask(:));
        currTerVeinMask(1,:) = 1;
        currTerVeinMask(end,:) = 1;
        currTerVeinCH = imfill(currTerVeinMask, 'holes'); % simulate convex hull by filling the holes inside the vein boundary
        currTerVeinCHShort = currTerVeinCH(2:end-1,:);
        currTerVeinCH = padarray(currTerVeinCHShort, [1 0], 'replicate', 'both');
        areaRatioCH = areaCurrTerVein / sum(currTerVeinCH(:));
        if (areaRatioCH < areaRatioCHThresh) % splitting vein if the ratio is below a threshold
            imgAbnormalVeinMask(imgTerVeinLabel==i) = 1;
            imgTerVeinMask = double(xor(currVeinMask, imgTerVeinMask)); % remove splitting vein from tertiary vein mask
        else
            if (i == 1 || i == max(unique(imgTerVeinLabel))) % skip the first and the last veins
                continue;
            end
            % %-- check if the current vein is an incomplete vein
            currVeinLabel = bwlabel(currVeinMask); % connected component labeling for current tertiary vein
            currVeinStats = regionprops('table', logical(currVeinLabel), 'Extrema');
            currVeinExtremaY = currVeinStats.Extrema{1,1}(:,2);
            if ((max(currVeinExtremaY)-min(currVeinExtremaY)) < rowsOrig*0.9) % incomplete vein if it is short
                imgAbnormalVeinMask(currVeinMask==1) = 1;
                imgTerVeinMask = double(xor(currVeinMask, imgTerVeinMask)); % remove incomplete vein from tertiary vein mask
            end
        end
    end
    
    % %-- Compute distance between vertical veins
    veinDistList = [];
    imgBgMask = 1 - imgVertVeinMask;
    imgBgLabel = bwlabel(imgBgMask); % connected component labeling for BG regions
    bgLabelStats = regionprops('table', logical(imgBgLabel), 'Area', 'Extrema');
    bgLabelAreas = bgLabelStats.Area; % area for each BG region
    bgLabelExtrema = bgLabelStats.Extrema; % extrema points for each BG region
    numBgRegion = length(bgLabelAreas); % number of BG region
    for i = 1:numBgRegion
        bgRegionExtrema = bgLabelExtrema{i,1};
        bgRegionHeight = max(bgRegionExtrema(:,2)) - min(bgRegionExtrema(:,2)); % BG region height
        veinDist = round(bgLabelAreas(i) / bgRegionHeight); % distance between veins
        veinDistList = [veinDistList veinDist];
    end
    veinDistMean = mean(veinDistList); % mean distance between veins
    
    % %-- Visualize the detected veins
    imgVeinLabel = uint8(zeros(size(imgGray))); % label image for detected veins
    imgVeinLabel(imgSecVeinMask==1) = 255; % secondary vein label: 255
    imgVeinLabel(imgTerVeinMask==1) = 192; % tertiary vein label: 192
    imgVeinLabel(imgQuatVeinMask==1) = 96; % quaternary vein label: 96
    imgVeinLabel(imgAbnormalVeinMask==1) = 48; % abnormal vein label: 48
    imgSecVeinOverlay = MarkFG(img, imgSecVeinMask, [0 255 0], [0.5 0.5]); % secondary veins (green)
    imgTerVeinOverlay = MarkFG(imgSecVeinOverlay, imgTerVeinMask, [255 0 0], [0.5 0.5]); % tertiary veins (red)
    imgTerVeinOverlay = MarkFG(imgTerVeinOverlay, imgAbnormalVeinMask, [255 255 0], [0.5 0.5]); % abnormal veins (yellow)
    imgAllVeinOverlay = MarkFG(imgTerVeinOverlay, imgQuatVeinMask, [0 0 255], [0.5 0.5]); % vertical veins + quaternary veins (blue)
    if (DEBUG)
        figure(1); subplot(2,2,4); imshow(imgAllVeinOverlay); title('Color-labeled Veins');
    end
    
    % %-- Compute number of abnormal cases
    imgAbnormalVeinLabel = bwlabel(imgAbnormalVeinMask); % connected component labeling for abnormal veins
    numAbnormalCases = max(unique(imgAbnormalVeinLabel)); % number of abnormal cases
    
    % %-- Compute number of tertiary veins (including abnormal)
    % %-- And compute width/thickness for tertiary veins including abnormal veins
    numTerVein = 0;
    numPatches = 5; % number of mask patches to check vein number
    patchHeight = ceil(rowsOrig/numPatches); % patch height
    patchStartRowIdx = [1:patchHeight:rowsOrig, rowsOrig];
    imgTerVeinLabel = bwlabel(imgTerVeinMask); % connected component labeling for tertiary veins
    terVeinLabelStats = regionprops('table', logical(imgTerVeinLabel), 'Area', 'Extrema');
    terVeinLabelAreas = terVeinLabelStats.Area; % area for each tertiary vein
    terVeinLabelExtrema = terVeinLabelStats.Extrema; % extrema points for each tertiary vein
    for i = 1:max(unique(imgTerVeinLabel)) % process each tertiary vein
        %-- count number of veins
        currTerVeinMask = zeros(size(imgGrayNorm));
        currTerVeinMask(imgTerVeinLabel==i) = 1; % mask for the current tertiary vein
        currTerVeinNumVec = [];
        for j = 1:numPatches
            currTerVeinMaskPatch = currTerVeinMask(patchStartRowIdx(j):patchStartRowIdx(j+1)-1, :);
            currTerVeinMaskPatchLabel = bwlabel(currTerVeinMaskPatch);
            currTerVeinNumVec = [currTerVeinNumVec, max(unique(currTerVeinMaskPatchLabel))];
        end
        numCurrVein = max(1, mode(currTerVeinNumVec));
        numTerVein = numTerVein + numCurrVein;
        %-- compute vein thickness
        veinExtrema = terVeinLabelExtrema{i,1};
        veinHeight = max(veinExtrema(:,2)) - min(veinExtrema(:,2)); % vein height
        veinWidth = round(terVeinLabelAreas(i) / veinHeight); % vein width/thickness
        veinWidth = veinWidth / numCurrVein; % divided by number of veins
        vertVeinWidthList = [vertVeinWidthList veinWidth];
    end
    abnormalVeinLabelStats = regionprops('table', logical(imgAbnormalVeinLabel), 'Area', 'Extrema');
    abnormalVeinLabelAreas = abnormalVeinLabelStats.Area; % area for each abnormal vein
    abnormalVeinLabelExtrema = abnormalVeinLabelStats.Extrema; % extrema points for each abnormal vein
    for i = 1:numAbnormalCases % process each abnormal vein
        %-- count number of veins
        currVeinMask = zeros(size(imgGrayNorm));
        currVeinMask(imgAbnormalVeinLabel==i) = 1; % mask for the current abnormal vein
        currVeinNumVec = [];
        for j = 1:numPatches
            currVeinMaskPatch = currVeinMask(patchStartRowIdx(j):patchStartRowIdx(j+1)-1, :);
            currVeinMaskPatchLabel = bwlabel(currVeinMaskPatch);
            currVeinNumVec = [currVeinNumVec, max(unique(currVeinMaskPatchLabel))];
        end
        numCurrVein = max(1, mode(currVeinNumVec));
        numTerVein = numTerVein + numCurrVein;
        %-- compute vein thickness
        veinExtrema = abnormalVeinLabelExtrema{i,1};
        veinHeight = max(veinExtrema(:,2)) - min(veinExtrema(:,2)); % vein height
        veinWidth = round(abnormalVeinLabelAreas(i) / veinHeight); % vein width/thickness
        veinWidth = veinWidth / numCurrVein; % divided by number of veins
        vertVeinWidthList = [vertVeinWidthList veinWidth];
    end
    vertVeinWidthMean = mean(vertVeinWidthList); % average thickness of vertical veins
    
    % %-- Compute vein density
    vertVeinLabelStats = regionprops('table', logical(imgVertVeinLabel), 'Extrema');
    vertVeinLabelExtrema = vertVeinLabelStats.Extrema; % extrema points for each vertical vein
    numVertVein = length(vertVeinLabelExtrema); % number of vertical veins
    for i = 1:numVertVein
        veinExtrema = vertVeinLabelExtrema{i,1};
        if (i == 1)
            firstVeinIdxX = min(veinExtrema(:,1)); % X position of left edge of the 1st vein
        end
        if (i == numVertVein)
            lastVeinIdxX = max(veinExtrema(:,1)); % X position of right edge of the last vein
        end
    end
    veinDensity = (numSecVein + numTerVein) / (lastVeinIdxX - firstVeinIdxX); % vein density: # veins / dist(1st vein, last vein)
    roiWidth = lastVeinIdxX - firstVeinIdxX; % width for image ROI
    
    % %-- Compute stats for quaternary veins
    imgQuatVeinLabel = bwlabel(imgQuatVeinMask); % connected component labeling
    quatVeinLabelStats = regionprops('table', logical(imgQuatVeinLabel), 'Area');
    numQuatVein = length(quatVeinLabelStats.Area);
    
    % %-- Export stats to a csv file
    if (EXPORT_CSV)
        % %-- export vein thickness to csv file
        if exist(csvNameVeinWidth, 'file')
            fid1 = fopen(csvNameVeinWidth, 'a');
            fprintf(fid1, '%s,', imgRawName(1:end-4));
            dlmwrite(csvNameVeinWidth, vertVeinWidthList, '-append');
            fclose(fid1);
        else
            fid1 = fopen(csvNameVeinWidth, 'w');
            fprintf(fid1, '%s,', imgRawName(1:end-4));
            dlmwrite(csvNameVeinWidth, vertVeinWidthList, '-append');
            fclose(fid1);
        end
        % %-- export overall stats to csv file
        if exist(csvNameAllStats, 'file')
            fid1 = fopen(csvNameAllStats, 'a');
            fprintf(fid1, '%s\n', ' ');
            fprintf(fid1, '%s,%.2f,%.6f,%d,%.6f,%.6f,%d,%d,%d,%d,%.2f,%.6f,%.2f,%.6f,%.2f,%.6f', ...
                imgRawName(1:end-4), roiWidth, roiWidth/SCALE_FACTOR/1000, numSecVein+numTerVein, veinDensity, ...
                veinDensity*SCALE_FACTOR*1000, numSecVein, numTerVein, numQuatVein, numAbnormalCases, ...
                vertVeinWidthMean, vertVeinWidthMean/SCALE_FACTOR, veinDistMean, veinDistMean/SCALE_FACTOR, ...
                mean(quatVeinDistList), mean(quatVeinDistList)/SCALE_FACTOR);
            fclose(fid1);
        else
            fid1 = fopen(csvNameAllStats, 'w');
            csvHeaders = {'Image ID', 'ROI (px)', 'ROI (mm)', 'No. Long v. (2+3)', 'Long v. density (vein/px)', ...
                'Long v. density (vein/mm)', 'No. 2 v.', 'No. 3 v.', 'No. 4 v.', 'No. Irreg v.', ...
                'Mean v. width (px)', 'Mean v. width (um)', 'Mean interv. dist. (px)', 'Mean interv. dist. (um)', ...
                'Mean 4 v. interval (px)', 'Mean 4 v. interval (um)'};
            fprintf(fid1, '%s,', csvHeaders{1, 1:end});
            fprintf(fid1, '%s\n', ' ');
            fprintf(fid1, '%s,%.2f,%.6f,%d,%.6f,%.6f,%d,%d,%d,%d,%.2f,%.6f,%.2f,%.6f,%.2f,%.6f', ...
                imgRawName(1:end-4), roiWidth, roiWidth/SCALE_FACTOR/1000, numSecVein+numTerVein, veinDensity, ...
                veinDensity*SCALE_FACTOR*1000, numSecVein, numTerVein, numQuatVein, numAbnormalCases, ...
                vertVeinWidthMean, vertVeinWidthMean/SCALE_FACTOR, veinDistMean, veinDistMean/SCALE_FACTOR, ...
                mean(quatVeinDistList), mean(quatVeinDistList)/SCALE_FACTOR);
            fclose(fid1);
        end
    end
    
    % %-- Export vein segmentation mask and color-coded result images
    if (EXPORT_IMG)
        if ~exist(strcat(savePath, 'Label'), 'dir')
            mkdir(savePath, 'Label');
        end
        saveImgLabelName = strcat(savePath, '/Label/', imgRawName(1:end-4), '_label.png');
        imwrite(imgVeinLabel, saveImgLabelName);
        if ~exist(strcat(savePath, 'Viz'), 'dir')
            mkdir(savePath, 'Viz');
        end
        saveImgVizName = strcat(savePath, '/Viz/', imgRawName(1:end-4), '_viz.jpg');
        imwrite(imgAllVeinOverlay, saveImgVizName);
    end
    
end









