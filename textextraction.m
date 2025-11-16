clc;        
clear;      
close all;  
warning off;  

% 1.Load image
disp(' 1. Loading Image ');
try
    k = imread("C:\Users\mahad\OneDrive\Desktop\Input\WhatsApp Image 2025-11-16 at 16.53.36_2e3e9a97.jpg");
catch
    disp('Error: Image file not found. Please check the path.');
    return; % if image fails to load then stop
end

sam = k; % store original color image 
figure; 
imshow(k); 
title('1. Original Image');

% 2.Pre-process: contrast adjustment
disp(' 2. Adjusting Contrast ');
m_adj = imadjust(im2gray(k)); % Convert to grayscale and apply contrast adjustment

figure; 
imshow(m_adj); 
title('2. High Contrast Grayscale Image');

%3. Pre-process: skew correction
disp(' 3. Detecting & Correcting Skew ');

% Canny edge and Hough transform to find the dominant line angle
[r, c] = size(m_adj); 
BW = edge(m_adj, 'canny');
[H, T, R] = hough(BW);
P = houghpeaks(H, 1, 'threshold', ceil(0.3 * max(H(:)))); 
lines = houghlines(BW, T, R, P, 'FillGap', 0.5 * c, 'MinLength', 40);

% rotation angle 
if ~isempty(lines)
    angle = lines(1).theta; 
    
    %  Hough angle to rotation angle
    if angle < 0
        rot_angle = (90 - abs(angle));
    else
        rot_angle = (angle - 90);
    end
    
    disp(['Rotating by ', num2str(rot_angle), ' degrees.']);
    
    % Rotate both grayscale and color images
    deskewed_img = imrotate_white(m_adj, rot_angle);
    deskewed_color = imrotate_white(sam, rot_angle);
else
    % no lines found, assume no skew
    disp('No dominant lines found; proceeding without rotation.');
    deskewed_img = m_adj;
    deskewed_color = sam;
end

figure;
subplot(1, 2, 1); imshow(m_adj); title('Before Deskew');
subplot(1, 2, 2); imshow(deskewed_img); title('3. After Deskew');

%4. Detect MSER candidates
disp('4. Detecting MSER Candidates ');
% Blur and invert image for MSER
H = fspecial('gaussian', [3 3], 1);
gray_filt = imfilter(deskewed_img, H, 'replicate');
gray_inv = imcomplement(gray_filt); 

% Detect all candidate regions
[regions, mserConnComp] = detectMSERFeatures(gray_inv, ...
    'RegionAreaRange', [15 2000]);

%  properties for filtering
stats = regionprops(mserConnComp, 'BoundingBox', 'Area', 'PixelIdxList');

%  all found candidates
figure; 
imshow(deskewed_img); 
hold on;
plot(regions, 'showPixelList', true, 'showEllipses', false);
title('4. All MSER Candidates (Before Filtering)');
hold off;

%5. Filter 1: geometric pruning
disp('5. Pruning Candidates (Geometry) ');
% Remove candidates based on simple geometry (Area, Aspect Ratio, Extent)
minArea = 15;
maxArea = 2000;
minAspectRatio = 0.1;
maxAspectRatio = 1.5;
minExtent = 0.25;
maxExtent = 0.9;

keptIndices = []; % store indices of regions to keep

for i = 1:length(stats)
    bbox = stats(i).BoundingBox;
    if bbox(4) == 0, aspectRatio = 0; else, aspectRatio = bbox(3) / bbox(4); end
    bboxArea = bbox(3) * bbox(4);
    if bboxArea == 0, extent = 0; else, extent = stats(i).Area / bboxArea; end

    % Check if region is within all thresholds
    if stats(i).Area > minArea && ...
       stats(i).Area < maxArea && ...
       aspectRatio > minAspectRatio && ...
       aspectRatio < maxAspectRatio && ...
       extent > minExtent && ...
       extent < maxExtent
       
        keptIndices(end+1) = i;
    end
end

% only the kept regions
stats_geom_filtered = stats(keptIndices);
disp(['    ... ' num2str(length(stats_geom_filtered)) ' candidates remaining.']);

%6. Filter 2: stroke width variance pruning
disp(' 6. Pruning Candidates (Stroke Width) ');
% Filter based on stroke width consistency (text has low variance)
MAX_STROKE_VARIATION = 0.6; 
finalStats = struct([]);    

for i = 1:length(stats_geom_filtered)
    currentStat = stats_geom_filtered(i);
    originalIndex = keptIndices(i); 
    pixelIdxList = mserConnComp.PixelIdxList{originalIndex};
    
    %  mask of the single component
    componentMask = false(size(deskewed_img));
    componentMask(pixelIdxList) = true;
    
    % crop the mask to its bounding box
    bbox = round(currentStat.BoundingBox);
    bbox = [max(1, bbox(1)-1), max(1, bbox(2)-1), bbox(3)+2, bbox(4)+2]; 
    bbox(3) = min(bbox(3), size(componentMask, 2) - bbox(1) + 1);
    bbox(4) = min(bbox(4), size(componentMask, 1) - bbox(2) + 1);
    
    croppedMask = imcrop(componentMask, bbox);
    
    if isempty(croppedMask) || all(croppedMask(:) == 0)
        continue;
    end
    
    % Calculate distance transform to find stroke radii
    distMap = bwdist(~croppedMask);
    strokeRadii = distMap(croppedMask);
    meanRadius = mean(strokeRadii);
    stdRadius = std(strokeRadii);
    
    if meanRadius == 0 
        continue; 
    end
    
    % coefficient of variation
    coeffOfVariation = stdRadius / meanRadius;
    
    % keep if variation is low
    if coeffOfVariation < MAX_STROKE_VARIATION
        if isempty(finalStats)
            finalStats = currentStat;
        else
            finalStats(end+1) = currentStat;
        end
    end
end

disp(['    ... ' num2str(length(finalStats)) ' candidates remaining.']);

% 7. Display final detected text regions
disp('7. Visualizing Final Candidates ');
figure; 
imshow(deskewed_img); 
hold on;
title('7. Final Detected Text (After All Pruning)');

if isempty(finalStats)
    disp('No text regions found after filtering. Cannot proceed.');
    return; % Stop 
end

% green boxes around final character candidates
finalBBoxes = vertcat(finalStats.BoundingBox);
for i = 1:length(finalStats)
    rectangle('Position', finalBBoxes(i,:), 'EdgeColor', 'g', 'LineWidth', 1);
end
hold off;

%8. Group characters into text lines
disp(' 8. Grouping Characters into Lines ');
% individual characters into lines based on vertical length
y_centers = finalBBoxes(:,2) + finalBBoxes(:,4)/2;
avgHeight = mean(finalBBoxes(:,4));

% Cluster boxes based on Y-center distance
line_indices = clusterdata(y_centers, 'Cutoff', avgHeight * 0.5, 'Criterion', 'distance');
numLines = max(line_indices);
lineBBoxes = zeros(numLines, 4); 
[h, w] = size(deskewed_img); 
padding = 5; % Add padding to line boxes

for i = 1:numLines
    boxes_in_this_line = finalBBoxes(line_indices == i, :);
    
    if isempty(boxes_in_this_line)
        continue;
    end
    
    % Find the min/max boundaries for the entire line
    min_x = min(boxes_in_this_line(:,1));
    min_y = min(boxes_in_this_line(:,2));
    max_x = max(boxes_in_this_line(:,1) + boxes_in_this_line(:,3));
    max_y = max(boxes_in_this_line(:,2) + boxes_in_this_line(:,4));
    
    % Add padding
    min_x = max(1, min_x - padding);
    min_y = max(1, min_y - padding);
    max_x = min(w, max_x + padding);
    max_y = min(h, max_y + padding);
    
    % Store the final line box
    lineBBoxes(i,:) = [min_x, min_y, max_x - min_x, max_y - min_y];
end

% Display the grouped text lines
figure;
imshow(deskewed_img);
hold on;
title('8. Grouped Text Lines');
for i = 1:numLines
    if all(lineBBoxes(i,:) == 0) 
        continue;
    end
    rectangle('Position', lineBBoxes(i,:), 'EdgeColor', 'r', 'LineWidth', 2);
end
hold off;

% 9. Run ocr on the grouped lines
disp('9. Recognizing Text (OCR) ');
% Run OCR, constrained to the line boxes we found
ocrResults = ocr(deskewed_img, lineBBoxes);

% Display the recognized text in the command window
disp('The following text was recognized:');
disp('================================');
for i = 1:length(ocrResults)
    text = strtrim(ocrResults(i).Text);
    if ~isempty(text)
        disp(['Line ', num2str(i), ' detected: ''', text, '''']);
    end
end
disp('================================');

% 10. final summary
disp('10. Summary ');
% Print filtering statistics
disp(['    Initial MSER candidates: ', num2str(length(stats))]);
disp(['    After geometric pruning: ', num2str(length(stats_geom_filtered))]);
disp(['    After stroke width pruning: ', num2str(length(finalStats))]);
disp(['    Final text lines found: ', num2str(numLines)]);

%11. Display final recognized words
disp('--- 11. Generating Final Output Image ---');
figure;
imshow(deskewed_color); % Use the deskewed color image
hold on;
title('11. Final Recognized Words (Confidence > 50%)');

for i = 1:length(ocrResults)
    wordBBoxes = ocrResults(i).WordBoundingBoxes;
    words = ocrResults(i).Words;
    wordConfidences = ocrResults(i).WordConfidences;
    
    % ROBUSTNESS FIX 
    % Prevents error if OCR finds different # of boxes vs words
    numValidWords = min([size(wordBBoxes, 1), length(words), length(wordConfidences)]);
    
    for j = 1:numValidWords
        % Only show words with > 50% confidence
        if wordConfidences(j) > 0.5 
            
            % Draw the word's bounding box
            rectangle('Position', wordBBoxes(j,:), 'EdgeColor', 'cyan', 'LineWidth', 2);
            
            % Add the recognized word as a text label
            text_label = ['''', words{j}, ''''];
            text(wordBBoxes(j,1), wordBBoxes(j,2) - 10, text_label, ...
                 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end
hold off;
disp('--- Pipeline Complete ---');

%helper function(s) 

function rotated_image = imrotate_white(image, rot_angle_degree)
    % Rotates image and fills the new background areas with white
    
    % Determine fill value based on image type
    if islogical(image)
        fill_val = 0; % Black for logical
    elseif isinteger(image)
        fill_val = 255; % White for uint8
    else
        fill_val = 1.0; % White for double/single
    end
    
    tform = affine2d([cosd(rot_angle_degree)   -sind(rot_angle_degree)    0; ...
                      sind(rot_angle_degree)    cosd(rot_angle_degree)    0; ...
                      0                         0                         1]);

    % Crop the rotated image to the original size
    outputView = imref2d(size(image));
    
    % Rotate the image
    rotated_image = imwarp(image, tform, 'OutputView', outputView, ...
                           'Interp', 'cubic', 'FillValues', fill_val);
end