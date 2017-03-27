%% Set Image Path
% the folder where all the images from a single animal are kept
% you only need to run this if you want to work on another folder during a
% session
folder = uigetdir;
files = dir(folder);
files(1:2) = [];
files([files.isdir]) = [];
file_names = {files.name};
%% Create a structure to enter you data
% T becomes the structure array. 
% You can save this and reopen it using open at any time
% If you open a structure you do not need to create a new one
T = struct;
%% 1:
    % Load image 
    % brings up a list dialog that you can use to load a single image
    % splits channels into red, green, and blue
[s,ok] = listdlg('ListString',file_names, 'ListSize', [200 600]);
if ok
    img = imread(fullfile(folder, file_names{s}));
end

% add the file name, what image in file, genotype, and stain
T(end+1).name = file_names{s};
T(end).id = s;


%% Run for additional ROIS
    % adds in the information from the file without having to open the file
    % again
    % add the file name, what image in file, genotype, and stain
    
T(end+1).name = file_names{s};
T(end).id = s;

% Can alter anything in purple (string) to fit your experiment

%figure, imshow(img);
%% Create ROIs
% Select tissue to analyze 
% brings an image so you can trace the tissue of interest
% double click on line when you are done tracing
% split channels

green = img(:,:,2);
img_rot = imrotate(green,180);
green_adjust = imadjust(img_rot);
figure, imshow(green_adjust);
close all
%%
figure,
gra = imadjust(green_region,[0.01 0.2],[]);
imshow(green_final);
scaleBarWidth = 800;
    scaleBarHeight = 10;
    xPos = size(green_region,2)*0.9 - scaleBarWidth;
    yPos = size(green_region,1)*0.95 - scaleBarHeight;
    rectPosition = [xPos, yPos, scaleBarWidth, scaleBarHeight];
    hRect = rectangle('Position', rectPosition);
    set(hRect,'EdgeColor','w');
    set(hRect,'FaceColor','w');
%%
% ROI code
%img_adjusted = imadjust(img);
figure, 
h_img = imshow(green);
display('outline regions to remove from analysis')
h_roi_tissue = impoly(gca); % create a draggable, resizable polygon
position = wait(h_roi_tissue); % wait for a double click on polygon
display('tissue roi captured')
tissue = createMask(h_roi_tissue,h_img);
green_region = green;
green_region(~tissue) = 0;
figure, imshow(green_region), title('green')
%
% Identify the coordinates of the polygon (min and max)

x = position(:,1);
y = position(:,2);
min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);
width = abs(max_y - min_y+1);
height = abs(max_x - min_x+1);
%
% Crop the image to the coordinates of the polygon (imcrop)

green_region = imcrop(green_region,[min_x min_y width height]);

% Allows you to add the region you selected to array
prompt = 'What is the region? ';
region_name = input(prompt);
T(end).region = region_name;
% Calculates the area of the ROI selected
region = regionprops(tissue, 'basic');
area = region.Area; %in pixels
scale = 1./810; %(ENTER mm/pixel: 1/810 for 25x image)
tissue_thickness = 0.03; %(mm thick)
%tissue_area = area*scale*tissue_thickness;
tissue_area = ((sqrt(area)*scale)^2)*tissue_thickness;%(mm3)
T(end).img_area = tissue_area;
% closes all images

pause(1)
close all
%%
tissue_area = tissue_area/2;

%% 1:
    % Green channel (GFP) need to split large objects
    % Adjustments:
    %               n_big cutoff for too big = 70 pixels
    %               erode to find center of object to count = area of 2
n = green_region;
figure, imshow(n);
    [select_x, select_y] = ginput(1)
clear D_v
D_v = struct;
S = struct;
disk = 0
figure,
for i=1:25
    D_v(end+1).disk = disk;
    background = imopen(n,strel('disk',disk));
    n_adjusted = n - background;
    %figure, imshow(n_adjusted), title('after adjustment');
    %figure, imshow(n), title ('original');
% Threshold with set value
% Threshold green channel with set value
% creates a binary image
    thresh = 0.08;
    D_v(end).thresh = thresh
    n_thresh = im2bw(n_adjusted,thresh);
    green_final = n_thresh;
    %figure, imshowpair(n, n_thresh), title('threshold image')
% Determine connected object distribution
    objects = bwconncomp(n_thresh);
    num_objects = objects.NumObjects;
    numObjects = cellfun(@numel, objects.PixelIdxList);
    object_size = mean(numObjects);
    %plot_objects = histogram(numObjects);
    D_v(end).object_size = object_size;
    index_small=numObjects<5; 
    sum_small = sum(index_small);
    p_small = (sum_small/num_objects)*100;
    S(end+1).p_small = p_small;
    index_large =numObjects>50;
    sum_large = sum(index_large);
    p_large = (sum_large/num_objects)*100;
    S(end).p_large = p_large;

% Remove small objects
    n_elim = bwareaopen(n_thresh, 5);
    
    %figure, imshowpair(green_final, green_region);

% Remove large objects
    large = bwareaopen(n_elim, 50);
    green_final = n_elim - large;


% Count cells
    green_count = bwconncomp(green_final);
    green_num = green_count.NumObjects;
    green_density = green_num/tissue_area;

% Determine average cell size
    numPixels = cellfun(@numel,green_count.PixelIdxList);
   %figure,plot_size = histogram(numPixels);
    cell_size_g = mean(numPixels);


% Count cells
    % Density of Oligos
    %T(end).density_green = green_density;
    D_v(end).density = green_density;
    %T(end).count_green = green_num;
    D_v(end).count = green_num;
    %T(end).cell_size_g = cell_size_g;
    D_v(end).size = cell_size_g;
    image{i} = green_final
    h= subplot(5,5,i)
    sub_green = imcrop(green_region, [select_x, select_y, 100, 100]);
    subImage{i} = imcrop(image{i}, [select_x, select_y, 100, 100]);
    final_img = imshowpair(subImage{i},sub_green);
    scaleBarWidth = 80;
    scaleBarHeight = 5;
    xPos = size(sub_green,2)*0.9 - scaleBarWidth;
    yPos = size(sub_green,1)*0.9 - scaleBarHeight;
    rectPosition = [xPos, yPos, scaleBarWidth, scaleBarHeight];
    hRect = rectangle('Position', rectPosition);
    set(hRect,'EdgeColor','w');
    set(hRect,'FaceColor','w');
    title(['Disk: ',num2str(disk)])
    %pause(1)
    %close all
    disk = disk+1;
end

    
%% 1:
    % Green channel (GFP) need to split large objects
    % Adjustments:
    %               n_big cutoff for too big = 70 pixels
    %               erode to find center of object to count = area of 2
n = green_region;
figure, imshow(n);
    [select_x, select_y] = ginput(1)
clear T_v
T_v = struct;
thresh = 0.00;
S2 = struct;
figure, 
for i=1:10
    disk = 20;
    T_v(end+1).disk = disk;
    background = imopen(n,strel('disk',disk));
    n_adjusted = n - background;
    %figure, imshow(n_adjusted), title('after adjustment');
    %figure, imshow(n), title ('original');
% Threshold with set value
% Threshold green channel with set value
% creates a binary image
    
    T_v(end).thresh = thresh;
    n_thresh = im2bw(n_adjusted,thresh);
    green_final = n_thresh;
    %figure, imshowpair(n, n_thresh), title('threshold image')
% Determine connected object distribution
    objects = bwconncomp(n_thresh);
    num_objects = objects.NumObjects;
    numObjects = cellfun(@numel, objects.PixelIdxList);
    object_size = mean(numObjects);
    %plot_objects = histogram(numObjects);
    T_v(end).object_size = object_size;
 index_small=numObjects<5; 
    sum_small = sum(index_small);
    p_small = (sum_small/num_objects)*100;
    S(end+1).p_small = p_small;
    index_large =numObjects>50;
    sum_large = sum(index_large);
    p_large = (sum_large/num_objects)*100;
    S(end).p_large = p_large;

% Remove small objects
    n_elim = bwareaopen(n_thresh, 5);
    
    %figure, imshowpair(green_final, green_region);

% Remove large objects
    large = bwareaopen(n_elim, 50);
    green_final = n_elim - large;


% Count cells
    green_count = bwconncomp(green_final);
    green_num = green_count.NumObjects;
    green_density = green_num/tissue_area;

% Determine average cell size
    numPixels = cellfun(@numel,green_count.PixelIdxList);
   %figure,plot_size = histogram(numPixels);
    cell_size_g = mean(numPixels);

% Count cells
    % Density of Oligos
    %T(end).density_green = green_density;
    T_v(end).density = green_density;
    %T(end).count_green = green_num;
    T_v(end).count = green_num;
    %T(end).cell_size_g = cell_size_g;
    T_v(end).size = cell_size_g;
    image{i} = green_final

    subplot(2,5,i)
    sub_green = imcrop(green_region, [select_x, select_y, 100, 100]);
    subImage{i} = imcrop(image{i}, [select_x, select_y, 100, 100]);
    final_img = imshowpair(subImage{i},sub_green);
    scaleBarWidth = 80;
    scaleBarHeight = 10;
    xPos = size(sub_green,2)*0.9 - scaleBarWidth;
    yPos = size(sub_green,1)*0.9 - scaleBarHeight;
    rectPosition = [xPos, yPos, scaleBarWidth, scaleBarHeight];
    hRect = rectangle('Position', rectPosition);
    set(hRect,'EdgeColor','w');
    set(hRect,'FaceColor','w');
    title(['Threshold: ',num2str(thresh)])
    
    thresh = thresh+0.02
    %pause(1)
    %close all
end

%% Run with best values
disk = 9;
thresh = 0.14;
T(end).disk = disk;
T(end).thresh = thresh;

    %T(end).disk = disk;
    n = sub_green;
    background = imopen(n,strel('disk',disk));
    n_adjusted = n - background;
    figure, imshow(n_adjusted), title('after adjustment');
    figure, imshow(n), title ('original');
% Threshold with set value
% Threshold green channel with set value
% creates a binary image
   
    %T(end).thresh = thresh;
    n_thresh = im2bw(n_adjusted,thresh);
    figure, imshowpair(n, n_thresh), title('threshold image')
% Determine connected object distribution
    objects = bwconncomp(n_thresh);
    num_objects = objects.NumObjects;
    numObjects = cellfun(@numel, objects.PixelIdxList);
    object_size = mean(numObjects);
    figure, plot_objects = histogram(numObjects);
    T(end).object_size = object_size;
    index_small=numObjects<10; 
    sum_small = sum(index_small);
    p_small = (sum_small/num_objects)*100;
    T(end).p_small = p_small;
    index_large =numObjects>50;
    sum_large = sum(index_large);
    p_large = (sum_large/num_objects)*100;
    T(end).p_large = p_large;
    
% Remove small objects
    n_elim = bwareaopen(n_thresh, 5);
    
    figure, imshowpair(green_final, n);

% Remove large objects
    large = bwareaopen(n_elim, 70);
    green_final = n_elim - large;


% Count cells
    green_count = bwconncomp(green_final);
    green_num = green_count.NumObjects;
    green_density = green_num/tissue_area;

% Determine average cell size
    numPixels = cellfun(@numel,green_count.PixelIdxList);
    figure,plot_size = histogram(numPixels);
    cell_size_g = mean(numPixels);


% Count cells
    % Density of Oligos
    %T(end).density_green = green_density;
    T(end).density = green_density;
    %T(end).count_green = green_num;
    T(end).count = green_num;
    %T(end).cell_size_g = cell_size_g;
    T(end).size = cell_size_g;

    pause(1)
    %close all


%% Manual cell counts in region selected for threshold determination
figure, imshow(sub_green);
pts = getpts
cell = length(pts);
%% Save table
%%
S = struct2table(T);
writetable(S, fullfile(folder, '170320_WT_values_methods.csv'));