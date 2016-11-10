%% Callin Switzer

%% 06/29/2016
% processing videos of kalmia to remove background and visualize
% the pollen trajectory
%%

% choose file
clear all
close all

[FileName,PathName,FilterIndex] = uigetfile('*.avi')

% choose file
cd(PathName);
%%
% Part 1 Setup
videos = dir2cell('*.avi'); % function was downloaded 
% videos(1) = [] % gets rid of ampcheck

charNum = 15; % number of characters to extract from title

%% read in each video one at a time preprocess them all
% make black and white and save as .mp4
for ii =  1:3%length(videos)
    nums = ii; 
    vidName = char(videos(nums));
    mm = VideoReader(vidName);
    
    % Load in an image, display, and select bottom of stamen
    
    im = read(mm, 1);
    
    % info on the image
    % whos im
    
    % remove the color channel
    figure
    imagesc(im) % images in matlab are y, x
    
    % Get some info on the video
    NumberOfFrames  = mm.NumberOfFrames;
    Width           = mm.Width;
    Height          = mm.Height;
    % closing figure seems to help
    % For speed, preallocate the array
    frames = uint8(zeros(Height, Width, NumberOfFrames));
    
    % Load in the frames into the preallocated array
    
    for kk=1:NumberOfFrames
        tmp = read(mm,kk); % load in the kk'th image
        %imagesc(tmp);
        frames(:,:,kk) = tmp();    % save the reduced x,y portion of the image
        %imagesc(tmp);
        kk
    end
    
    display('frames loaded -- subtracting background')
    bkg = median(frames, 3);
    
    %imagesc(bkg);
    imwrite(bkg, 'bkg.tif');
    
    display('backgound subtracted')
    % thresholded + mediand + mean light intensity Pollen video
    lightintensity = squeeze(sum(sum(frames)))/(size(frames,1)*size(frames,2));
    % colormap gray
    % figure(1)
    thresh = 0.9;
  
    fyle = [strcat(vidName(1:charNum),'_noTail')];
    
    video = VideoWriter(fyle,'MPEG-4');
    video.FrameRate = 30;
    video.Quality = 100;
    
    
    open(video);
    
    cmap = gray(2);
    display('saving black and white video')
    for kk=1: NumberOfFrames
        normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
        normalizedbkg = double(bkg)/mean(bkg(:));
        bwimage = ((normalizedimage./normalizedbkg < thresh));
        imshow(bwimage);
        drawnow;
        tmp = cmap(bwimage + 1, :);
        tmp = reshape(tmp, size(frames,1), size(frames,2), 3);
        writeVideo(video, tmp);
        kk
    end
    
    close(video)
end
%% write video to show tails
for ii =  28%length(videos)
    nums = ii; 
    vidName = char(videos(nums));
    mm = VideoReader(vidName);
    
    % Load in an image, display, and select bottom of stamen
    
    im = read(mm, 1);
    imshow(im); 
    pause(3); 
    close all
    
    % info on the image
    % whos im
    
    % Get some info on the video
    NumberOfFrames  = mm.NumberOfFrames;
    Width           = mm.Width;
    Height          = mm.Height;
    % closing figure seems to help
    % For speed, preallocate the array
    frames = uint8(zeros(Height, Width, NumberOfFrames));
    
    % Load in the frames into the preallocated array
    
    for kk=1:NumberOfFrames
        tmp = read(mm,kk); % load in the kk'th image
        %imagesc(tmp);
        frames(:,:,kk) = tmp();    % save the reduced x,y portion of the image
        %imagesc(tmp);
        kk
    end
    
    display('frames loaded -- subtracting background')
    bkg = median(frames, 3);
    
    %imagesc(bkg);
    imwrite(bkg, 'bkg.tif');
    
    display('backgound subtracted')
    % thresholded + mediand + mean light intensity Pollen video
    lightintensity = squeeze(sum(sum(frames)))/(size(frames,1)*size(frames,2));
    % colormap gray
    % figure(1)
    thresh = 0.9;
  
    fyle = [strcat(vidName(1:charNum),'_Tail')];
    
    video = VideoWriter(fyle,'MPEG-4');
    video.FrameRate = 30;
    video.Quality = 100;
    
    
    open(video);
    
    cmap = gray(2);
    display('saving black and white video')
    for kk=1: NumberOfFrames
        normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
        normalizedbkg = double(bkg)/mean(bkg(:));
        if kk == 1
            bwimage = ((normalizedimage./normalizedbkg < thresh));
        else bwimage = bwimage + ((normalizedimage./normalizedbkg < thresh));
        end
        % remove digital noise by dialating and eroding
        seD = strel('disk',1);
        BW1 = imerode(bwimage,seD);
        bwimage = imdilate(BW1, seD); 

        imshow(bwimage);
        drawnow;
    %     tmp = cmap(bwimage + 1, :);
    %     tmp = reshape(bwimage, size(frames,1), size(frames,2), 3);
        bw_double = double(bwimage >= 1);
         writeVideo(video, bw_double);
        kk
    end

    
    close(video)
end


%% use regionprops to track pollen

close all

thresh = 0.9;

fyle = [vidName(1:charNum)];

video = VideoWriter(fyle,'MPEG-4');
video.FrameRate = 30;
video.Quality = 100;


%open(video);
clearvars centroids_saved
cmap = gray(2);
display('saving black and white video')
for kk=500:700%300:400%: NumberOfFrames
    normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
    normalizedbkg = double(bkg)/mean(bkg(:));

    bwimage = ((normalizedimage./normalizedbkg < thresh));
  
    imshow(bwimage);
    hold on; 
    drawnow;
    
    % detect blobs that are pollen
   
    % remove border blobs
    BWnobord = imclearborder(bwimage, 4);
    %imshow(BWnobord), title('cleared border image');

    % erode images
    seD = strel('square',2);
    BWfinal = imerode(BWnobord,seD);
    %imshow(BWfinal); 

    % dilate
    BWfinal = imdilate(BWfinal,seD);
    %imshow(BWfinal);
   
    
    
    %
    %
    bw = BWfinal >0;


    % get region properties
    s = regionprops(bw,'Centroid','MajorAxisLength',...
        'MinorAxisLength', 'Eccentricity', 'Area');
    centroids = cat(1, s.Centroid);

    imshow(bw)
    hold on
    plot(centroids(:,1),centroids(:,2), 'b*')
    hold off

    % manually set cutoff points, based on region props
    % histogram([s.Eccentricity])
    bigEnough = [s.MajorAxisLength] > 3 & [s.MajorAxisLength] < 50 ...
        & [s.MinorAxisLength] > 3 ;
    %show centroids
    
%     imshow(bw)
%     hold on
%     plot(centroids(bigEnough,1),centroids(bigEnough,2), 'r*')
%     hold off

    % keep specific regions (single pollen grains)
    labeledImage = bwlabel(bw, 8);
    %imshow(labeledImage);


    %imshow(labeledImage, []);  % Show the gray scale image.
    %title('Labeled Image, from bwlabel()');

    keeperIndexes = find(bigEnough);
    keeperBlobsImage = ismember(labeledImage, keeperIndexes); 

    newLabeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
    % Now we're done.  We have a labeled image of blobs that meet our specified criteria.
    imshow(newLabeledImage);
    title('"Keeper" blobs');
    drawnow; 
    %hold on
    %plot(centroids(bigEnough,1),centroids(bigEnough,2), 'r*');
    %hold off

    % save centroids
    if exist('centroids_saved', 'var') == 1; 
        centroids_saved = vertcat(centroids_saved, centroids(bigEnough, 1:2));
    else centroids_saved =  centroids(bigEnough, 1:2);
    end
   
%     tmp = cmap(bwimage + 1, :);
%     tmp = reshape(bwimage, size(frames,1), size(frames,2), 3);
    %bw_double = double(bwimage >= 1);
     %writeVideo(video, bw_double);
     %imshow(bw_double);
    kk
end

hold on
plot(centroids_saved(:,1), centroids_saved(:,2), 'b*');
hold off
 
 
 %close(video)






%% Use regionprops to track pollen tragectories
thresh = 0.9;

fyle = [vidName(1:charNum)];

video = VideoWriter(fyle,'MPEG-4');
video.FrameRate = 30;
video.Quality = 100;


open(video);

cmap = gray(2);
display('saving black and white video')
for kk=1: NumberOfFrames
    normalizedimage = double(frames(:,:,kk))/double(lightintensity(kk));
    normalizedbkg = double(bkg)/mean(bkg(:));
    if kk == 1
        bwimage = ((normalizedimage./normalizedbkg < thresh));
    else bwimage = bwimage + ((normalizedimage./normalizedbkg < thresh));
    end
    imshow(bwimage);
    drawnow;
    
    % detect edges
    % http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
    [~, threshold] = edge(I, 'sobel');
    fudgeFactor = 0.5;
    BWs = edge(I,'sobel', threshold * fudgeFactor);
    imshow(BWs), title('binary gradient mask');

    se90 = strel('disk', 10);

    BWsdil = imdilate(BWs, [se90]);
    imshow(BWsdil), title('dilated gradient mask'); 


    BWdfill = imfill(BWsdil, 'holes');
    imshow(BWdfill);
    title('binary image with filled holes');

    % remove border blobs
    BWnobord = imclearborder(BWdfill, 4);
    imshow(BWnobord), title('cleared border image');

    % erode images
    seD = strel('disk',10);
    BWfinal = imerode(BWnobord,seD);
    %BWfinal = imerode(BWfinal,seD);
    %BWfinal = imerode(BWfinal,seD);
    imshow(BWfinal), title('segmented image');


    subplot(1,2,1); 
    imshow(flr); 
    subplot(1,2,2); 
    imshow(BWfinal); 

%     tmp = cmap(bwimage + 1, :);
%     tmp = reshape(bwimage, size(frames,1), size(frames,2), 3);
    bw_double = double(bwimage >= 1);
     writeVideo(video, bw_double);
    kk
end

close(video)







%%
