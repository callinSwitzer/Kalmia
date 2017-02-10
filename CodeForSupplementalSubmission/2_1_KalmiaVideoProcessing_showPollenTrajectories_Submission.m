%% Callin Switzer

%% 06/29/2016
% revisited 10 Aug 2016
% processing videos of kalmia to remove background and visualize
% the pollen trajectory
%%

% choose files
clear all
close all

% these videos are not included in supplemental data, 
% because they are too big (10GB)!
PathName = '/Users/callinswitzer/Dropbox/ExperSummer2016/Kalmia/KalmiaExamples/KalmiaManualTrig_VidsToProcess';

% choose file
cd(PathName);
%%
% Part 1 Setup
videos = dir2cell('*.avi'); % function was downloaded 
% videos(1) = [] % gets rid of ampcheck

charNum = 15; % number of characters to extract from title

%% read in each video one at a time preprocess them all
% make black and white and save as .mp4
% this version doesn't show the tail on the video
for ii =  1:length(videos)
    nums = ii; 
    vidName = char(videos(nums));
    mm = VideoReader(vidName);
    
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
        %imshow(bwimage);
        %drawnow;
        tmp = cmap(bwimage + 1, :);
        tmp = reshape(tmp, size(frames,1), size(frames,2), 3);
        writeVideo(video, tmp);
        kk
    end
    
    close(video)
end
%% write video to show tails
for ii =  1:length(videos)
    nums = ii; 
    vidName = char(videos(nums));
    mm = VideoReader(vidName);
    
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
    thresh = 0.95;
  
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

%         imshow(bwimage);
%         drawnow;
    %     tmp = cmap(bwimage + 1, :);
    %     tmp = reshape(bwimage, size(frames,1), size(frames,2), 3);
        bw_double = double(bwimage >= 1);
         writeVideo(video, bw_double);
        kk
    end

    
    close(video)
end


