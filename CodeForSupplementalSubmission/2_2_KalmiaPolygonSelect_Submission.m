%% choose a file
% in the directory of files that you want to digitize
[FileName,PathName,FilterIndex] = uigetfile
%%
foo = strsplit(FileName, '.')

file = char(strcat(PathName, foo(1), '.mp4'));

vidObj = VideoReader(strcat(PathName, FileName));

% get info
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

video = read(vidObj,Inf);

imshow(video); 

% outline all of the pollen (exclude the anther, when it has no pollen)
[BW, xi, yi] = roipoly(video);


% look at non BW video to identify parts of flower
bar = FileName(1:15);
file = char(strcat(PathName, bar, '.avi'));

vidObj = VideoReader(file);

% get info
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;


video = read(vidObj,Inf);

imshow(video); 
% select the four parts of the flower for scaling
[xk,yk] = ginput(4)

close figure 1; 

xvals = vertcat(xi, xk); 
yvals = vertcat(yi, yk); 


writetable(table(xvals, yvals),strcat('~/Desktop/', bar, '_polygon.txt'))



%% Do same thing, except process a whole directory of _Tail videos

dirName = '/Users/callinswitzer/Dropbox/KalmiaProject/KalmiaExamples/KalmiaManualTrig_VidsToProcess/';              %# folder path
files = dir( fullfile(dirName,'*_Tail.mp4') );   %# list all *.xyz files
files = {files.name}'; 
cd(dirName);

files


for ii = 27:length(files)
    
    
    disp(char(strcat(num2str(ii), '_of_', num2str(length(files)))))
    % get info for filename
    foo = strsplit(char(files(ii)) , '.')
    

    % read in last frame of video
    vidObj = VideoReader(char(files(ii)));
    %vidHeight = vidObj.Height; % get info
    %vidWidth = vidObj.Width;
    video = read(vidObj,Inf);

    % display image
    %imshow(video,'InitialMagnification',130); 

    % outline all of the pollen , include anther then pollen is in it.
    [BW, xi, yi] = roipoly(video);

    % look at non Binary image to identify parts of flower
    FileName = char(files(ii));
    bar = FileName(1:15);
    file = char(strcat(bar, '.avi'));

    vidObj = VideoReader(file);
    

    % get info
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;


    video1 = read(vidObj,Inf);

    imshow( imadjust( video1,  [], [1; 0.0]), 'InitialMagnification',130); 
    % select the four parts of the flower for scaling
    
    [xk,yk] = ginput(9);

    close figure 1; 
   

    xvals = vertcat(xi, xk); 
    yvals = vertcat(yi, yk); 

    % save talbe
    writetable(table(xvals, yvals),strcat( bar, '_polygon.txt'));

end
    
