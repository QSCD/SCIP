inputImage = 'F:/Teststack.tif';
[~,fileName,~]=fileparts(inputImage);

folder=sprintf('%s\\Results_%s',fileparts(mfilename('fullpath')),fileName);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
addpath(sprintf('%s\\Fiji.app\\scripts',fileparts(mfilename('fullpath'))))
javaaddpath (sprintf('%s\\Fiji.app\\java\\mij.jar',fileparts(mfilename('fullpath'))))
javaaddpath (sprintf('%s\\Fiji.app\\java\\ij.jar',fileparts(mfilename('fullpath'))))
import java.util.HashMap
import ij.*
import fiji.plugin.trackmate.*
 Miji(false);
%% load image stack
disp('Loading image...');
 I = tiffread2(inputImage);
 Im = [];
for i = 1: numel(I)
    Im(:,:,i) = I(i).data;
end

 xRes = I(1).x_resolution; %resolution of the image
xRes = (xRes(2)/xRes(1));
%% detect spots on projection
disp('Detecting spots (Round 1)...');
spots = spotDetection(inputImage,xRes,folder);
 
%% filter spots and determine z coodinate
disp('Determine Z coordinate (Round 1)...');
s = getZtoXY(spots,Im);
disp('Filtering spots (Round 1)...');
[s,p]= excludeDistantCells(s,3);
disp('Writing out cropped image...');
outputIm = (sprintf('%s\\stackCropped.tif',folder));
writeOutCroppedImage(Im,p,outputIm);

%% redo analysis
disp('Loading cropped image...');
I = tiffread2(outputIm);
 Im = [];
for i = 1: numel(I)
    Im(:,:,i) = I(i).data;
end
%% detect spots on projection
disp('Detecting spots (Round 2)...');
 spots = spotDetection(outputIm,1,folder);
 
%% filter spots and determine z coodinate
disp('Determine Z coordinate (Round 2)...');
s = getZtoXY(spots,Im);
disp('Filtering Spots (Round 2)...');
s= excludeDistantCells(s,4);
h1=figure;set(gcf,'Visible', 'off');  Im = imread(sprintf('%s\\MAX_proj.tif',folder)); warning 'off'; imshow(Im);hold on;
scatter([s.x], [s.y],'go');
saveas(h1,(sprintf('%s\\MAX_proj.tif',folder)));
xlswrite((sprintf('%s\\coordinates.xlsx',folder)),{'X','Y','Z'},1,'A');
xlswrite((sprintf('%s\\coordinates.xlsx',folder)),[[s.x]',[s.y]',[s.z]'],1,'A2');

function spots=spotDetection(inputIm,xRes,folder)
import ij.*
imp = ImagePlus(inputIm);
 imp.show()
MIJ.run('Z Project...', 'projection=[Max Intensity]');
maxfile = (sprintf('%s\\MAX_proj.tif',folder));
MIJ.run('Save', sprintf('Tiff..., path=[%s]',maxfile));
MIJ.closeAllWindows()
imp = ImagePlus(maxfile);    

%----------------------------
% Create the model object now
%----------------------------
    
% Some of the parameters we configure below need to have
% a reference to the model at creation. So we create an
% empty model now.
model = fiji.plugin.trackmate.Model();
    
% Send all messages to ImageJ log window.
model.setLogger(fiji.plugin.trackmate.Logger.IJ_LOGGER)
       
%------------------------
% Prepare settings object
%------------------------
       
settings = fiji.plugin.trackmate.Settings();
settings.setFrom(imp)
       
% Configure detector - We use a java map
settings.detectorFactory = fiji.plugin.trackmate.detection.LogDetectorFactory();
map = java.util.HashMap();
map.put('DO_SUBPIXEL_LOCALIZATION', true);
map.put('RADIUS', 8*xRes);
map.put('TARGET_CHANNEL', 1);
map.put('THRESHOLD', 0.4);
map.put('DO_MEDIAN_FILTERING', false);
settings.detectorSettings = map;
    
% Configure spot filters - Classical filter on quality
filter1 = fiji.plugin.trackmate.features.FeatureFilter('QUALITY', 0.0, true);
settings.addSpotFilter(filter1)
     
% Configure tracker - We want to allow splits and fusions
settings.trackerFactory  = fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory();
settings.trackerSettings = fiji.plugin.trackmate.tracking.LAPUtils.getDefaultLAPSettingsMap(); % almost good enough
% settings.trackerSettings.put('ALLOW_TRACK_SPLITTING', false);
% settings.trackerSettings.put('ALLOW_TRACK_MERGING', false);
    
% Configure track analyzers - Later on we want to filter out tracks 
% based on their displacement, so we need to state that we want 
% track displacement to be calculated. By default, out of the GUI, 
% not features are calculated. 
    
% The displacement feature is provided by the TrackDurationAnalyzer.
settings.addTrackAnalyzer(fiji.plugin.trackmate.features.track.TrackDurationAnalyzer())
    
% Configure track filters - We want to get rid of the two immobile spots at 
% the bottom right of the image. Track displacement must be above 10 pixels.
filter2 = fiji.plugin.trackmate.features.FeatureFilter('TRACK_DISPLACEMENT', 2.0, false);
settings.addTrackFilter(filter2)
    
    
%-------------------
% Instantiate plugin
%-------------------
    
trackmate = fiji.plugin.trackmate.TrackMate(model, settings);
       
%--------
% Process
%--------
    
ok = trackmate.checkInput();
if ~ok
    display(trackmate.getErrorMessage())
end
 
ok = trackmate.process();
MIJ.run('Close')
if ~ok
    display(trackmate.getErrorMessage())
end
% get all the spots and their pixel values

it=model.getSpots.iterator(0);
it2=model.getSpots.iterator(0);
spots=struct;
while (it.hasNext())
x=double(it.next().getFeature('POSITION_X'));
y=double(it2.next().getFeature('POSITION_Y'));
spots(end+1).x=x/xRes;
spots(end).y=y/xRes;
end
spots(1) =[];
end


function  [sreturn] = getZtoXY(s,Im)
%returns z values to x,y coordianted spots (s) on the 3d image Im
cellRadiusXY = 2.5;
lowInt = true(numel(s),1);
for i = 1: numel(s)
   x = round(s(i).x)+1;
   y = round(s(i).y)+1;
   % !! switched x and y to compensate switched I.data
   ztube = (Im(max(y-round(cellRadiusXY/2),1):min(y+round(cellRadiusXY/2), size(Im,2)),max(x-round(cellRadiusXY/2), 1):min(x+round(cellRadiusXY/2), size(Im,1)),:));% to prevent that the tube is outside the image
   ground = min(mean(mean(ztube)));
   zs =[zeros(10,1); squeeze(mean(mean(ztube)))-ground; zeros(10,1)];
   
    % fit gaussian to the intensity profile in z
   f  = fit([1:numel(zs)]',zs, 'gauss1');

   if (max(zs) > 5 && (4*f.c1 > 5 || abs(f.b1-10 -size(Im,3)) < 3 || abs(f.b1-10) < 3 ))
       % compensate the zeros set at the beginning of the z array
       if  4*f.c1 > 40
           [~,id]=max(zs);
           s(i).z = id-10;
       else
       mid = min(numel(zs)-10,round(f.b1));
       [~,maxIDAroundMid] =  max(zs(max(mid-5,1):mid+5));
       s(i).z =mid-6+maxIDAroundMid-10; 
       end
       % the mean of gauss fit could be outside the z axis
       [~,maid] = max(zs);
       if(s(i).z < 1)
           s(i).z = maid-10;
       elseif (s(i).z > size(Im,3))
           s(i).z = maid-10;
       end
       % cell diameter = 10mum * mum/pxl
       ys = max(squeeze(Im(max(x-round(2*cellRadiusXY),1):min(x+round(2*cellRadiusXY),size(Im,1)),y,max(1,round(s(i).z)-5):min(round(s(i).z)+5,size(Im,3))))')';
       ys = [zeros(10,1); ys; zeros(10,1)];
       xs = max(squeeze(Im(x,max(y-round(2*cellRadiusXY),1):min(y+round(2*cellRadiusXY),size(Im,2)),max(1,round(s(i).z)-5):min(round(s(i).z)+5,size(Im,3))))');
       xs = [zeros(10,1); xs'; zeros(10,1)];
       fx = fit([1:numel(xs)]',xs, 'gauss1');
       fy = fit([1:numel(ys)]',ys, 'gauss1');

       % exclude too small (in x,y dim) cells 
       if (4*fx.c1 < cellRadiusXY && 4*fy.c1 < cellRadiusXY )
           lowInt(i) = false;
       end
   else
     lowInt(i) = false;
     %i;
   end
   
end
% correct for cells beeing too close to each other
s=s(lowInt);
uniq = true(numel(s),1);
  for i = 1: numel(s)
      for j = i+1: numel(s)
          if(sqrt((s(i).x-s(j).x)^2+(s(i).y-s(j).y)^2+(s(i).z-s(j).z)^2) < 1.0*cellRadiusXY && uniq(j))
              % keep higher intensity cell
              if(Im(round(s(i).y)+1,round(s(i).x)+1,round(s(i).z)) < Im(round(s(j).y)+1,round(s(j).x)+1,round(s(j).z)))
                  uniq(i) = false;
              else
                  uniq(j) = false;
              end
              
          end
          
      end
  end
  sreturn = s(uniq);
end

function [sreturn,p]= excludeDistantCells(s,dimension)

sreturn = s;

%% EM to exclude outliers
p=polyfitn2([[s.x]',[s.y]'],[s.z]',dimension);
while abs(0-max(p.RMSE.*p.Sign)) > 16 ||  abs(0-mean(p.RMSE.*p.Sign)) > 0.1
    [~,idx] = max(p.RMSE);
    sreturn(idx) = [];
    p=polyfitn2([[sreturn.x]',[sreturn.y]'],[sreturn.z]',dimension);
end
end

function writeOutCroppedImage(Im,p,outputFileName)
%% set values too far from plane to zero/bg value
imageSize = size(Im,1);
maxDist2plane = 5*5;
syms X1 X2;
[xg,yg]=meshgrid(1:1:imageSize);
zg = polyvaln(p,[xg(:),yg(:)]);
% background = mean of all 8 corners of the image
bg = round(mean ([Im(1,1,1), Im(1,1,size(Im,3)), Im(1, size(Im,2),1),...
    Im(1, size(Im,2), size(Im,3)), Im(size(Im,1), 1,1), Im(size(Im,1),1,size(Im,3)), ...
    Im(size(Im,1) , size(Im,2), 1), Im(size(Im,1), size(Im,2), size(Im,3))]));
ImOut = Im;
for i=1:size(Im,1)
    for j=1:size(Im,2)        
        z = round(zg((i-1)*imageSize + j));
        zmin = round(z - maxDist2plane);
        zmax = round(z + maxDist2plane);
        % !! switched x,y see getZtoXY.m line 11
        if (zmin >  size(Im,3) || zmax < 1)
            ImOut(j,i, 1:size(Im,3)) = bg;
        else
            ImOut(j,i, 1:max(1,zmin)) = bg;
            ImOut(j,i, min(size(Im,3), zmax): size(Im,3)) = bg;
        end
    end
    
    
end
ImOut = ImOut/255;
%write out tiff file

imwrite(ImOut(:, :, 1), outputFileName,'Compression', 'none');
for i=2:length(ImOut(1, 1, :))
   imwrite(ImOut(:, :, i), outputFileName, 'WriteMode', 'append','Compression', 'none');
end
end