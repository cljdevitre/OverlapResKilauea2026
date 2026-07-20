%% Charlie's basic script
% Containing some useful simple processing and plotting options

% Updated May 2025 for mtex 6.0.0 (using updated grain calculation code)

close all
clear all

% Start up mtex
addpath("/Users/cljd/Documents/MATLAB/mtex-6.0.0")
startup_mtex

%% Set up

% % crystal symmetry %% Turns out this is not necessary if you don't
% specify CS!! (Charlie and I realized this on 5/15/2025)

% CS = {... 
%   'notIndexed',...
%   crystalSymmetry('-1', [8.2 13 7.1], [93.44,116.21,90.23]*degree, 'X||a*', 'Z||c', 'mineral', 'Andesine', 'color', [0.53 0.81 0.98]),...
%   crystalSymmetry('12/m1', [8.6 13 7.2], [90,116.073,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Orthoclase', 'color', [0.56 0.74 0.56]),...
%   crystalSymmetry('6/m', [9.5 9.5 6.9], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Apatite', 'color', [0.85 0.65 0.13]),...
%   crystalSymmetry('-3m1', [4.9 4.9 5.5], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Quartz-new', 'color', [0.94 0.5 0.5]),...
%   crystalSymmetry('m-3m', [8.4 8.4 8.4], 'mineral', 'Magnetite', 'color', [0 0 0.55]),...
%   crystalSymmetry('12/m1', [6.6 8.7 7.4], [90,119.72,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Titanite', 'color', [0 0.39 0]),...
%   crystalSymmetry('12/m1', [9.9 18 5.3], [90,105.2,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Hornblende', 'color', [0.55 0 0]),...
%   crystalSymmetry('4/mmm', [6.6 6.6 5.9], 'mineral', 'Zircon', 'color', [0.58 0 0.83]),...
%   crystalSymmetry('-3', [5.1 5.1 14], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ilmenite', 'color', [0.5 0.5 0.5]),...
%   crystalSymmetry('12/m1', [5.4 9.3 10], [90,100.222,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Biotite', 'color', [0.85 0.75 0.85]),...
%   crystalSymmetry('12/m1', [5.4 9.3 10], [90,100.161,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Annite', 'color', [0.5 1 0.83]),...
%   crystalSymmetry('12/m1', [5.4 9.3 14], [90,96.282,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Clinochlore 1MIa', 'color', [0.24 0.7 0.44]),...
%   crystalSymmetry('-3m1', [5 5 17], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Calcite', 'color', [1 0.63 0.48]),...
%   crystalSymmetry('-1', [7.9 7.6 7.2], [76.55,104.38,66.87]*degree, 'X||a*', 'Z||c', 'mineral', 'Microcline', 'color', [0.27 0.51 0.71]),...
%   crystalSymmetry('12/m1', [9.7 9 5.3], [90,105.63,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside   CaMgSi2O6', 'color', [0.58 0.44 0.86]),...
%   crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77', 'color', [1 0.89 0.71]),...
%   crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.42 0.56 0.14]),...
%   crystalSymmetry('12/m1', [43 9.2 7.3], [90,91.6,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Antigorite superstructure', 'color', [0.53 0.81 0.98]),...
%   crystalSymmetry('12/m1', [5.3 9.2 7.3], [90,9.877,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Lizardite', 'color', [0.56 0.74 0.56]),...
%   crystalSymmetry('-1', [5.2 8.9 7.4], [91.928,105.044,89.792]*degree, 'X||a*', 'Z||c', 'mineral', 'Kaolinite 1A', 'color', [0.85 0.65 0.13])};
%

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Set path & file names
pname = ["/Users/cljd/pCloud Drive/WORK-GENERAL/POSTDOC-UCB/BERKELEY-VIBE/Documents/Projects/PrePostCaldera_Kil2025/Github_repo/Data_processing_notebooks/Data_reprocessing_Literature/Data_Lit/K23_EBSDmap"];

specific_name="/K23_wholemount_052225 (montaged map) - EBSD Data.ctf";
fname = pname + specific_name;

% This turns off the X Y labels on pole figures
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);    


%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');

 % Charlie's EBSD correction
    % rotate everything 180 degrees around x to get the correct map
    ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
    % then rotate only the orientations 180 degrees around z
    ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepXY');

%% Check that map looks correct
plot(ebsd)

%% Calculate grains
%Charlie's grain definition, 10 degrees for true grains
[grains, ebsd.grainId] = calcGrains(ebsd,'alpha',2.2,'angle',10*degree,'minPixel',10);

% This is from Penny's script, kuwahara 5 neighbor and smoorthing
ebsd(grains(grains.grainSize<2)).phase=0; %removes misindexed pixels 
F = KuwaharaFilter; 
F.numNeighbours = 5;
ebsd= smooth(ebsd('indexed'), F, 'fill',grains);

% Penny's grain definition

% This is used to pull out what counts as a separate grain 
% when plotting each grain separatly
subgrain_angle=0.3; 
%
grain_div_angle=3; 
color_angle=3;
len_filter=1; % filter out shorter grain boundaries than this um

[grains,ebsd.grainId]= calcGrains(ebsd,'alpha',2.2,'angle',[subgrain_angle*degree,grain_div_angle*degree],'minPixel',10); % Want this big, so it pulls out each grain separatly. 


%% This is to save data and make the grain plot to match with Raman
% We want to save the data now and export
% Sub only the forsterite and count grains
foGrains = grains('fo');

iter=[1,5,8];
for Z=1:length(iter)
    foGrains=smooth(foGrains, iter(Z));
end
foGrains=foGrains(foGrains.area >5000);

numGrains = numel(foGrains);
% Set the variables we want to export
grainID            = foGrains.id;
GOS                = foGrains.GOS ./ degree;
grainSize          = foGrains.grainSize;
area               = foGrains.area;
subBoundaryLength  = foGrains.subBoundaryLength;
equivRadius        = foGrains.equivalentRadius;

% Create empty mats to store the GB lengths.
Tot_length  = zeros(numGrains,1);
prop_length = zeros(numGrains,1);

% Calculate the GB lengths per grain
for k = 1:numGrains
    gb_C_m = foGrains(k).innerBoundary;
    Tot_length(k) = sum(gb_C_m.segLength);
    prop_length(k) = Tot_length(k) / sqrt(grainSize(k));
end

% Placeholder values
placeholderCols = zeros(size(grainID));
subgrain_angle_export  = subgrain_angle*(ones(size(grainID)));  
grain_div_angle_export = grain_div_angle*(ones(size(grainID)));  

Filename=repmat(string(specific_name), numel(grainID), 1);

% Create table
resultsTable = table(Filename,grainID, GOS, grainSize, area, subBoundaryLength, ...
    equivRadius, prop_length, Tot_length,...
    placeholderCols, placeholderCols, placeholderCols, ...
    placeholderCols, placeholderCols, ...
    subgrain_angle_export, grain_div_angle_export, ...
    'VariableNames', {'Filename','grainID',	'GOS' ,'Grain Size (pixels)', 'Grain Size (um2)', ...
    'subBoundaryLength (um)', 'equivalentRadius (um)','GB Length/Sqrt Size', 'Total GB length', ...
    'Tilt length',	'Twist length',	'Perc Tilt',	'Perc Twist', 'Perc unclassified', 'subgrain_angle', 'grain_div_angle'});


% Export to Excel
writetable(resultsTable, pname+'/grain_results.xlsx');

% Now we plot to match grains
plot(foGrains)

hold on
grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 8);
end
legend off
hold off
saveas(gcf,'Grains.pdf')

%% All options past this are optional. Run sections individually.
% This is set up for opx. To change the phase, replace 'fo' with the
% first few letters of your phase.

% You can call the crystal symmetry either by using ebsd('fo').CS or by
% specifying a symmetry like so:
csFo = crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77');
csFo = crystalSymmetry('mmm', [4.7560 10.2070 5.9800], 'mineral', 'Forsterite');

%% --- Basic maps ---

%% Nice band contrast map

plot(ebsd,ebsd.bc)
colormap gray 
mtexColorbar


%% Euler map

colorKey = BungeColorKey(ebsd('fo').CS);

plot(ebsd('fo'),colorKey.orientation2color(ebsd('fo').orientations),'micronbar','off')

%% IPF map

% make the colour key
ipfKey_fo = ipfColorKey(ebsd('fo').CS);
% specify reference direction (e.g., xvector, yvector, zvector)
ipfKey_fo.inversePoleFigureDirection = zvector; 

plot(ebsd('fo'),ipfKey_fo.orientation2color(ebsd('fo').orientations),'micronbar','off')

%% --- Pole figures ---
% For all pole figures: 'hkl' for poles to planes
%                       'uvw' for crystallographic directions

%                       'lower', 'upper' or 'complete' for hemispheres

%                       'noSymmetry' to remove antipodal symmetry

% Can choose to specify a set of planes or axes like so:
axes = Miller({1,0,0},{0,1,0},{0,0,1},'uvw',csFo);
planes = Miller({1,0,0},{0,1,0},{2,1,0},{1,1,1},'hkl',csFo);

%% Plain pole figures (point per pixel)

%plotPDF(ebsd('fo').orientations,Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','MarkerColor','k','MarkerSize',2)

plotPDF(ebsd('fo').orientations,axes,'lower','MarkerColor','k','MarkerSize',2)

%% Plain pole figures (point per grain)
plotPDF(grains('fo').meanOrientation,axes,'lower','MarkerColor','k','MarkerSize',2)

%% Axes coloured by ipfkey
plotPDF(ebsd('fo').orientations,ipfKey_fo.orientation2color(ebsd('fo').orientations),Miller({1,0,0},{0,1,0},{0,0,1},'uvw',csFo),'lower','MarkerSize',2)

%% Planes coloured by ipfkey
plotPDF(ebsd('fo').orientations,ipfKey_fo.orientation2color(ebsd('fo').orientations),Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','MarkerSize',2)

%% Smoothed pole figure

% Smoothed
plotPDF(ebsd('fo').orientations,Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','smooth','halfwidth',5*degree)
mtexColorMap LaboTeX
setColorRange('equal')
mtexColorbar('location','eastoutside','title','multiples of uniform distribution')
   

%% --- Internal distortion ---

%% Axis-angle colouring
% Good for visualising subtle distortion 
% Difficult to interpret the colours intuitively

% Set up colour key
AxisAngleKey = axisAngleColorKey(ebsd('fo'));
AxisAngleKey.oriRef = grains(ebsd('fo').grainId).meanOrientation;
AxisAngleKey.maxAngle = 3*degree;

color = AxisAngleKey.orientation2color(ebsd('fo').orientations);
plot(ebsd('fo'),color,'micronbar','off')

hold on
plot(grains.boundary,'linewidth',0.5)
hold off

%% Angular intragranular misorientation in degrees (GROD)

% Calculate misorientation of each pixel relative to the average
% orientation of its grain
grod = ebsd('fo').calcGROD(grains);

plot(ebsd('fo'),grod.angle./degree,'micronbar','off')

mtexColorbar('title','GROD angle in degrees')
setColorRange([0 3])
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',0.5)
legend off
hold off

%% Basic KAM map (to visualise subgrain boundaries)

% Tell mtex that your pixels are in a grid
ebsd = ebsd.gridify;

kam = ebsd.KAM / degree;

plot(ebsd,kam,'micronbar','off')
setColorRange([0,3]) % can fiddle with this or remove it
mtexColorbar
mtexColorMap LaboTeX

hold on
plot(grains.boundary,'lineWidth',0.5)

%% --- CRYSTALSHAPES ---

% To design a crystalShape for your mineral I recommend 
% midat.org to find likely habits
% https://www.smorf.nl/draw.php to design a shape (use 'crystallographic'
% setting)

% Some minerals (e.g. olivine) have shapes pre-loaded into mtex


% define the crystal shape of Forsterite and store it in the variable cS
cS = crystalShape.olivine(ebsd('Forsterite').CS);

%plot(cS,'colored')

%% Plot of crystalShapes of grains

plot(ebsd('fo'),'facecolor','lightgrey','micronbar','off')
hold on
plot(grains.boundary,'lineWidth',0.5)
plot(grains('fo'),0.7*cS,'linewidth',0.1,'colored')
hold off



