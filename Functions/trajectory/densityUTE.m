clearvars;
close all;

NPro = 51480; %Number of projections.
beta = 4; % For psf calculation, the same as in grid3; ratio of maximum expansion factor and current expansion factor; alpha_x=[1 2 4 ...], alpha_y=[1 2 4...].
Matrix = 128; %Matrix size of final images.
NPoints = 128; %Number of points along each projection, = Matrix/2 + NPShift + gradient ramp points.
fidpoints = 128; %Number of points in each FID ( = NPoints + numbe of zero-filling points).
NPShift = 0; %number of points shifted to compensate for signal build-up and ringing at the beginning of FID.



rampoints = NPoints - NPShift - Matrix/2; %Number of points sampled during gradient ramp.


RespMode = 'measured';          %How was the trajetory generated? Theoretical, or %Measured?
RealNPoints = NPoints - NPShift;




realpath = fullfile(path,num2str(fileNumber(i))); trajfilename = strcat('traj_',RespMode); DCFfilename = strcat('DCF_',RespMode,'.raw');
% load coordinates
fid = fopen(fullfile(realpath,trajfilename));
tmp = squeeze(fread(fid,inf, 'double')); fclose(fid);
tmp = reshape(tmp,[3,NPoints(i),NPro(i)]);
crds = tmp(:,1:RealNPoints,:);%cut ending points along one spoke; r = sqrt(crds(1,RealNPoints,:).^2 + crds(2,RealNPoints,:).^2 + crds(3,RealNPoints,:).^2);
crds = crds./max(r(:))/2;
crds = crds./max(r(:))/2/beta; % Normalized k-space to [-0.5 0.5]. disp(['generating DCF for ',realpath]);
'   SDC params:'
numIter
= 25; %Number of iternations to iteratively calculate sampling %density compensation.
effMtx
osf
verbose = 1;
= (RealNPoints - rampoints(i)) * 2; %Matrix size of images. = 2.1; %oversampling ratio of intermediate grid.
'   start SDC calc'
DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf); %Call the sampling density %compensation package from Zwart et al.
DCF = single(DCF); % float32
'   write output DCF'
tmp = DCF;
fid = fopen(fullfile(realpath,DCFfilename),'w');
148
fwrite(fid,tmp,'float32'); %Save sampling density compensation function %into a text file.
fclose(fid);
clear tmp DCF crds r trajfilename DCFfilename;
%Begin regridding and reconstructing. Regridding package by Zwart et al was used.
disp(['Regridding and Reconstructing ',realpath]); grid3(NPoints,NPro,rampoints,fidpoints,NPShift,RespMode,realpath);
%End