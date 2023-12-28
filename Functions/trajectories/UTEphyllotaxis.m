%% test

clc;
clearvars;
close all;

%% 3D filling from scanner

data = load('/Users/Gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps example data/Retrospective/UTE/873/16621/lut_20130.txt');
coor = reshape(data,[3,length(data)/3]);

rr = 1:5000;

X1 = coor(1,rr);
Y1 = coor(2,rr);
Z1 = coor(3,rr);

scatter3(X1,Y1,Z1)

%% 3D golden means

phi1 = 0.4656;
phi2 = 0.6823;

Xm = [];
Ym = [];
Zm = [];
r = 32767;

for m = 1:20130

    a = pi*mod(m*phi2,1);
    b = 2*pi*mod(m*phi1,1);

    [Xm(m),Ym(m),Zm(m)] = sph2cart(a,b,32767);

end

scatter3(Xm,Ym,Zm)




%% Phyllotaxis 

tiny_golden_angles = [111.24611, 68.75388, 49.75077, 38.97762, 32.03967, 27.19840, 23.62814, 20.88643, 18.71484, 16.95229];

Xp = [];
Yp = [];
Zp = [];
rp = 32767;
pmax = 80; % number of phylotaxes
imax = 183*3;

cnt = 0;

for p = pmax:-1:0
        cnt = cnt + 1;
        Zp(cnt) = p/pmax; 
        phi = sqrt(2*pi*pmax/2)*asin(Zp(cnt));
        Yp(cnt) = sin(phi)*sqrt(1-Zp(cnt)^2);
        Xp(cnt) = cos(phi)*sqrt(1-Zp(cnt)^2);
end
for p = 1:pmax
        cnt = cnt + 1;
        Zp(cnt) = -p/pmax; 
        phi = sqrt(2*pi*pmax/2)*asin(Zp(cnt));
        Yp(cnt) = sin(phi)*sqrt(1-Zp(cnt)^2);
        Xp(cnt) = cos(phi)*sqrt(1-Zp(cnt)^2);
end

nrpoints = cnt;

for i = 1:imax
 
     phi1 = (i-1) * tiny_golden_angles(1);
     phi2 = (i-i) * tiny_golden_angles(1);

     for j = 1:nrpoints
 
        r = [Xp(j),Yp(j),Zp(j)]*roty(phi1,'deg')*rotx(phi1,'deg');
      
        Xp(nrpoints*(i-1)+j) = r(1);
        Yp(nrpoints*(i-1)+j) = r(2);
        Zp(nrpoints*(i-1)+j) = r(3);

     end

 
end

scatter3(Xp,Yp,Zp)




%% Normalize


Xp = round(Xp*32767/max(Xp(:)));
Yp = round(Yp*32767/max(Yp(:)));
Zp = round(Zp*32767/max(Zp(:)));




%% Show point by point

figure(1);
plot1 = scatter3(Xp(1),Yp(1),Zp(1),'s');

r = 1000:2000;

for i=r
    plot1.XData = Xp(r(1):i);
    plot1.YData = Yp(r(1):i);
    plot1.ZData = Zp(r(1):i);
    pause(0.0000000000001);
end




%% Save as file

pnts = [Xp;Yp;Zp];
pnts = reshape(pnts,[3*length(Xp),1]);

filename = 'lut_32000ph.txt';
fileID = fopen(filename,'w');

maxpoints = 32000;

for i = 1:maxpoints*3
   
    fprintf(fileID,num2str(pnts(i,1),'%6.f'));
    fprintf(fileID,'\n');
    
end

fclose(fileID);




