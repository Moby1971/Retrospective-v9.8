% ----------------------------
%  Pseudo radial k-space trajectory
%  For MR Solutions custom 3D k-space (exLUT)
%
%  Gustav Strijkers
%  Oct 2023
%  
% ----------------------------

%% clear all

clc;
clearvars;
close all force;


%% Initialization

outputdir = pwd;

dimy = 80;              % k-space y dimension (no_views)
dimz = 80;              % k-space z dimension (no_views_2)
order = 1;              % 1 = back and forth, 2 = one direction
angle_nr = 10;          % golden angle number (see list below)
display = false;        % show result true / false

tiny_golden_angles = [111.24611, 68.75388, 49.75077, 38.97762, 32.03967, 27.19840, 23.62814, 20.88643, 18.71484, 16.95229];


%% fill the list

ry = round(dimy/2)-0.5;
rz = round(dimz/2)-0.5;
step = 1/max([dimy, dimz]);
radius = max([dimy, dimz]);
number_of_spokes = max([dimz,dimy])*2;
angle = 0;
kspacelist=[];

for ns = 1:number_of_spokes
    
    nr=1;
    clear c;
    
    if order == 1
        
        if rem(ns,2)
            
            for i=-1:step:1-step
                
                c(nr,:) = [floor(ry * i * cos(angle*pi/180)),floor(rz * i * sin(angle*pi/180))]; %#ok<*SAGROW> 
                nr = nr + 1;
                
            end
            
        else
            
            for i=1-step:-step:-1
                
                c(nr,:) = [floor(ry * i * cos(angle*pi/180)),floor(rz * i * sin(angle*pi/180))];
                nr = nr + 1;
                
            end
            
        end
        
    else
        
        for i=-1:step:1-step
            
            c(nr,:) = [floor(ry * i * cos(angle*pi/180)),floor(rz * i * sin(angle*pi/180))];
            nr = nr + 1;
            
        end
        
    end
    
    c = unique(c,'Rows','Stable');
    
    kspacelist = [kspacelist;c]; %#ok<*AGROW> 
    
    angle = angle + tiny_golden_angles(angle_nr);
end



%% export matrix

kspacelist = kspacelist(1:dimy*dimz,:);
disp(length(kspacelist))

if order == 1 
    ord = 'r';
else
    ord = 'o';
end

filename = strcat(outputdir,filesep,'Radial_y=',num2str(dimy),'_z=',num2str(dimz),'_a=',num2str(round(tiny_golden_angles(angle_nr),2)),'_',ord,'.txt');
fileID = fopen(filename,'w');

for i = 1:length(kspacelist)
   
    fprintf(fileID,num2str(kspacelist(i,1)));
    fprintf(fileID,'\n');
    fprintf(fileID,num2str(kspacelist(i,2)));
    fprintf(fileID,'\n');
    
end

fclose(fileID);



%% Display the trajectory

if display
    
    figure(1); %#ok<*UNRCH> 
    plot1 = scatter(kspacelist(1,1),kspacelist(1,2),'s');
    
    for i=1:length(kspacelist)
        plot1.XData = kspacelist(1:i,1);
        plot1.YData = kspacelist(1:i,2);
        pause(0.000001);
    end
    
end
