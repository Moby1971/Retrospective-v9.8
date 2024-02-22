function res = D(image)

res=zeros(size(image,1),size(image,2),size(image,3),size(image,4),4);

% image = 3D + time image
%
% This function computes the finite difference transform of the image
%


Dz = image([2:end,1],:,:,:) - image;
Dy = image(:,[2:end,1],:,:) - image;
Dx = image(:,:,[2:end,1],:) - image;
Dt = image(:,:,:,[2:end,1]) - image;          

res(:,:,:,:,1) = 0;  
res(:,:,:,:,2) = 0;
res(:,:,:,:,3) = 0;
res(:,:,:,:,4) = Dt;  % finite difference in time dimension only

end