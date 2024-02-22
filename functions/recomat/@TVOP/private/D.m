function res = D(image)

res=zeros(size(image,1),size(image,2),size(image,3),3);

% image = a 2D + time image
%
% This function computes the finite difference transform in time of the image

Dt = image(:,:,[2:end,1]) - image;          

res(:,:,:,1) = 0;  
res(:,:,:,2) = 0;
res(:,:,:,3) = Dt;  % finite difference in time dimension only

end