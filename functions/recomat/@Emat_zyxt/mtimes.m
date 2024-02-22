function res = mtimes(a,b)

if a.adjoint
    
    % x=H'*y (x=res,b=y), x: object, y: multi-coil k-space data
    
    % multi-coil data in image domain
    for ch = 1:size(b,5)
        
        x_array(:,:,:,:,ch) = ifft3c_mri(b(:,:,:,:,ch).*a.mask);
        
    end
    
    % multi-coil combination in the image domain
    for tt = 1:size(b,4)
        
        res(:,:,:,tt) = sum(squeeze(x_array(:,:,:,tt,:)).*conj(a.b1),4)./sum(abs((a.b1)).^2,4); 
        
    end
    
else
    
    % y=H*x (y=res,b=x), x: object, y: multi-coil k-space data
    
    % multi-coil image from object
    for tt=1:size(b,4)
        
        for ch=1:size(a.b1,4)
            
            x_array(:,:,:,tt,ch) = b(:,:,:,tt).*a.b1(:,:,:,ch); 
            
        end
        
    end
    
    % multi-coil image to k-space domain
    res = fft3c_mri(x_array);
    
    % apply sampling mask
    for ch=1:size(a.b1,4)
        
        res(:,:,:,:,ch) = res(:,:,:,:,ch).*a.mask;
        
    end
    
end