function res = adjD(y)

res = zeros(size(y,1),size(y,2),size(y,3),3);

% Only TV in time dimension
res = adjDt(y(:,:,:,3));

end

function res = adjDt(x)

res = x(:,:,[1,1:end-1]) - x;
res(:,:,1) = x(:,:,end);
res(:,:,end) = x(:,:,end-1);

end
