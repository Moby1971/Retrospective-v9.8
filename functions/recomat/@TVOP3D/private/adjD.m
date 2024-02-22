function res = adjD(y)

res = zeros(size(y,1),size(y,2),size(y,3),size(y,4),4);

% only gradient in time dimension
res = adjDt(y(:,:,:,:,4));

end

function res = adjDz(x)

res = x([1,1:end-1],:,:,:) - x;
res(1,:,:,:) = x(end,:,:,:);
res(end,:,:,:) = x(end-1,:,:,:);

end

function res = adjDy(x)

res = x(:,[1,1:end-1],:,:) - x;
res(:,1,:,:) = x(:,end,:,:);
res(:,end,:,:) = x(:,end-1,:,:);

end

function res = adjDx(x)

res = x(:,:,[1,1:end-1],:) - x;
res(:,:,1,:) = x(:,:,end,:);
res(:,:,end,:) = x(:,:,end-1,:);

end

function res = adjDt(x)

res = x(:,:,:,[1,1:end-1]) - x;
res(:,:,:,1) = x(:,:,:,end);
res(:,:,:,end) = x(:,:,:,end-1);

end