function [ f ] = DataDensityPlot( x, y, levels )

map = dataDensity(x, y, 256, 256);
%levels = round(max(max(map)));
%map = map - min(min(map));
map = floor(map ./ max(max(map)) * (levels-1));

image(map);
colormap(jet(levels));
set(gca, 'XTick', [1 256]);
set(gca, 'XTickLabel', [min(x) max(x)]);
set(gca, 'YTick', [1 256]);
set(gca, 'YTickLabel', [min(y) max(y)]);


    function [ dmap ] = dataDensity( x, y, width, height, limits, fudge )

        if(nargin == 4)
            limits(1) = min(x);
            limits(2) = max(x);
            limits(3) = min(y);
            limits(4) = max(y);
        end
        deltax = (limits(2) - limits(1)) / width;
        deltay = (limits(4) - limits(3)) / height;
        if(nargin < 6)
            fudge = sqrt(deltax^2 + deltay^2);
        end
        dmap = zeros(height, width);
        for ii = 0: height - 1
            yi = limits(3) + ii * deltay + deltay/2;
            for jj = 0 : width - 1
                xi = limits(1) + jj * deltax + deltax/2;
                dd = 0;
                for kk = 1: length(x)
                    dist2 = (x(kk) - xi)^2 + (y(kk) - yi)^2;
                    dd = dd + 1 / ( dist2 + fudge);
                end
                dmap(ii+1,jj+1) = dd;
            end
        end

    end



end