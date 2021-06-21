function [vx,vy] = vertexcrop(x,y,vx,vy)
%%This function removes any points outside the domain and replaces them
%%with points intersection the boundary for floe initialization

nnx=0.25*abs(max(x)-min(x));
nny=0.25*abs(max(y)-min(y));
for ii = length(vx):-1:1
    if vx(1,ii) < min(x)-nnx && vx(2,ii) < min(x)-nnx
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif isnan(vx(1,ii)) == 1 || isnan(vx(2,ii)) == 1
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif vx(1,ii) < min(x)-nnx
        ynew = interp1(vx(:,ii),vy(:,ii),min(x)-nnx);
        vx(1,ii) = min(x)-nnx;
        vy(1,ii) = ynew;
    elseif vx(2,ii)<min(x)-nnx
        ynew = interp1(vx(:,ii),vy(:,ii),min(x)-nnx);
        vx(2,ii) = min(x)-nnx;
        vy(2,ii) = ynew;
    end
end
for ii = length(vx):-1:1
    if vx(1,ii) > max(x)+nnx && vx(2,ii) > max(x)+nnx
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif vx(1,ii) > max(x)+nnx
        ynew = interp1(vx(:,ii),vy(:,ii),max(x)+nnx);
        vx(1,ii) = max(x)+nnx;
        vy(1,ii) = ynew;
    elseif vx(2,ii) > max(x)+nnx
        ynew = interp1(vx(:,ii),vy(:,ii),max(x)+nnx);
        vx(2,ii) = max(x)+nnx;
        vy(2,ii) = ynew;
    end
end
for ii = length(vx):-1:1
    if vy(1,ii) < min(y)-nny && vy(2,ii) < min(y)-nny
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif isnan(vy(1,ii)) == 1 || isnan(vy(2,ii)) == 1
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif vy(1,ii) < min(y)-nny
        xnew = interp1(vy(:,ii),vx(:,ii),min(y)-nny);
        vy(1,ii) = min(y)-nny;
        vx(1,ii) = xnew;
    elseif vy(2,ii) < min(y)-nny
        xnew = interp1(vy(:,ii),vx(:,ii),min(y)-nny);
        vy(2,ii) = min(y)-nny;
        vx(2,ii) = xnew;
    end
end
for ii = length(vx):-1:1
    if vy(1,ii) > max(y)+nny && vy(2,ii) > max(y)+nny
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif isnan(vy(1,ii)) == 1 || isnan(vy(2,ii)) == 1
        vx(:,ii) = NaN;
        vy(:,ii) = NaN;
    elseif vy(1,ii) > max(y)+nny
        xnew = interp1(vy(:,ii),vx(:,ii),max(y)+nny);
        vy(1,ii) = max(y)+nny;
        vx(1,ii) = xnew;
    elseif vy(2,ii) > max(y)+nny
        xnew = interp1(vy(:,ii),vx(:,ii),max(y)+nny);
        vy(2,ii) = max(y)+nny;
        vx(2,ii) = xnew;
    end
end