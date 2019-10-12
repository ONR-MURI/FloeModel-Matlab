function [X,Y,c,cc, U,V, UU, VV] = create_eulerian_data(Floe, X,Y, c_fact, c2_boundary)


n_cutX=fix(size(X,1)/c_fact)*c_fact;
n_cutY=fix(size(X,2)/c_fact)*c_fact;
X=X(1:n_cutX,1:n_cutY);
Y=Y(1:n_cutX,1:n_cutY);
xx=X(1,:); yy=Y(:,1);
nx=length(xx); ny=length(yy);
inBoundaries=inpolygon(X(:),Y(:),c2_boundary(1,:),c2_boundary(2,:));

c=zeros(size(X));
U=zeros(size(X));
V=zeros(size(X));

xc=cat(1,Floe.Xi);
yc=cat(1,Floe.Yi);
alive=cat(1,Floe.alive);
rm=cat(1,Floe.rmax);

parfor j=1:length(Floe)
    if alive(j)

    x=Floe(j).Xg;   y=Floe(j).Yg;
    
    A=imrotate(Floe(j).A,-Floe(j).alpha_i/pi*180,'bilinear','crop');
    
    ind_x=logical( (xx> xc(j)-rm(j)) .* (xx< xc(j)+rm(j)));
    ind_y=logical( (yy> yc(j)-rm(j)) .* (yy< yc(j)+rm(j)));
            
    [xx_g, yy_g]=meshgrid(xx(ind_x), yy(ind_y));
    tmp=interp2(x+xc(j),y+yc(j),double(A),xx_g,yy_g); tmp(isnan(tmp))=0;
    Ag=zeros(ny,nx); Ag(ind_y,ind_x)=tmp; 
        
    c=c+Ag;    
    
    
    [Xg, Yg]=meshgrid(x,y); % grid centered around the ice floe
    
    [theta,rho] = cart2pol(Xg,Yg);
        
    Uice=Floe(j).Ui-rho*Floe(j).ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
    Vice=Floe(j).Vi+rho*Floe(j).ksi_ice.*cos(theta); % Y-dir velocity
       
    Uice(~A)=0; Vice(~A)=0;%to make zero where there is no floe
  
    tmp=interp2(x+xc(j),y+yc(j),double(Uice),xx_g,yy_g); tmp(isnan(tmp))=0;  

    ug=zeros(ny,nx); ug(ind_y,ind_x)=tmp;
    U=U+ug;
    
    tmp=interp2(x+xc(j),y+yc(j),double(Vice),xx_g,yy_g); tmp(isnan(tmp))=0;  
    vg=zeros(ny,nx); vg(ind_y,ind_x)=tmp;
    V=V+vg;
    

    end
    
end

c(~inBoundaries)=0;
V(~inBoundaries)=0; 
U(~inBoundaries)=0;
    
c(c>1)=1; % for concentration; otherwise it would be mass
 
cc=reshape(c,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
cc=squeeze(mean(mean(cc,1),3));

UU=reshape(U,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
UU=squeeze(mean(mean(UU,1),3));

VV=reshape(V,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
VV=squeeze(mean(mean(VV,1),3));
end

