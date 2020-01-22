function [X,Y,c,cc, U,V, UU, VV] = create_eulerian_data(Floe, X,Y, c_fact)

n_cutX=fix(size(X,1)/c_fact)*c_fact;
n_cutY=fix(size(X,2)/c_fact)*c_fact;
X=X(1:n_cutX,1:n_cutY);
Y=Y(1:n_cutX,1:n_cutY);

c=zeros(size(X));
U=zeros(size(X));
V=zeros(size(X));

parfor j=1:length(Floe)
    if Floe(j).alive
    xc=Floe(j).Xi; yc=Floe(j).Yi;
    x=Floe(j).Xg;   y=Floe(j).Yg;
    
    A=imrotate(Floe(j).A,-Floe(j).alpha_i/pi*180,'crop');

    Ag=interp2(x+xc,y+yc,double(A),X,Y); Ag(isnan(Ag))=0;   
    c=c+Ag;    
    
    [Xg, Yg]=meshgrid(x,y); % grid centered around the ice floe
    
    [theta,rho] = cart2pol(Xg,Yg);
        
    Uice=Floe(j).Ui-rho*Floe(j).ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
    Vice=Floe(j).Vi+rho*Floe(j).ksi_ice.*cos(theta); % Y-dir velocity
       
    Uice(~A)=0; Vice(~A)=0;%to make zero where there is no floe
  
    ug=interp2(x+xc,y+yc,double(Uice),X,Y); ug(isnan(ug))=0;   
    U=U+ug;
    vg=interp2(x+xc,y+yc,double(Vice),X,Y); vg(isnan(vg))=0;   
    V=V+vg;
    end
    
end

c(c>1)=1; % for concentration; otherwise it would be mass
 
cc=reshape(c,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
cc=squeeze(mean(mean(cc,1),3));


UU=reshape(U,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
UU=squeeze(mean(mean(UU,1),3));

VV=reshape(V,[c_fact,n_cutX/c_fact,c_fact,n_cutY/c_fact]);
VV=squeeze(mean(mean(VV,1),3));

end

