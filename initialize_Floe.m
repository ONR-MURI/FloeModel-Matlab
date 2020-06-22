function Floe=initialize_Floe(FloeShapes)


disp('Initializing Floes...');
thickness=2; % ice floe thickness in meters; save for all floes.
rho_ice=920;

floePixelSize=120; % in meters.

% floe shapes and locations are stored .mat file; run create_floe_pics.m to generate floes from images.
load('FloeShapes.mat','FloeShapes');

centr_mean=mean(cat(1,FloeShapes.Centroid));
ind=1;
total_multiple_contour_floes=0;
for j=1:length(FloeShapes)
    dist=max(abs(FloeShapes(j).Centroid-centr_mean))*floePixelSize;
    if dist < 65e3
        floe_position(ind,:)=(FloeShapes(j).Centroid-centr_mean)*floePixelSize;
        mask(ind).floe = FloeShapes(j).Image;
        
        A=mask(ind).floe;
        
        
        x=1:size(A,2); y=1:size(A,1);
        [x,y]=meshgrid(floePixelSize*x,floePixelSize*y); %% floePixelSize m is the pixel size from Earthdata NASA webpage
        xc=mean(x(A)); yc=mean(y(A));
        
        r=sqrt((x-xc).^2+(y-yc).^2); r_max=max(r(A))+2*floePixelSize; % pad with extra row/column of zeros around the floe
        
        mask(ind).rmax=r_max;
        mask(ind).area=sum(sum(A))*floePixelSize^2;
        mask(ind).mass=mask(ind).area*thickness*rho_ice; % total mass
        mask(ind).inertia_moment=sum(r(A).^2)*floePixelSize^2*thickness*rho_ice; % moment of inertia
        
        Ngr=2*fix(r_max/floePixelSize); % number of grid boxes over floe
        
        [X,Y]=meshgrid((-1:2/Ngr:1)*r_max, (-1:2/Ngr:1)*r_max); % floe grid
        mask(ind).X=X; mask(ind).Y=Y;
        mask(ind).Xg=(-1:2/Ngr:1)*r_max; mask(ind).Yg=mask(ind).Xg;
        
        A=griddata(x-xc,y-yc,double(A),X,Y); A(isnan(A))=0;
        A=imbinarize(A,'global');
        if sum(size(A))>0
            mask(ind).floe=A;
            
            ksi_ini=0; Ui_ini=0; Vi_ini=0; alpha_ini=0;
            %ksi_ini=1e-4*(2*rand(1)-1); Ui_ini=(2*rand(1)-1)*0.1; Vi_ini=(2*rand(1)-1)*0.1;alpha_ini=0;
            Xi=floe_position(ind,1);
            Yi=floe_position(ind,2);
            
            c0=contourc(mask(ind).Xg,mask(ind).Yg,double(A),[0.5 0.5]);
            
            if ~isempty(c0)
                if c0(2,1)< size(c0,2)-1,
                    multiple_contours=1;
                    %display('multiple contours');
                    total_multiple_contour_floes=total_multiple_contour_floes+1;
                else multiple_contours=0;
                end
                
                c0(:,c0(1,:)==0.5)=NaN;
                c0=c0(:,2:end);
                A_rot=[cos(alpha_ini) -sin(alpha_ini); sin(alpha_ini) cos(alpha_ini)]; %rotation matrix
                c_alpha=A_rot*c0;
                
                L_boundary=68e3;
                out_of_boundary=max(max(mask(ind).X(mask(ind).floe)))+Xi>L_boundary || min(min(mask(ind).X(mask(ind).floe)))+Xi<-L_boundary || max(max(mask(ind).Y(mask(ind).floe)))+Yi>L_boundary || min(min(mask(ind).Y(mask(ind).floe)))+Yi<-L_boundary;
                if (multiple_contours ||  out_of_boundary)
                    %disp('Ice floe is out of ocean grid bounds');
                else
                    
                    Floe(ind).area=mask(ind).area;
                    Floe(ind).mass=mask(ind).mass;
                    Floe(ind).inertia_moment=mask(ind).inertia_moment;
                    
                    Floe(ind).c0=c0; %contour
                    Floe(ind).c_alpha=c_alpha;
                    Floe(ind).A=uint8(A);
                    %           Floe(ind).A_alpha=imrotate(A,-alpha_ini/pi*180,'bilinear','crop');
                    Floe(ind).rmax=mask(ind).rmax;
                    Floe(ind).Xg=mask(ind).Xg;
                    Floe(ind).Yg=mask(ind).Yg;
                    Floe(ind).X=mask(ind).X;
                    Floe(ind).Y=mask(ind).Y;
                    
                    Floe(ind).Xi=Xi;
                    Floe(ind).Yi=Yi;
                    Floe(ind).alpha_i=alpha_ini;
                    
                    Floe(ind).Ui=Ui_ini;
                    Floe(ind).Vi=Ui_ini;
                    Floe(ind).ksi_ice=ksi_ini;
                    Floe(ind).alive=1;
                    
                    Floe(ind).dXi_p=0;
                    Floe(ind).dYi_p=0;
                    Floe(ind).dUi_p=0;
                    Floe(ind).dVi_p=0;
                    Floe(ind).dalpha_i_p=0;
                    Floe(ind).dksi_ice_p=0;
                    
                    ind=ind+1;
                end
            end
        end
    end
end

display(['Good floes: ' num2str(ind) ]);
display(['Out of boundary floes: ' num2str(length(FloeShapes)-ind-total_multiple_contour_floes)]);
display(['Multiple contour floes:' num2str(total_multiple_contour_floes)]);


end

