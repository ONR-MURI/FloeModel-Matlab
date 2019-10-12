function [P, worked] = sort_contact_points(P, c1, c2 )

%P_initial=P;
dist=zeros(1,length(P));
for i=1:length(P)
    dist(i)=betweenpoints(P(:,i), c1 );
end
[~, b]=sort(dist); P=P(:,b); % points are sorted by distance along the contour


N_contact=size(P,2)/2;

moved=0;
worked=0;


floe2_in_floe1=zeros(1,N_contact);
while ~worked && moved<2
    
    for k=1:N_contact
        
        p=P(1:2,2*k-1:2*k);
        
        pcenter(k,:)=mean(p,2); dp=p(:,1)-p(:,2);
        
        d=1e5; %checking intersections with countour at +-d distance from the center point.
        
        dirr=dp(1)/dp(2);
        if abs(dirr)>1e3;
            
            x1=pcenter(k,1);
            y1=pcenter(k,2) - d;
            
            x2=pcenter(k,1);
            y2=pcenter(k,2) + d;
            
        else
            x1=pcenter(k,1) + d/sqrt(1+dirr^2);
            y1=pcenter(k,2) - d/sqrt(1+dirr^2)*dirr;
            
            x2=pcenter(k,1) - d/sqrt(1+dirr^2);
            y2=pcenter(k,2) + d/sqrt(1+dirr^2)*dirr;
        end
        
        
        pp2=InterX(c2,[x1 x2; y1 y2]);
        [~, i_min] = sort(sqrt( ( (pp2(1,:)-pcenter(k,1)).^2 + (pp2(2,:)-pcenter(k,2)).^2)));
        pp2=pp2(:,i_min(1));
        
        
        floe2_in_floe1(k)=inpolygon(pp2(1),pp2(2),c1(1,:),c1(2,:));
    end
    %display(floe2_in_floe1);
    
    if sum(floe2_in_floe1) == N_contact,
        worked=1;
        
    elseif  moved==0
        P = [P(:,end) P(:,1:end-1)]; moved=1; % if still didn't work, shuffle points
        
    elseif moved==1
        %display('Contact points could not be sorted!');              
        %P=P_initial; %returning the initial guess.
        moved=2;
    end
    
end

% 
% if ~worked
%     figure;
%     plot(c1(1,:),c1(2,:),'r'); hold on;
%     plot(c2(1,:),c2(2,:),'b');
%     for k=1:N_contact, plot(P(1,2*k-1:2*k),P(2,2*k-1:2*k),'x'); end
% end



end

