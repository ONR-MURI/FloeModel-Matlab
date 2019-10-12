function [dist, sideNum]= betweenpoints(c, floeContour )

c1=floeContour(:,1:end-1);

c2=floeContour(:,2:end);

c=repmat(c,1,length(c1));

a=sum((c1-c).*(c2-c),1)./sqrt(sum((c1-c).*(c1-c),1))./sqrt(sum((c2-c).*(c2-c),1));

ind=1:length(c1); sideNum=ind(a<-1+1e-5);

if ~isempty(ind)
     dist=sqrt((c(:,sideNum)-c1(:,sideNum))'*(c(:,sideNum)-c1(:,sideNum)));
     if sideNum>1
         dist=dist+sum(sqrt(sum((c2(:,1:sideNum)-c1(:,1:sideNum)).*(c2(:,1:sideNum)-c1(:,1:sideNum)),1)));
     end
else sideNum=0; dist=0;
    
end


end

