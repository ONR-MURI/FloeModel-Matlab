function figg = plotAtmData(figg,atmdata,Time,ind)

clf(figg);

colormap default
%clf(figg),
imagesc([-1e6,1e6], [-1e6,1e6], atmdata,'CDataMapping','scaled');
colorbar

if(ind==1)
    title(['psihigh: Time = ' num2str(Time) ' s']);
elseif(ind==2)
    title(['Time = ' num2str(Time) ' s']);
elseif(ind==3)
    title(['qhigh: Time = ' num2str(Time) ' s']);
else
    title(['qlow: Time = ' num2str(Time) ' s']);
end

xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

end
