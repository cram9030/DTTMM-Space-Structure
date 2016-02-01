function capMov(X,L,filename)

writerObj = VideoWriter(filename);
writerObj.FrameRate = 30;
open(writerObj);

%axis([0 1000 min([min(min(X(1,:,:))) min(min(Y(:,1:2:40)))]) max([max(max(X(1,:,:))) max(max(Y(:,1:2:40)))])])
axis([0 sum(L) min([min(min(X(1,:,:)))]) max([max(max(X(1,:,:)))])])
%axis([0 1000 min([min(min(Y(:,1:2:40)))]) max([max(max(Y(:,1:2:40)))])])
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

[~,~,n] = size(X);

for k = 1:n 
   plot([0 cumsum(L)'], X(1,:,k))%,[0 cumsum(L)'],[0 Y(k,1:2:40)])
   frame = getframe;
   writeVideo(writerObj,frame);
end

close(writerObj);