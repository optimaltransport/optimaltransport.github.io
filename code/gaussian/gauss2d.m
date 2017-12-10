mu = [-2 0];
Sigma = .5*[1 -.5;-.5 1];

mu2 = [3 1];
Sigma2 = 2*[1 .25; .25 .5];

A= Sigma^(-1/2)*(Sigma^(1/2)*Sigma2*Sigma^(1/2))^(1/2)*Sigma^(-1/2)



x1 = -12:.05:12; x2 = -12:.05:12;

[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));

F2 = mvnpdf([X1(:) X2(:)],mu2,Sigma2);
F2 = reshape(F2,length(x2),length(x1));

x=mvnrnd(mu,1.7*Sigma,5)'; % sample with more variance
U=bsxfun(@plus,A*bsxfun(@minus,x,mu') ,mu2');
XX1=x(1,:); XX2=x(2,:)

figure('color','white');
%colormap([cool(64);hot(64)])
colormap(flipud(hot))
[C,h]=contourf(x1,x2,F,[.001 .01 .05:.05:.95 .99],'linewidth',1);
%clabel(C,h,'FontSize',15,'Color','red')

hold on
[C2,h2]=contourf(x1,x2,F2,[ .01 .02 .03 .05:.05:.95 .99],'linewidth',1);
%clabel(C,h2,'FontSize',15,'Color','red')
quiver(XX1,XX2,U(1,:)-XX1,U(2,:)-XX2,0,'color','black','linewidth',2,'MaxHeadSize',.1);
xlim([-5,7])
ylim([-3,4])
box off
grid off
%axis off


figure('color','white');
h=surf(X1,X2,100*reshape(F+.01,size(X1)),'FaceColor','interp','EdgeColor','none');
hold on
h2=surf(X1,X2,100*F2+.01,'FaceColor','interp','EdgeColor','none');
%quiver(XX1,XX2,U(1,:)-XX1,U(2,:)-XX2,0,'color','black','linewidth',2);
hh= @(x,Y) .5*(bsxfun(@minus,[x,Y],mu)*A*(bsxfun(@minus,[x;Y],mu'))) + mu2*[x;Y]+50;
h3=fsurf(hh,'ShowContours','on');
set(h3,'FaceAlpha',0.7,'EdgeColor','none');
xlim([-5,8])
ylim([-4,4])
zlim([0,150])
caxis([0,150])
view([26.5,20.4])
colormap([flipud(hot(64));cool(round(2.6*64))])
camlight(26,24)
%colormap(flipud(hot))
%box off
%grid off
%axis off
h3=gca;
set(gca,'ZTick',[]);

%axis tight


