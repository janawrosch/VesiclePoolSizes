function [imgv]=thresh_img(img,n_sig,max_man, show_control_figures)
img=double(img);
[hy,hx]=hist(img(:),500);
shy=smooth(hy);
if show_control_figures==1
figure;plot(hx,shy);
title('smooth histogram')
end
if max_man==0
[mx,my]=max(shy(1:end-round(0.2*length(shy))));
else
[none,my]=min(abs(hx-max_man));
mx=shy(my);
end
gy=[shy(1:my);flipud(shy(1:my));zeros((size(hy,2)-(2*my)),1)];
gx=hx';

% gx und gy gleich lang machen
n=min(length(gy),length(gx));
gx=gx(1:n);
gy=gy(1:n);

fit_gauss = fittype('a1*exp(-((x-b1)/c1)^2)');
a1_start = mx;
b1_start = hx(my);
c1_start = my;%/2;
fit_options = fitoptions('Method' , 'NonlinearLeastSquares' ,'Robust', 'off' , 'Algorithm', 'Trust-Region', 'DiffMinChange',1e-08,'DiffMaxChange',0.1,'TolFun',1e-6,'TolX',1e-6,'Startpoint', [a1_start b1_start c1_start], 'Lower' ,[0 0 0],'Upper',[inf inf inf]);
[c_gauss] = fit(gx,gy,fit_gauss, fit_options );
a1 = c_gauss.a1;
b1 = c_gauss.b1;
c1 = c_gauss.c1;
imgv.bg_sigma = c_gauss.c1;
for i=1:size(gx,1)
yf(i)=a1*exp(-((gx(i)-b1)/c1)^2);
end
%diff=shy-yf;
if show_control_figures==1
figure; plot(hx,shy); hold on; plot(gx,yf,'g');
title('Histgram with fit')
end
imgv.bg_val=hx(my);
clear x
clear y
imgv.img_bgs=img-hx(my);
imgv.img_bgs=(double(imgv.img_bgs>0)).*imgv.img_bgs;
% figure;imagesc(img);colormap(gray);
if show_control_figures==1
figure;imagesc(imgv.img_bgs);colormap(gray);
title('Diff Image with background removed')
end
%imgv.sg_val=mean(mean(imgv.img_bgs(333:343,243:248)));
imgv.img_over_bg=(img>(imgv.bg_val+(n_sig*imgv.bg_sigma))).*img;
imgv.threshold=imgv.bg_val+(n_sig*imgv.bg_sigma);
if show_control_figures==1
figure; imagesc(imgv.img_over_bg);colormap(gray);
title('Values above threshold')
end

