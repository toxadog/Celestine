%%
cd 'D:\SP_step1\Mouse 9\Experiment 2022-04-05T 17.37.42 awake spont'
W_S_P0 = load('WavestartPoint.mat');
XY0 = W_S_P0.WavestartPoint;
cd 'D:\SP_step1\Mouse 9\Experiment 2022-04-05T 18.13.34 awake spont\Converted brain'
BlueImage0 = double(fliplr(rot90(read(Tiff('Image_000201.tif')))));
NewImage = BlueImage0;

cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.34.55 1.5% spont'
W_S_P1 = load('WavestartPoint.mat');
XY1 = W_S_P1.WavestartPoint;
cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.34.55 1.5% spont\Converted brain'
BlueImage1 = double(fliplr(rot90(read(Tiff('Image_000201.tif')))));
Reference = BlueImage1;
cd 'D:\SP_step1\Mouse 9\trace_waves'
[~,tform] = registerImages(Reference,NewImage);
XY1_c = zeros(size(XY1));
for t = 1:size(BlueImage1,2)
    XYt = [XY1(1,t), XY1(2,t), 1]*(tform.T.*[1 1 1;1 1 1;1 -1 1]);
    XY1_c(:,t) = [XYt(1);XYt(2)];
end

cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.13.30 2%  spont'
W_S_P2 = load('WavestartPoint.mat');
XY2 = W_S_P2.WavestartPoint;
cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.13.30 2%  spont\Converted brain'
BlueImage2 = double(fliplr(rot90(read(Tiff('Image_000201.tif')))));
Reference = BlueImage2;
cd 'D:\SP_step1\Mouse 9\trace_waves'
[~,tform] = registerImages(Reference,NewImage);
XY2_c = zeros(size(XY2));
for t = 1:size(BlueImage2,2)
    XYt = [XY2(1,t), XY2(2,t), 1]*(tform.T.*[1 1 1;1 1 1;1 -1 1]);
    XY2_c(:,t) = [XYt(1);XYt(2)];
end

cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 17.12.00 3% spont'
W_S_P3 = load('WavestartPoint.mat');
XY3 = W_S_P3.WavestartPoint;
cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 17.12.00 3% spont\Converted brain'
BlueImage3 = double(fliplr(rot90(read(Tiff('Image_000201.tif')))));
Reference = BlueImage3;
cd 'D:\SP_step1\Mouse 9\trace_waves'
[~,tform] = registerImages(Reference,NewImage);
XY3_c = zeros(size(XY3));
for t = 1:size(BlueImage3,2)
    XYt = [XY3(1,t), XY3(2,t), 1]*(tform.T.*[1 1 1;1 1 1;1 -1 1]);
    XY3_c(:,t) = [XYt(1);XYt(2)];
end

%% making histogram manually 

sub_area = 8;
h0 = zeros(304/sub_area);
for Xcnt = 0:sub_area:304-sub_area
   for Ycnt = 0:sub_area:304-sub_area
    
    h0(Xcnt/sub_area+1,Ycnt/sub_area+1) = sum((XY0(1,:)<Xcnt+sub_area & XY0(1,:)> Xcnt).*(XY0(2,:)<Ycnt+sub_area & XY0(2,:)> Ycnt));
       
   end
end 

h1 = zeros(304/sub_area);
for Xcnt = 0:sub_area:304-sub_area
   for Ycnt = 0:sub_area:304-sub_area
    
    h1(Xcnt/sub_area+1,Ycnt/sub_area+1) = sum((XY1_c(1,:)<Xcnt+sub_area & XY1_c(1,:)> Xcnt).*(XY1_c(2,:)<Ycnt+sub_area & XY1_c(2,:)> Ycnt));
       
   end
end 

h2 = zeros(304/sub_area);
for Xcnt = 0:sub_area:304-sub_area
   for Ycnt = 0:sub_area:304-sub_area
    
    h2(Xcnt/sub_area+1,Ycnt/sub_area+1) = sum((XY2_c(1,:)<Xcnt+sub_area & XY2_c(1,:)> Xcnt).*(XY2_c(2,:)<Ycnt+sub_area & XY2_c(2,:)> Ycnt));
       
   end
end 

h3 = zeros(304/sub_area);
for Xcnt = 0:sub_area:304-sub_area
   for Ycnt = 0:sub_area:304-sub_area
    
    h3(Xcnt/sub_area+1,Ycnt/sub_area+1) = sum((XY3_c(1,:)<Xcnt+sub_area & XY3_c(1,:)> Xcnt).*(XY3_c(2,:)<Ycnt+sub_area & XY3_c(2,:)> Ycnt));
       
   end
end 


% ----
figure();
subplot(1,2,1)
imagesc(rot90(h0,+1))
pbaspect([1 1 1])
colorbar

subplot(1,2,2)
scatter(XY0(1,:),XY0(2,:))
pbaspect([1 1 1])
xlim([1 304]);ylim([1 304])
sgtitle('start - Awake','FontSize',20, ...
'FontName','Times New Roman');

figure();
subplot(1,2,1)
imagesc(rot90(h1,+1))
pbaspect([1 1 1])
colorbar

subplot(1,2,2)
scatter(XY1_c(1,:),XY1_c(2,:))
pbaspect([1 1 1])
xlim([1 304]);ylim([1 304])
sgtitle('start - 1.5%','FontSize',20, ...
'FontName','Times New Roman');

figure();
subplot(1,2,1)
imagesc(rot90(h2,+1))
pbaspect([1 1 1])
colorbar

subplot(1,2,2)
scatter(XY2_c(1,:),XY2_c(2,:))
pbaspect([1 1 1])
xlim([1 304]);ylim([1 304])
sgtitle('start - 2%','FontSize',20, ...
'FontName','Times New Roman');

figure();
subplot(1,2,1)
imagesc(rot90(h3,+1))
pbaspect([1 1 1])
colorbar

subplot(1,2,2)
scatter(XY3_c(1,:),XY3_c(2,:))
pbaspect([1 1 1])
xlim([1 304]);ylim([1 304])
sgtitle('start - 3%','FontSize',20, ...
'FontName','Times New Roman');

%% making histogram using hist function
% [X,Y] = find(ones(304,304));
% nbin = 16;
% figure();
% subplot(1,4,1)
% hist3([X(:)',XY0(1,:);Y(:)',XY0(2,:)]','Nbins',[nbin nbin],'CDataMode','auto','FaceColor','interp')
% % [~, c0] = hist3([XY0(1,:);XY0(2,:)]','Nbins',[nbin nbin]);
% % H0 = hist3([XY0(2,:);XY0(1,:)]','CdataMode','auto','Nbins',[nbin nbin]);
% % % H0(H0>10) = 0;
% % surf(c0{1},c0{2},H0);
% caxis([361 515])
% title('AW')
% colormap(jet)
% view(2)
% pbaspect([1 1 1])
% colorbar
% xlim([1 304]);ylim([1 304])
% 
% subplot(1,4,2)
% [r,c] = find(XY1_c==[0;0]);
% XY1_c(:,c)=[];
% hist3([X(:)',XY1_c(1,:);Y(:)',XY1_c(2,:)]','Nbins',[nbin nbin],'CDataMode','auto','FaceColor','interp')
% % [~, c1] = hist3([XY1_c(2,:);XY1_c(1,:)]','Nbins',[nbin nbin]);
% % H1 = hist3([XY1_c(2,:);XY1_c(1,:)]','CdataMode','auto','Nbins',[nbin nbin]);
% % % H1(H1>10) = 0;
% % surf(c1{1},c1{2},H1);
% caxis([361 380])
% title('1.5%')
% colormap(jet)
% view(2)
% pbaspect([1 1 1])
% colorbar
% xlim([1 304]);ylim([1 304])
% 
% subplot(1,4,3)
% [r,c] = find(XY2_c==[0;0]);
% XY2_c(:,c)=[];
% hist3([X(:)',XY2_c(1,:);Y(:)',XY2_c(2,:)]','Nbins',[nbin nbin],'CDataMode','auto','FaceColor','interp')
% % [~, c2] = hist3([XY2_c(2,:);XY2_c(1,:)]','Nbins',[nbin nbin]);
% % H2 = hist3([XY2_c(2,:);XY2_c(1,:)]','CdataMode','auto','Nbins',[nbin nbin]);
% % % H2(H2>30) = 0;
% % surf(c2{1},c2{2},H2);
% caxis([361 380])
% title('2%')
% colormap(jet)
% view(2)
% pbaspect([1 1 1])
% colorbar
% xlim([1 304]);ylim([1 304])
% 
% subplot(1,4,4)
% [r,c] = find(XY3_c==[0;0]);
% XY3_c(:,c)=[];
% hist3([X(:)',XY3_c(1,:);Y(:)',XY3_c(2,:)]','Nbins',[nbin nbin],'CDataMode','auto','FaceColor','interp')
% % [~, c3] = hist3([XY3_c(2,:);XY3_c(1,:)]','Nbins',[nbin nbin]);
% % H3 = hist3([XY3_c(2,:);XY3_c(1,:)]','CdataMode','auto','Nbins',[nbin nbin]);
% % % H3(H3>10) = 0;
% % surf(c3{1},c3{2},H3);
% caxis([361 380])
% title('3%')
% colormap(jet)
% view(2)
% pbaspect([1 1 1])
% colorbar
% xlim([1 304]);ylim([1 304])