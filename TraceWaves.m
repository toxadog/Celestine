tiffdir_Cortex = fullfile('Converted brain');
tiff_files_Cortex = dir(fullfile(tiffdir_Cortex,'*.tif'));
tiff_files_Cortex = {tiff_files_Cortex.name};
%%
cnt=1;
clear('A');
for i=2:2:800
    t2B = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{i}));
    imageCortexBlue = read(t2B);
    A(:,:,cnt)=imageCortexBlue;
    cnt=cnt+1;
end
[n1,n2]=size(imageCortexBlue);
%%
STDVideo=std(double(A)/max(max(max(double(A)))),[],3);
%%
STDVideoF=medfilt2(STDVideo,[101 101]);
MapMax=max(max(STDVideoF));
MapMin=min(min(STDVideoF));
MapRange=MapMax-MapMin;
%%
thr=0.5;
MapThr=MapMin+thr*MapRange;
%%
mask=STDVideoF>MapThr;
%%
figure;
subplot(1,2,1)
imagesc(mean(A,3))
subplot(1,2,2)
imagesc(mean(A,3).*mask)
%%
save('Mask.mat','mask')
%%
connectivity=6;
thrH=0.1e-3;
thrW=0.1e-2;
%% Create F0
lengthF0=100;
BlueArray=zeros(304,304,lengthF0);
t2B_0 = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{2}));
imageCortexBlue0 = read(t2B_0);
BlueArray(:,:,1) = double(imageCortexBlue0);

for i=2:lengthF0
    t2B = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{i*2}));  
    imageCortexBlue = read(t2B);
%     imageCortexBlue=imregister(imageCortexBlue,imageCortexBlue0,'rigid',optimizer,metric);
    BlueArray(:,:,i) = double(imageCortexBlue);
end

GreenArray=zeros(304,304,lengthF0);
t2G_0 = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{1}));
imageCortexGreen0 = read(t2G_0);
GreenArray(:,:,1) = double(imageCortexGreen0);

for i=2:lengthF0
    t2G = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{i*2-1}));  
    imageCortexGreen = read(t2G);
%     imageCortexGreen=imregister(imageCortexGreen,imageCortexGreen0,'rigid',optimizer,metric);
    GreenArray(:,:,i) = double(imageCortexGreen);
end

imageCortexF0=(mean(BlueArray,3))./(mean(GreenArray,3));
%%
sc1=20;
sc2=30;
Rmin=20;
PeaksArray=[];
PeaksArrayPrevious=[];
PeaksArrayInd=[];
PeaksArrayIndPrevious=[];
Waves=[];
wavecounter=0;
figure;
set(gcf,'Position',[3 260 1803 718]);
v = VideoWriter('Waves.avi');
v.FrameRate=100;
open(v);
for i=1:floor(length(tiff_files_Cortex)/2)-10
    t2B = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{i*2}));
    imageCortexBlue = read(t2B);
    t2G = Tiff(horzcat(tiffdir_Cortex,'\',tiff_files_Cortex{i*2-1}));
    imageCortexGreen = read(t2G);
    
    

    if i>round(lengthF0/2)
        if mod(i+round(lengthF0/2),lengthF0)~=0
            BlueArray(:,:,mod(i+round(lengthF0/2),lengthF0)) = imageCortexBlue;
            GreenArray(:,:,mod(i+round(lengthF0/2),lengthF0)) = imageCortexGreen;
            imageCortexF0=(mean(BlueArray,3))./(mean(GreenArray,3));
          
            
            
        else
            BlueArray(:,:,lengthF0) = imageCortexBlue;
            GreenArray(:,:,lengthF0) = imageCortexGreen;
            imageCortexF0=(mean(BlueArray,3))./(mean(GreenArray,3));
            
            
            
        end
    end
    imageCortex=double(imageCortexBlue)./double(imageCortexGreen);
    A=(imageCortex-imageCortexF0)./imageCortexF0;
    
    Im1=repmat((double(imageCortexGreen)/400*1.2),[1,1,3]);
    Im2=zeros(304,304,3);
    Im2(:,:,2)=mask.*A/0.1*0.5;
    A = A.*mask + (~mask)*mean(A(logical(mask)));
    A=waveletDecompose(A,20,30);
    WHmax=imextendedmax(A,thrH,connectivity);
    WHmaxCon=bwconncomp(WHmax);
    % Measure blob properties
    Peaks = regionprops(WHmaxCon,A,'WeightedCentroid','Centroid','PixelList','PixelValues');
    PeaksArray=zeros(size(Peaks,1),2);
    cnt=1;
    subplot(1,2,1)
    imshow(fliplr(rot90(Im1+Im2)));
    xticklabels([])
    yticklabels([])
    ylim([1 304])
    xlim([1 304])
    hold on
        
    for j=1:size(Peaks,1)
        [~,nmax]=max(Peaks(j).PixelValues);
        if (mask(Peaks(j).PixelList(nmax,2),Peaks(j).PixelList(nmax,1)))&&(A(Peaks(j).PixelList(nmax,2),Peaks(j).PixelList(nmax,1))>thrW)
            PeaksArray(cnt,:)=Peaks(j).PixelList(nmax,:);
            plot(305-Peaks(j).PixelList(nmax,2),305-Peaks(j).PixelList(nmax,1),'o','markerfacecolor',[1,1,1],'markeredgecolor',[0,0,0])

            cnt=cnt+1;
        end
    end
    hold off


    subplot(1,2,2)
    imagesc(fliplr(rot90(A)))
    caxis([-0.01 0.01]);
    axis equal
    xticklabels([])
    yticklabels([])
    ylim([1 304])
    xlim([1 304])
    hold on
        
    for j=1:size(Peaks,1)
        [~,nmax]=max(Peaks(j).PixelValues);
        if (mask(Peaks(j).PixelList(nmax,2),Peaks(j).PixelList(nmax,1)))&&(A(Peaks(j).PixelList(nmax,2),Peaks(j).PixelList(nmax,1))>thrW)
            PeaksArray(cnt,:)=Peaks(j).PixelList(nmax,:);
            plot(305-Peaks(j).PixelList(nmax,2),305-Peaks(j).PixelList(nmax,1),'o','markerfacecolor',[0,0,0],'markeredgecolor',[0,0,0]);
        end
    end
    hold off
    PeaksArray=PeaksArray(1:cnt-1,:);


   
    
    
    if isempty(PeaksArrayPrevious)
      for j=1:size(PeaksArray,1)
        wavecounter=wavecounter+1;
        Waves.(horzcat('p',int2str(wavecounter))).x=PeaksArray(j,1);
        Waves.(horzcat('p',int2str(wavecounter))).y=PeaksArray(j,2);
        Waves.(horzcat('p',int2str(wavecounter))).t0=i;
        Waves.(horzcat('p',int2str(wavecounter))).tn=i;
        Waves.(horzcat('p',int2str(wavecounter))).Amp=A(PeaksArray(j,2),PeaksArray(j,1));
        PeaksArrayInd=[PeaksArrayInd wavecounter];
      end
    else
        DistMatrix=createDistMatrix(PeaksArrayPrevious,PeaksArray);
        MinDistPos = checkMatr(DistMatrix,Rmin);
        for j=1:length(MinDistPos)
            if ~isnan(MinDistPos(j))
                Waves.(horzcat('p',int2str(PeaksArrayIndPrevious(MinDistPos(j))))).x=...
                    [Waves.(horzcat('p',int2str(PeaksArrayIndPrevious(MinDistPos(j))))).x;PeaksArray(j,1)];
                Waves.(horzcat('p',int2str(PeaksArrayIndPrevious(MinDistPos(j))))).y=...
                    [Waves.(horzcat('p',int2str(PeaksArrayIndPrevious(MinDistPos(j))))).y; size(A,2)+1-PeaksArray(j,2)];
                Waves.(horzcat('p',int2str(PeaksArrayIndPrevious(MinDistPos(j))))).tn=i;
                Waves.(horzcat('p',int2str(wavecounter))).Amp=...
                    [Waves.(horzcat('p',int2str(wavecounter))).Amp;A(PeaksArray(j,2),PeaksArray(j,1))];
                PeaksArrayInd=[PeaksArrayInd PeaksArrayIndPrevious(MinDistPos(j))];
            else
                wavecounter=wavecounter+1;
                Waves.(horzcat('p',int2str(wavecounter))).x=PeaksArray(j,1);
                Waves.(horzcat('p',int2str(wavecounter))).y=size(A,2)+1-PeaksArray(j,2);
                Waves.(horzcat('p',int2str(wavecounter))).t0=i;
                Waves.(horzcat('p',int2str(wavecounter))).tn=i;
                Waves.(horzcat('p',int2str(wavecounter))).Amp=A(PeaksArray(j,2),PeaksArray(j,1));          
                PeaksArrayInd=[PeaksArrayInd wavecounter];
            end
            
        end

    end
    PeaksArrayIndPrevious=PeaksArrayInd;
    PeaksArrayPrevious=PeaksArray;
    PeaksArrayInd=[];
    PeaksArray=[];
    
    
    
    for j=1:size(Peaks,1)
        [~,nmax]=max(Peaks(j).PixelValues);
        Peaks(j).MaxCoord=Peaks(j).PixelList(nmax,:);
        if mask(Peaks(j).PixelList(nmax,2),Peaks(j).PixelList(nmax,1))
        end
    end


  frame = getframe(gcf);
  writeVideo(v,frame);   
  if ~mod(i,100)
        display(i);
  end

end
close(v)
%%
WaveNames=fieldnames(Waves);
figure;
hold on;
for i=1:numel(WaveNames)
    plot(305-Waves.(WaveNames{i}).y,305-Waves.(WaveNames{i}).x);
end
%%
WaveNames=fieldnames(Waves);
NumColors=numel(WaveNames);
figure;imshow(fliplr(rot90(Im1*1.0)));
hold on;
for i=1:numel(WaveNames)
    if length(Waves.(WaveNames{i}).x)>2
        %         plot((Waves.(WaveNames{i}).x),305-(Waves.(WaveNames{i}).y),'color',[1-1e-4*i 1-3e-4*i 1-3e-4*i]);
        tempx=medfilt1((Waves.(WaveNames{i}).x),2);
        tempy=medfilt1((Waves.(WaveNames{i}).y),2);
        tempx=tempx(2:end);
        tempy=tempy(2:end);
        plot(305-tempy,305-tempx,'color',[1-0.3*i/NumColors 1-1*i/NumColors 1-1*i/NumColors]);
    end
end
%%
WaveNames=fieldnames(Waves);
MaxIntens=0;
MinIntens=1e10;
for i=1:numel(WaveNames)
    if abs(mean(Waves.(WaveNames{i}).Amp))<MinIntens
       MinIntens = abs(mean(Waves.(WaveNames{i}).Amp));
    end
    if abs(mean(Waves.(WaveNames{i}).Amp))>MaxIntens
       MaxIntens = abs(mean(Waves.(WaveNames{i}).Amp));
    end    
end

%%
figure;imshow(flipud(Im1));
hold on;
for i=1:numel(WaveNames)
    ColorInd = abs(mean(Waves.(WaveNames{i}).Amp))/(MaxIntens-MinIntens)-MinIntens/(MaxIntens-MinIntens);
    if length(Waves.(WaveNames{i}).x)>2
        %         plot((Waves.(WaveNames{i}).x),305-(Waves.(WaveNames{i}).y),'color',[1-1e-4*i 1-3e-4*i 1-3e-4*i]);
        tempx=medfilt1((Waves.(WaveNames{i}).x),2);
        tempy=medfilt1((Waves.(WaveNames{i}).y),2);
        tempx=tempx(2:end);
        tempy=tempy(2:end);
        plot(tempx,305-tempy,'color',[1-0.3*ColorInd 1-1*ColorInd 1-1*ColorInd]);
    end
end
%%
Amp = zeros(1,wavecounter);
Times = zeros(1,wavecounter);
Length = zeros(1,wavecounter);
for i=1:wavecounter
Amp(i)=mean(Waves.(horzcat('p',int2str(i))).Amp);
Times(i)=length(Waves.(horzcat('p',int2str(i))).x);
Length(i)=sqrt((Waves.(horzcat('p',int2str(i))).x(1)-Waves.(horzcat('p',int2str(i))).x(end)).^2+(Waves.(horzcat('p',int2str(i))).y(1)-Waves.(horzcat('p',int2str(i))).y(end)).^2);
end
%%
Thr=1e-4;
q=0.001;
WaveNames=fieldnames(Waves);
NumColors=numel(WaveNames);
figure;imshow(fliplr(rot90(Im1*1)));
hold on;
for i=1:NumColors
    if (Amp(i)>Thr)&&(L(i)>1)
        %         plot((Waves.(WaveNames{i}).x),305-(Waves.(WaveNames{i}).y),'color',[1-1e-4*i 1-3e-4*i 1-3e-4*i]);
        CurrentWave=Waves.(WaveNames{i});
        CurrentWaveFilt=WaveKalman(CurrentWave,q,R);
        tempx=(CurrentWaveFilt.x);
        tempy=(CurrentWaveFilt.y);
        
        plot(305-tempy,305-tempx,'color',[1-1*i/NumColors 1-0.3*i/NumColors 1-0.1*i/NumColors]);
    end
end

