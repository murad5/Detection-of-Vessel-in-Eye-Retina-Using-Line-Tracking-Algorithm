image=imread('01_test.tif');
%Split into RGB Channels
Red = image(:,:,1);
Green = image(:,:,2);
Blue = image(:,:,3);
%Get histValues for each channel
[yRed, x] = imhist(Red);
[yGreen, x] = imhist(Green);
[yBlue, x] = imhist(Blue);
%Plot them together in one plot
%figure;
%plot(x, yRed, 'Red', x, yGreen, 'Green', x, yBlue, 'Blue');

%I = rgb2gray(image);
I=Green;
%figure;
%imhist(I);
%figure;
%imshow(I);

% preprosesing 
edgeRetina = edge(I,'sobel',0.15);
% figure, imshow(edgeRetina);
seD = strel('disk',5);
dilateEdge = imdilate(edgeRetina, seD);
%figure,imshow(dilateEdge);title('Original dialate image');


[row,column,z]=size(image);
k = uint8(ones(row,column));



    for i = 1: row
        for j = 1:column;
            
               k(i,j) = image(i,j,2);  
              % disp(I(i,j));
           
        end
    end
figure;imshow(k);title('k');


%%%%% Local Normalization with green channel image
mymean = mean(double(k));
mystd = std(double(k));
mymean2=transpose(mymean);

normImage = uint8(ones(row,column));
   for i = 1: row
        for j = 1:column;
            
           normImage(i,j) = (k(i,j)*uint8(mymean(j)))+uint8(mystd(j));
           
        end
    end

lnfim=localnormalize(k,4,4);
lnfim2=mat2gray(lnfim);
figure,imshow(lnfim); title('Local Normalization in uint8');
%figure,imshow(lnfim2); title('Local Normalization in double');
lnSub=k-lnfim;
figure,imshow(lnSub); title('Substract from real image');

%de noise with median filter
   denoise = medfilt2(lnSub);
   
   %Adjust the contrast (Make it darker)
   figure, imshow(denoise); title('MedFiltered Grayscale Image');
   %contrast adjustment for better visibility
   darkImg = imadjust (denoise);
  
   %figure; imshow(darkImg); title('contrast adjusted')

   %we need to compliment the image as the vessels are dark now
   compliment = imcomplement(darkImg);
   %figure; imshow(compliment); title('complement')
  %histogram equilization
  hist = adapthisteq(compliment); 
  %figure; imshow(hist); title ('adjust');
  
  
  %structuring element ball shaped as the eye is ball shaped
  se = strel('ball',9,9);
  % opening
  fuOpen = imopen(hist,se);
  %figure; imshow(fuOpen); title ('opening')
  %new structuring element
  se2 = strel('line',1 ,4);
  fuuOpen = imerode(fuOpen, se2);
  %figure, imshow(fuuOpen); title('disk')
  
  
  % Remove Optic Disk around eye
  
  opDisk = hist - fuuOpen; 
  %figure; imshow(opDisk); title ('removing edge disc');
  
  %2D Median Filter
  %again de noising
  medFilt = medfilt2(opDisk); 
  %figure; imshow(medFilt); title ('medFilt')


  
  %adjusting contrast
  I3 = imadjust(medFilt); 
  
  BW = edge(I3);
  %figure,imshow(BW);title('Edges of Original Image');
    
  %line tracking masks
  
  kernelHorz=[-1,-1,-1;
      2,2,2;
      -1,-1,-1];
  
  kernelVert=[-1,2,-1;-1,2,-1; -1 ,2, -1];
  
  kernelDiag=[-1, -1 , 2;  -1, 2 , -1; 2 , -1 , -1];
  
  kernelDiag2=[2,-1,-1; -1, 2, -1; -1, -1 ,2 ];
  
  [r,c]=size(I3);
R=zeros(r,c,'uint8');

 R1 = imfilter(I3 ,kernelHorz);
R2 = imfilter(I3 ,kernelVert);
R3 = imfilter(I3 ,kernelDiag);
R4 = imfilter(I3 ,kernelDiag2);



  for i = 1:r
    for j = 1:c
      
    R(i,j)= max([abs(R1(i,j)), abs(R2(i,j)), abs(R3(i,j)), abs(R4(i,j))]);

 
    end
  end
  
 
 %figure,imshow(R); title('Masking');
 
 
 se3 = strel('line',10 ,4);
 fuuOpen2 = imdilate(R, se3);
 %figure, imshow(fuuOpen2); title('Dialate mask')
 
 maskSub = fuuOpen2-R;
 %figure, imshow(maskSub); title('Substract mask')
 
 %post processing into binary image for better visual of extrracted
  %vessels
  
  % Gray Thresholding
  level2 = graythresh(fuuOpen2);
  level = graythresh(R);
  
   % Binarization
  bw = imbinarize(R,level);  
  bw2 = imbinarize(fuuOpen2,level2);  
  %   this basically removes the small objects  in Image Containing Fewer Than 50 Pixels
  
  %figure, imshow(bw); title('graythresh 1 mask')
  %figure, imshow(bw2); title('graythresh 2 mask')
  
  bw = bwareaopen(bw,50);   
  
 subU=bw;
 subEdge = bw-dilateEdge;
 figure,imshow(subEdge); title('edge subtract');
 
 seDline = strel('line',1,3);
 

 BW2 = imdilate(subEdge,seDline); 
 
 
 figure,imshow(BW2); title('vessels'); 
  
 localMean = imboxfilt(BW2,7);
 figure,imshow(localMean); title('vessels using box filter');
 newk = double(ones(row,column));
 for i = 1: row
    for j = 1:column;
            
        if localMean(i,j)>0.2
           if (localMean(i,j)+0.2)>1.0
              newk(i,j)=1.0; 
           else
              newk(i,j)= localMean(i,j)+0.2;
           end
        else
           newk(i,j)=0;  
        end
           
    end
end
figure,imshow(newk); title('vessels after threshholding'); 
 %histogram equilization
  hist1 = adapthisteq(localMean); 
  figure,imshow(hist1); title('vessels histogram equilization'); 
  b = imsharpen(hist1,'Radius',3,'Amount',15);
figure, imshow(b)
title('Sharpened Image');
  
 %{
 newk = uint8(ones(row,column));
 for i = 1: row
    for j = 1:column;
            
        if localMean(i,j)>0.4
           newk(i,j)=255; 
        else
           newk(i,j)=0;  
        end
           
    end
end
figure,imshow(newk); title('vessels after binarizing'); 

b = imsharpen(localMean,'Radius',2,'Amount',1);
figure, imshow(b)
title('Sharpened Image');

%contrast adjustment for better visibility
   darkImg = imadjust (localMean);
   figure,imshow(darkImg); title('vessels darker'); 
%histogram equilization
  hist1 = adapthisteq(darkImg); 
  figure,imshow(hist1); title('vessels histogram equilization'); 
 %figure,imshow(bw); title('vessels');
 
 %Laplacian
 filter=[-1 -1 -1;-1 9 -1; -1 -1 -1];
result=hist1;
for i=2:r-1
    for j=2:c-1
        sum=0;
        row=0;
        col=1;
        
        for k=i-1:i+1
            row=row+1;
            col=1;
            for l=j-1:j+1
                sum = sum+hist1(k,l)*filter(row,col);               
                col=col+1;
            end
        end
      result(i,j)=sum;      
    end
end
lapAdd = result+hist1;
selap = strel('line',2 ,4);
lapErode = imerode(lapAdd, selap);
figure,imshow(lapErode); title('vessels after laplacian errosion'); 
 
erSub = lapErode-hist1;
figure,imshow(lapErode); title('vessels after substract errosion'); 
 
figure;
imhist(k);
%}
