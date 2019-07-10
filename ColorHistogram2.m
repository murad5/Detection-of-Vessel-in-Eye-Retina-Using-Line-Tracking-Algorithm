clear all;
image=imread('01_test.tif');
%Split into RGB Channels
Red = image(:,:,1);
Green = image(:,:,2);
Blue = image(:,:,3);

GreenChannelImage=Green;

% preprosesing 
edgeRetina = edge(GreenChannelImage,'sobel',0.15);
% figure, imshow(edgeRetina);
seD = strel('disk',7);
dilateEdge = imdilate(edgeRetina, seD);
figure,imshow(dilateEdge);title('Original dialate image');


[row,column,z]=size(image);
k=GreenChannelImage;
figure;imshow(k);title('k');


%%%%% Local Normalization with green channel image
localNormalized=localnormalize(GreenChannelImage,4);

smoothenedImage=k-localNormalized;
figure,imshow(localNormalized); title('Substract from real image');

%de noise with median filter
denoise1 = medfilt2(smoothenedImage);
figure, imshow(denoise1); title('MedFiltered Grayscale Image');

%Adjust the contrast (Make it darker)
%contrast adjustment for better visibility
darkImg = imadjust (denoise1);
%figure; imshow(darkImg); title('contrast adjusted')

%we need to compliment the image as the vessels are dark now
complement = imcomplement(darkImg);
%figure; imshow(compliment); title('complement')

%histogram equilization
histEqu = adapthisteq(complement); 
figure; imshow(histEqu); title ('adjust');
  
  
%structuring element ball shaped as the eye is ball shaped
se = strel('ball',9,9);
% opening
fuOpen = imopen(histEqu,se);
%figure; imshow(fuOpen); title ('opening')
%new structuring element
se2 = strel('line',1 ,4);
fuuOpen = imerode(fuOpen, se2);
%figure, imshow(fuuOpen); title('disk')
  
  
% Remove Optic Disk around eye
opDisk = histEqu - fuuOpen; 
figure; imshow(opDisk); title ('removing edge disc');
  
%2D Median Filter
%again de noising
denoise2 = medfilt2(opDisk); 
%figure; imshow(medFilt); title ('medFilt')

%adjusting contrast
I3 = imadjust(denoise2); 
BW = edge(I3);
%figure,imshow(BW);title('Edges of Original Image');
    
%line tracking masks
kernelHorz=[-1,-1,-1;
             2, 2, 2;
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
  
 
figure,imshow(R); title('Masking');

 
%post processing into binary image for better visual of extrracted
%vessels
  
% Gray Thresholding

level = graythresh(I3);
  
% Binarization
bw = imbinarize(I3,level);  

%   this basically removes the small objects  in Image Containing Fewer Than 50 Pixels
  
  figure, imshow(bw); title('graythresh 1 mask')
  %figure, imshow(bw2); title('graythresh 2 mask')
  
bw2 = bwareaopen(bw,50);   
subU=bw;
subEdge = bw2-dilateEdge; %Substract edge
figure,imshow(subU); title('edge subtract');
 
%seDline = strel('line',1,3);
%BW2 = imdilate(subEdge,seDline); 
%figure,imshow(BW2); title('vessels'); 
  
localMean = imboxfilt(subEdge,7);
figure,imshow(localMean); title('vessels using box filter');
 
%histogram equilization
hist1 = adapthisteq(localMean); 
figure,imshow(hist1); title('vessels histogram equilization'); 

b = imsharpen(hist1,'Radius',3,'Amount',15);
figure, imshow(b)
title('Sharpened Image');


 