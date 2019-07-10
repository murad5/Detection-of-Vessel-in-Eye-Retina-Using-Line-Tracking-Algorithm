function varargout = BloodVesselGui(varargin)
% BLOODVESSELGUI MATLAB code for BloodVesselGui.fig
%      BLOODVESSELGUI, by itself, creates a new BLOODVESSELGUI or raises the existing
%      singleton*.
%
%      H = BLOODVESSELGUI returns the handle to a new BLOODVESSELGUI or the handle to
%      the existing singleton*.
%
%      BLOODVESSELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLOODVESSELGUI.M with the given input arguments.
%
%      BLOODVESSELGUI('Property','Value',...) creates a new BLOODVESSELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BloodVesselGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BloodVesselGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BloodVesselGui

% Last Modified by GUIDE v2.5 10-Apr-2019 23:10:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BloodVesselGui_OpeningFcn, ...
                   'gui_OutputFcn',  @BloodVesselGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BloodVesselGui is made visible.
function BloodVesselGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BloodVesselGui (see VARARGIN)

%Datasets


% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BloodVesselGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BloodVesselGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in preButton.
function preButton_Callback(hObject, eventdata, handles)
% hObject    handle to preButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GreenChannelImage;handles.dilateEdge;handles.localNormalized;handles.smoothenedImage;
%cla(handles.axes4);

axes(handles.axes4);
imshow(handles.GreenChannelImage);
title('Green Channel');
axes(handles.axes5);
imshow(handles.dilateEdge);
title('Edge Detection');
axes(handles.axes6);
imshow(handles.localNormalized);
title('Local Normalization');
axes(handles.axes7);
imshow(handles.smoothenedImage);
title('Smoothened Image');


axes(handles.axes8);
imshow(handles.black);
axes(handles.axes9);
imshow(handles.black);
axes(handles.axes10);
imshow(handles.black);


% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in browseButton.
function browseButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[basefilename,path]= uigetfile({'*.tif'},'Open Tif Image File');
filename= fullfile(path, basefilename);
handles.image = imread (filename);
% if I = [MxNx4]
if(size(handles.image,3)==4)
    handles.image(:,:,4)=[]; % convert to I = [MxNx3]
end

% if I = [MxN]
if(size(handles.image,3)==1)
    [handles.image]=gray2rgb(handles.image); % convert to I = [MxNx3]
%     figure;imshow(I);
end
[row,column,z]=size(handles.image);
handles.black=uint8(255*ones(row,column));

axes(handles.axes2);
imshow(handles.image);

% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in analyzeButton.
function analyzeButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Split into RGB Channels
handles.Red = handles.image(:,:,1);
handles.Green = handles.image(:,:,2);
handles.Blue = handles.image(:,:,3);

handles.GreenChannelImage=handles.Green;

% preprosesing 
handles.edgeRetina = edge(handles.GreenChannelImage,'sobel',0.15);
% figure, imshow(edgeRetina);
handles.seD = strel('disk',7);
handles.dilateEdge = imdilate(handles.edgeRetina, handles.seD);
%figure,imshow(dilateEdge);title('Original dialate image');


[row,column,z]=size(image);
handles.k = handles.GreenChannelImage;
%figure;imshow(k);title('k');


%%%%% Local Normalization with green channel image
handles.localNormalized=localnormalize(handles.GreenChannelImage,4);

handles.smoothenedImage=handles.k- handles.localNormalized;
%figure,imshow(smoothenedImage); title('Substract from real image');

%de noise with median filter
handles.denoise1 = medfilt2(handles.smoothenedImage);
%figure, imshow(denoise1); title('MedFiltered Grayscale Image');

%Adjust the contrast (Make it darker)
%contrast adjustment for better visibility
handles.darkImg = imadjust (handles.denoise1);
%figure; imshow(darkImg); title('contrast adjusted')

%we need to compliment the image as the vessels are dark now
handles.complement = imcomplement(handles.darkImg);
%figure; imshow(compliment); title('complement')

%histogram equilization
handles.histEqu = adapthisteq(handles.complement); 
%figure; imshow(histEqu); title ('adjust');
  
  
%structuring element ball shaped as the eye is ball shaped
handles.se = strel('ball',9,9);
% opening
handles.fuOpen = imopen(handles.histEqu,handles.se);
%figure; imshow(fuOpen); title ('opening')
%new structuring element
handles.se2 = strel('line',1 ,4);
handles.fuuOpen = imerode(handles.fuOpen, handles.se2);
%figure, imshow(fuuOpen); title('disk')
  
  
% Remove Optic Disk around eye
handles.opDisk = handles.histEqu - handles.fuuOpen; 
%figure; imshow(opDisk); title ('removing edge disc');
  
%2D Median Filter
%again de noising
handles.denoise2 = medfilt2(handles.opDisk); 
%figure; imshow(medFilt); title ('medFilt')

%adjusting contrast
handles.I3 = imadjust(handles.denoise2); 
handles.BW = edge(handles.I3);
%figure,imshow(BW);title('Edges of Original Image');
    
%line tracking masks
handles.kernelHorz=[-1,-1,-1;
             2, 2, 2;
            -1,-1,-1];
  
handles.kernelVert=[-1,2,-1;-1,2,-1; -1 ,2, -1];
  
handles.kernelDiag=[-1, -1 , 2;  -1, 2 , -1; 2 , -1 , -1];
  
handles.kernelDiag2=[2,-1,-1; -1, 2, -1; -1, -1 ,2 ];
  
[handles.r,handles.c]=size(handles.I3);
handles.R=zeros(handles.r,handles.c,'uint8');

handles.R1 = imfilter(handles.I3 ,handles.kernelHorz);
handles.R2 = imfilter(handles.I3 ,handles.kernelVert);
handles.R3 = imfilter(handles.I3 ,handles.kernelDiag);
handles.R4 = imfilter(handles.I3 ,handles.kernelDiag2);



for i = 1:handles.r
   for j = 1:handles.c 
        handles.R(i,j)= max([abs(handles.R1(i,j)), abs(handles.R2(i,j)), abs(handles.R3(i,j)), abs(handles.R4(i,j))]);
   end
end
  
 
%figure,imshow(R); title('Masking');

 
%post processing into binary image for better visual of extrracted
%vessels
  
% Gray Thresholding

handles.level = graythresh(handles.R);
  
% Binarization
handles.bw = imbinarize(handles.R,handles.level);  

%   this basically removes the small objects  in Image Containing Fewer Than 50 Pixels
  
  %figure, imshow(bw); title('graythresh 1 mask')
  %figure, imshow(bw2); title('graythresh 2 mask')
  
handles.bw2 = bwareaopen(handles.bw,50);   
handles.subU= handles.bw;
handles.subEdge = handles.bw2 - handles.dilateEdge; %Substract edge
%figure,imshow(subU); title('edge subtract');
 
%seDline = strel('line',1,3);
%BW2 = imdilate(subEdge,seDline); 
%figure,imshow(BW2); title('vessels'); 
  
handles.localMean = imboxfilt(handles.subEdge,7);
%figure,imshow(localMean); title('vessels using box filter');
 
%histogram equilization
handles.hist1 = adapthisteq(handles.localMean); 
%figure,imshow(hist1); title('vessels histogram equilization'); 

handles.b = imsharpen(handles.hist1,'Radius',3,'Amount',15);
%figure, imshow(b)
%title('Sharpened Image');

axes(handles.axes3);
imshow(handles.b);


axes(handles.axes4);
imshow(handles.black);
axes(handles.axes5);
imshow(handles.black);
axes(handles.axes6);
imshow(handles.black);
axes(handles.axes7);
imshow(handles.black);

axes(handles.axes8);
imshow(handles.black);
axes(handles.axes9);
imshow(handles.black);
axes(handles.axes10);
imshow(handles.black);

% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in adjustButton.
function adjustButton_Callback(hObject, eventdata, handles)
% hObject    handle to adjustButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GreenChannelImage;handles.dilateEdge;handles.localNormalized;handles.smoothenedImage;
%cla(handles.axes4);

axes(handles.axes4);
imshow(handles.denoise1);
title('Noise Removal');
axes(handles.axes5);
imshow(handles.darkImg);
title('Darker Image');
axes(handles.axes6);
imshow(handles.complement);
title('Inverse Image');
axes(handles.axes7);
imshow(handles.histEqu);
title('Equalized Image');

axes(handles.axes8);
imshow(handles.black);
axes(handles.axes9);
imshow(handles.black);
axes(handles.axes10);
imshow(handles.black);


% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in remButton.
function remButton_Callback(hObject, eventdata, handles)
% hObject    handle to remButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GreenChannelImage;handles.dilateEdge;handles.localNormalized;handles.smoothenedImage;
%cla(handles.axes4);

axes(handles.axes4);
imshow(handles.fuOpen);
title('Ball Filtering');
axes(handles.axes5);
imshow(handles.fuuOpen);
title('Line Filtering');
axes(handles.axes6);
imshow(handles.opDisk);
title('Removing Edge Disk');

axes(handles.axes7);
imshow(handles.black);
axes(handles.axes8);
imshow(handles.black);
axes(handles.axes9);
imshow(handles.black);
axes(handles.axes10);
imshow(handles.black);

% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in lineButton.
function lineButton_Callback(hObject, eventdata, handles)
% hObject    handle to lineButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes4);
imshow(handles.denoise2);
title('Noise Removal');
axes(handles.axes5);
imshow(handles.I3);
title('Adjustment');
axes(handles.axes6);
imshow(handles.R);
title('Line Tracking');

axes(handles.axes7);
imshow(handles.black);
axes(handles.axes8);
imshow(handles.black);
axes(handles.axes9);
imshow(handles.black);
axes(handles.axes10);
imshow(handles.black);


% Choose default command line output for BloodVesselGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in postButton.
function postButton_Callback(hObject, eventdata, handles)
% hObject    handle to postButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes4);
imshow(handles.level);
title('Grey Thresholding');

axes(handles.axes5);
imshow(handles.bw);
title('Binarization');

axes(handles.axes6);
imshow(handles.bw2);
title('');

axes(handles.axes7);
imshow(handles.subEdge);
title('Edge Removal');

axes(handles.axes8);
imshow(handles.localMean);
title('Double Edge Removal');

axes(handles.axes9);
imshow(handles.hist1);
title('Equalization');

axes(handles.axes10);
imshow(handles.b);
title('Sharpened Image');
