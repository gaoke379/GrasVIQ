% % %*********************************************************************************
% * Filename: MarkFG.m
% *
% * Created by
% * Dr. Filiz Bunyak
% * Department of Electrical Engineering and Computer Science
% * 219 Naka Hall
% * University of Missouri-Columbia
% * Columbia, MO 65211-2060
% * bunyak@missouri.edu
% *
% % %*********************************************************************************

function [img_overlay]= MarkFG(img,mask,color,varargin)

if (isempty(varargin))
   mix_ratio=[0 1];
else
   mix_ratio=double(varargin{1});
   mix_ratio=mix_ratio/sum(mix_ratio);
end

color=double(color);
[~,~,channels]=size(img);
max_im=max(max(max(img)));
if (max_im>2)
    img=double(img)/255;
end
if (channels==1)
    R=img;
    G=img;
    B=img;
else
    R=img(:,:,1);
    G=img(:,:,2);
    B=img(:,:,3);
end
R(mask==1)=mix_ratio(1)*R(mask==1)+ mix_ratio(2)*color(1);
G(mask==1)=mix_ratio(1)*G(mask==1)+ mix_ratio(2)*color(2);
B(mask==1)=mix_ratio(1)*B(mask==1)+ mix_ratio(2)*color(3);

img_overlay(:,:,1)=uint8(R*255);
img_overlay(:,:,2)=uint8(G*255);
img_overlay(:,:,3)=uint8(B*255);

      
