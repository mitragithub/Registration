function W = detect_padding_in_nissl(I)

% because several pixels are 255,255,255
% close all


% I = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/719/MD719&718-N37-2019.05.04-03.34.59_MD719_1_0109.tif');
% in this example there are lots of places with padding, but actually only
% the left stripe is padding
% I = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/719/MD719&718-N147-2019.05.04-18.26.10_MD719_2_0440.tif');
% in this example, the padding is not a whole row

% saveprefix = 'mask_ex1_'

if strcmp(class(I),'uint8')
    I = double(I)/255.0;
end
% figure;
% imagesc(I)
% axis image
% title Image
% saveas(gcf,[saveprefix 'image.png'])

W = all(I==1,3);
% figure
% imagesc(W)
% axis image
% title Mask


% the idea is I only want to call you padding if you are a "big rectangle"
% or a "big line"
% let's look for regions that have 10 in a row
big = 100;
O = ones(size(W));
Orb = convn(O,ones(big,1),'same');
Ocb = convn(O,ones(1,big),'same');
Wrb = convn(W,ones(big,1),'same') == Orb;
Wcb = convn(W,ones(1,big),'same') == Ocb;


small = 10;
Os = convn(O,ones(small,small),'same');
Ws = convn(W,ones(small,small),'same') == Os;




% figure;
% imagesc(Wrb)
% axis image
% title Wr


Wnew = ((W==Wrb) & W) | ((W==Wcb) & W ) | ((Ws==W) & W);
% Wnew = ((W==Wrb) & W) | ((W==Wcb) & W ) ;

% figure;
% imagesc(Wnew)
% axis image
% 
% figure;
% imagesc(cat(3,W,Wnew,W*0))
% axis image
% saveas(gcf,[saveprefix 'masks.png'])
% this is pretty good
% the only thing its missing are inside corners
% is there a way to get

W = Wnew;
