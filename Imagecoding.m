clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Created By Malinda Sulochana Silva               %
           % manawaya@hotmail.com                             %
           % Faculty of Electrical and Electronic Engineering % 
           % University of Peradeniya                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

help Huffman;
% Reading the Image data
    I_original = imread('Parrots-680x680.jpg');
    [raws,cols,numofcolorChanels] = size(I_original);
    %figure;
    %imshow(I_original);
% Crop the Image (Starting Point = 208)
    %___________________________________________________
    % change these parameters to Custamize your cropping
   
    width = 16;
    length = 16;
    startpos = 208;
    %___________________________________________________
    
    I_cropped = Huffman.h_crop(length, width ,startpos, I_original);

% Converting into gray image
    I_cropped = rgb2gray(I_original);
    %figure;
    %imshow(I_cropped);
% Quantizing the Cropped image
    %______________________________________________________________
    % Change this parameter to change the no.of Quantization levels
    quantization_levels =7;
    %______________________________________________________________
    % defining 7 threshould levels
    thresh = multithresh(I_cropped,quantization_levels);
    I_quantized =imquantize(I_cropped,thresh); % Quantizing into_
                                               % threshold lvls
                                               
% Finding the probability of each symbol distribution,
    [probability, symbol] = hist(I_quantized,unique(I_quantized));

% Finding each symbol probability
	symbol_probability = (sum(probability,2)/sum(sum(probability,2)));
   
% Generate huffman dictionary
    dict = Huffman.huffman_dict(symbol_probability, symbol);   
%__________________________________________________________________________
% Encode image with generated Huffman codes (cropped image)
    I_encoded_cropped = Huffman.h_encode(I_quantized, width, length, ...
                            quantization_levels, dict);
% Decode image with generated Huffman codes (cropped image)       
    I_decoded_cropped = Huffman.h_decode(I_encoded_cropped, width, length, ...
                                            dict,thresh,quantization_levels);
    %figure;
    %imshow(I_decoded_cropped);
%__________________________________________________________________________    
% Encode image with generated Huffman codes (Original image)
    %Quantizing the Original image
    I_quantized1 =imquantize(I_original,thresh); % Quantizing into_
                                               % threshold lvls
  
    I_encoded_cmplete = Huffman.h_encode(I_quantized1, raws, cols, ...
                            quantization_levels, dict);
                        
% Decode image with generated Huffman codes (Encoded image)       
    I_decoded_cmplete = Huffman.h_decode(I_encoded_cmplete, raws, ...
                                 cols, dict,thresh,quantization_levels);
    %figure;
    %imshow(I_decoded_cmplete);
%__________________________________________________________________________    
% Finding the entropy of the source image
    entropy_original = Huffman.h_entropy(I_original);
    fprintf('Entropy of the Original image\t\t :%f\n', entropy_original);
% Finding the entropy of the decoded image
    entropy_decoded = Huffman.h_entropy(I_decoded_cmplete);
    fprintf('Entropy of the decoded image\t\t :%f\n', entropy_decoded);
% Finding the entropy of the decoded image
    entropy_cropped = Huffman.h_entropy(I_decoded_cropped);
    fprintf('Entropy of the decoded_cropped image\t\t :%f\n', ...
                                                  entropy_cropped);
%__________________________________________________________________________    
% Evaluating the PSNR
    PSNR_original = Huffman.h_PSNR(rgb2gray(I_original), ...
                                                 rgb2gray(I_original));
    fprintf('PSNR of the Original imag\t\t :%f\n', PSNR_original);
    
    PSNR_decoded = Huffman.h_PSNR(rgb2gray(I_original),I_decoded_cmplete);
    fprintf('PSNR of the decoded image\t\t :%f\n', PSNR_decoded);
%__________________________________________________________________________    
% Evaluating the compression ratio
    Compression_ratio = Huffman.h_compratio(I_original,I_encoded_cmplete);
    fprintf('Compression Ratio\t\t :%f\n', Compression_ratio);
%__________________________________________________________________________
% Evaluating the average length of the encoded_cropped image
    avglen = 0;
    for i = 1:width
        for j = 1:length
             avglen = avglen + size(I_encoded_cropped{i,j},2);
        end
    end
    fprintf('average Length of encoded_cropped image\t\t :%f\n', avglen);
    
    