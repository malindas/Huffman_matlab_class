%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Created By Malinda Sulochana Silva               %
           % manawaya@hotmail.com                             %
           % Faculty of Electrical and Electronic Engineering % 
           % University of Peradeniya                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%==========================================================================
%                              Huffman Class
%==========================================================================
%
%  === METHODS >>
%-- [obj] = huffman_dict(probability, symbol) will return a
%   structre with huffman dictionary values.
%   probability and symbols dimensions should be same;
% 
%-- [ret] = h_sort(huff_tree) will sort the huff_tree structure 
%   in assending order considering the probability values
%
%-- [I_cropped] = h_crop(length, width ,startpos, I_original) will 
%   crop the image
%
%-- [I_encoded] = h_encode(I_quantized, width, length, ...
%                                            qantization_levels,dict)
%    will encode the Image with generated Huffman dictionary
%
%-- [I_decoded] = h_decode(I_encoded, width, length, ...
%                                            dict,thresh,quantization_levels)
%   will decode the encoded image
%
%-- [entropy] = h_entropy(Image) gives the entropy of the image
%
%-- [psnr] = h_PSNR(Original, Decoded) gives the PSNR of the Image
%
%-- [comp_ratio] = h_compratio(original, encoded) will gives the
%   compression ratio of the image
%
%  === STRUCTURES >>
%      huff_tree = struct('symbol', [], 'probability', [],...
%                                 'child', [], 'code', []); 
%
%____________________Copyright.All Rights Reserved_________________________ 
%==========================================================================



classdef Huffman < handle
     
     methods(Static,Access = public)
    
% @Breif    :sort Huffman table in assending order
% @param    :Huffman Table structure
% @ret      :Huffman Table structure      
        function ret = h_sort(huff_tree)
          
            for i = 1:length(huff_tree) 
                for j = 1:length(huff_tree) 
                    if huff_tree(i).probability<huff_tree(j).probability
                       temp = huff_tree(i);
                       huff_tree(i) = huff_tree(j);
                       huff_tree(j) = temp;
                    end
                end
            end          
          ret = huff_tree;
        end
        

% @Breif    : Create a Huffman dictionary
% @param    : symbol, probabilitiessymbol
% @ret      : Huffman dict structrue 
        function dict= huffman_dict(probability, symbol)
              % huffman structrue to store node data
              huff_tree = struct('symbol', [], 'probability', [],...
                                  'child', [], 'code', []); 
              
              if length(probability) ~= length(probability)
                  disp('Error Matrix dimensions mismatch');
                  return;
              end
              % Create a node for each                 
              for i = 1:length(probability)
                  huff_tree(i).symbol = symbol(i);
                  huff_tree(i).probability = probability(i);  
              end
              huff_tree = Huffman.h_sort(huff_tree);

              i=1;
              while length(huff_tree) > 1
                 
                        huff_tree = Huffman.h_sort(huff_tree);    
                        newNode = huff_tree(1); % prototyping structure; 
                        newNode.symbol = ['Node', num2str(i)];
                        huff_tree(1).code = 1;
                        huff_tree(2).code = 0;
                        newNode.child =[huff_tree(1),huff_tree(2)]; 
                        newNode.probability = huff_tree(1).probability+ ...
                                              huff_tree(2).probability;           
                        
                        huff_tree = [newNode huff_tree(3:end)];
                        huff_tree = Huffman.h_sort(huff_tree);       
                   i = i + 1 ;
              end   
              % calling the  genCode function 
               dict = Huffman.h_genCode(huff_tree,{},[],1);%full_tree;
               
              % sorting the dict table in assending order
               for i = 1:length(dict) 
                    for j = 1:length(dict)
                       if dict{i,1} < dict{j,1}
                       temp = dict{i,1};
                       tmp = dict{i,2};
                       dict{i,1} = dict{j,1};
                       dict{i,2} = dict{j,2};
                       dict{j,1} = temp;
                       dict{j,2} = tmp;
                       end
                    end    
               end
        end
        
        
% This function will run through the huffman tree and generates the code
% book.
% This is a recersive function. Used as a private function 
        function [code_book, acc, index] = h_genCode(huff_tree, ...
                                                code_book, acc, index)
            
            %if reaches the end of a branch
            if isempty(huff_tree.child)
                 for i=1:length(huff_tree)                         
                     acc =[acc, huff_tree.code]; 
                     % generate the code book @ an end of a branch
                     code_book{index,1} = huff_tree.symbol;
                     code_book{index,2} = acc;    
                     index = index + 1;
                     % left last code when retruning to a node
                     acc = acc(1:end-1); 
                 end  
            else
                % accumilate code through node paths
                for j=1:length(huff_tree.child)           
                  acc = [acc, huff_tree.code];        
                  [code_book, acc, index] = Huffman.h_genCode( ...
                      huff_tree.child(j),code_book, acc, index);  
                end 
            end 
          % left last code when returning to previous node
          acc = acc(1:end-1);
        end
        
% @Breif    : Crop the Image from a desired pixel
% @param    : length, width ,startpos, I_original
% @ret      : I_cropped       
        function I_cropped = h_crop(length, width ,startpos, I_original)
          I_cropped = uint8(zeros(width,length,3));
            for i = 1:width
                for j = 1:length
                    I_cropped(i,j,1) = I_original(i+startpos,j+startpos,1);
                    I_cropped(i,j,2) = I_original(i+startpos,j+startpos,2);
                    I_cropped(i,j,3) = I_original(i+startpos,j+startpos,3);
                end
            end
        end
        
% @Breif    : Encode the image with generated huffman code
% @param    : I_quantized, width, length, qantization_levels,dict
% @ret      : I_encoded       
        function I_encoded = h_encode(I_quantized, width, length, ...
                                            qantization_levels,dict)
                  
            for i = 1:width
                for j = 1:length
                    for k = 1:qantization_levels
                        if I_quantized(i,j) < dict{1,1} 
                             I_encoded{i,j} = dict{1,2};
                        elseif I_quantized(i,j) >= dict{k,1} && I_quantized(i,j) ...
                                         <= dict{k+1,1}
                             I_encoded{i,j} = dict{k,2};
                        end
                    end
                end
            end
        end
        
% @Breif    : Decode the image with generated huffman code
% @param    : I_encoded, width, length, dict,thresh,quantization_levels
% @ret      : I_decoded             
        function I_decoded = h_decode(I_encoded, width, length, ...
                                            dict,thresh,quantization_levels)
            for i = 1:width
                for j = 1:length
                    for k = 1:quantization_levels
                        if size(I_encoded{i,j}) == size(dict{k,2})
                            if I_encoded{i,j} == dict{k,2}
                             I_decoded(i,j) =uint8(thresh(k));
                            end                          
                        end
                    end
                end
            end
        end
        
% @Breif    : Gives the Entropy
% @param    : symbol_probability,thresh
% @ret      : entropy                   
        function entropy = h_entropy(Image)
                % calculate histogram counts
                P = imhist(Image(:));

                % remove zero entries in p 
                P(P==0) = [];

                % normalize p so that sum(p) is one.
                P = P ./ numel(Image);

            entropy = -sum(P.*log2(P));
        end
% @Breif    : Find the PSNR 
% @param    : symbol_probability,thresh
% @ret      : entropy   
        function psnr = h_PSNR(Original, Decoded)
            im2 = double(Original);
            im1 = double(Decoded);
            mse = sum((im1(:)-im2(:)).^2) / numel(size(im1));       
            psnr = 10*log10(255*255/mse);
        end
        
% @Breif    : Gives the compression ratio 
% @param    : original, encoded
% @ret      : comp_ratio         
        function comp_ratio = h_compratio(original, encoded)
            original_size = size(original,1)*size(original,2)* 8;
            encoded_size = 0;
            for i = 1:size(original,1)
                for j = 1:size(original,2)
                     encoded_size = encoded_size + size(encoded{i,j},2);
                end
            end
            comp_ratio = encoded_size/original_size;
        end
        
     end
end