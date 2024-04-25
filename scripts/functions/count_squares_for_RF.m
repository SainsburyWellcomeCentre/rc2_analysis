function [len_w, len_h] = count_squares_for_RF(binary_matrix)
%COUNT_SQUARES_FOR_RF Summary of this function goes here
%   Detailed explanation goes here
     w_sum = sum(binary_matrix, 1);
     h_sum = sum(binary_matrix, 2);
    
     len_w = max(w_sum);
     len_h = max(h_sum);
end

