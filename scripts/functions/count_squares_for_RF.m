function [len_w, len_h] = count_squares_for_RF(binary_matrix)
%COUNT_SQUARES_FOR_RF Get the shape of the receptive field
%
%    [LEN_W, LEN_H] = count_squares_for_RF(BINARY_MATRIX)
%    from a BINARY_MATRIX of positive STA responses, calculates the total number of 
%    true (==1) squares. Returns the shape as LEN_W, LEN_H.
%    Expects a matrix that has been generated with bwareaopen, in which true squares
%    are connected.
     w_sum = sum(binary_matrix, 1);
     h_sum = sum(binary_matrix, 2);
    
     len_w = max(w_sum);
     len_h = max(h_sum);
end

