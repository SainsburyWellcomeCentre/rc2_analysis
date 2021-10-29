function data = data_transform(data, offset, scale)
% transform each channel of data matrix (written to convert volts
% to cm/s etc.) by subtracting offset and multiplying by scale
% data - n x m matrix (channels along rows)
% offset, scale - n x 1 arrays

data = bsxfun(@minus, data, offset);
data = bsxfun(@times, data, scale);
