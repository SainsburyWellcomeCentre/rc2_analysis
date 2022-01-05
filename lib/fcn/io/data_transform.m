function data = data_transform(data, offset, scale)
% DATA_TRANSFORM Transforms the RC2 data from int16 values to data with a
% meaningful unit (e.g. cm/s)
%
%   VOLTAGE_DATA = data_transform(INT_DATA, OFFSET, SCALE)
%   transforms each channel of data matrix INT_DATA of size 
%   # samples x # channels. OFFSET and SCALE are both 1 x # channel arrays
%   specifying the amount of offset and scale to apply to each channel to
%   convert it into a value with a unit, e.g. cm/s.

data = bsxfun(@minus, data, offset);
data = bsxfun(@times, data, scale);
