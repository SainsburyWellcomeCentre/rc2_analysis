function cameras = format_cameras(camera_dir)
%%FORMAT_CAMERAS
%   CAMERAS = FORMAT_CAMERAS(CAMERA_DIR)
%
%   Currently, simply loads in the .csv files of "global motion" and stores
%   in a structure.

VariableDefault('camera_dir', []);

if isempty(camera_dir)
    fprintf('No camera files selected. Returning empty array\n');
    cameras = [];
    return
end

for sess_i = length(camera_dir) : -1 : 1
    
    if isempty(camera_dir{sess_i})
        cameras(sess_i).camera0 = [];
        cameras(sess_i).camera1 = [];
        continue
    end
    
    fname = fullfile(camera_dir{sess_i}, 'camera0.csv');
    if exist(fname, 'file')
        t = read_camera_csv(fname);
        cameras(sess_i).camera0 = t.e00;
    else
        cameras(sess_i).camera0 = [];
    end
    
    fname = fullfile(camera_dir{sess_i}, 'camera1.csv');
    if exist(fname, 'file')
        t = read_camera_csv(fname);
        cameras(sess_i).camera1 = t.e00;
    else
        cameras(sess_i).camera1 = [];
    end
end