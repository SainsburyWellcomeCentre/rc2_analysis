function join_pdfs(fnames, join_fname, rm, overwrite)

VariableDefault('rm', false);
VariableDefault('overwrite', false);

pad_fnames = cellfun(@(x)(['"', x, '"']), fnames, 'uniformoutput', false);
fname_list = strjoin(pad_fnames, ' ');

join_fname = ['"', join_fname, '"'];

this_dir = fileparts(mfilename('fullpath'));
py_file = ['"', fullfile(this_dir, 'join_pdfs.py'), '"'];

if overwrite
    cmd = sprintf('python %s -f -o %s %s', py_file, join_fname, fname_list);
else
    cmd = sprintf('python %s -o %s %s', py_file, join_fname, fname_list);
end
exit_status = system(cmd);

if rm && ~exit_status
    for file_i = 1 : length(fnames)
        delete(fnames{file_i});
    end
end
