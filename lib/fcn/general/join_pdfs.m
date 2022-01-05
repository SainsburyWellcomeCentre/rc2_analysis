function join_pdfs(fnames, join_fname, rm, overwrite)
% JOIN_PDFS Join a set of pdfs
%
%   join_pdfs(SEPARATE_PDFS, JOINED_PDF, REMOVE, OVERWRITE) takes a set of
%   filenames refering to separate pdfs in SEPARATE_PDFS and creates a
%   single pdf, JOINED_PDF. SEPARATE_PDFS should be a cell array of strings
%   with the full paths to the pdfs to join.
%
%   REMOVE is a boolean with whether to remove the original separate pdfs
%   after joining (true). Default is false, to leave the pdfs as they are.
%   OVERWRITE is a boolean, with whether to prompt the user before
%   overwriting an existing file in JOINED_PDF. If set to true, the file
%   will be overwritten without asking. If false (default), the user is
%   prompted.
%
%   The method of merging pdfs was originally here:
%   https://uk.mathworks.com/matlabcentral/fileexchange/89127-merge-pdf-documents 

VariableDefault('rm', false);
VariableDefault('overwrite', false);

if ~overwrite && isfile(join_fname)
    user_ans = input(sprintf('%s already exists. Overwrite (Y/y)?', join_fname), 's');
    if ~strcmpi(user_ans, 'y')
        return
    end
end

mem_set = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;

cellfun(@(f) merger.addSource(f), fnames);
merger.setDestinationFileName(join_fname);
merger.mergeDocuments(mem_set);

if rm
    cellfun(@(f) delete(f), fnames);
end
