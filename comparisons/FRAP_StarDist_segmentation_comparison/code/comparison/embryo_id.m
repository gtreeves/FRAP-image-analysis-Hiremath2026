function id = embryo_id(filenum, file)
%EMBRYO_ID Short embryo ID used for all output folders and files.
%   id = embryo_id(filenum, file)  ->  "E<filenum>_<date>", e.g. "E03_2022-10-11"
%   file is the relative CZI path from embryo_list.csv (date folder first).
dt = regexp(string(file), '\d{4}-\d{1,2}-\d{1,2}', 'match', 'once');
id = string(sprintf('E%02d_%s', filenum, dt));
end
