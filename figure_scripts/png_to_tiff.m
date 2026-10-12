function png_to_tiff(folder)
% png_to_tiff(folder)
% Converts the figures written with "make_figures.py --tiff png" into the PLOS ONE format:
% <name>.png -> <name>.tif (RGB, LZW compression, 300 dpi). The PNG files are kept.
F = dir(fullfile(folder, '*.png'));
for i = 1:numel(F)
    I = imread(fullfile(folder, F(i).name));
    [~, name] = fileparts(F(i).name);
    imwrite(I, fullfile(folder, [name '.tif']), 'tif', 'Compression', 'lzw', 'Resolution', 300);
    fprintf('%s.tif  %d x %d px\n', name, size(I, 2), size(I, 1));
end
end
