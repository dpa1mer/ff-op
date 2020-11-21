function ExportMorse(filename, verts, faces, field)

% Write modified OFF file
nv = size(verts, 1);
nf = size(faces, 1);

writematrix('OFF', filename, 'FileType', 'text');
writematrix([nv nf], filename, 'WriteMode', 'append', 'Delimiter', ' ', 'FileType', 'text');
writematrix([verts field], filename, 'WriteMode', 'append', 'Delimiter', ' ', 'FileType', 'text');
writematrix([3*ones(nf, 1) (faces - 1)], filename, 'WriteMode', 'append', 'Delimiter', ' ', 'FileType', 'text');

end