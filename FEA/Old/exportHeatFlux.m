data = [centroids, q_w];
filename = 'heatflux_data.csv';
% Write header
header = {'x', 'y', 'z', 'q_w'};
fid = fopen(filename,'w');
fprintf(fid,'%s,%s,%s,%s\n', header{:});
fclose(fid);
writematrix(data, filename, 'WriteMode', 'append');
disp(['CSV file saved to ', filename])
