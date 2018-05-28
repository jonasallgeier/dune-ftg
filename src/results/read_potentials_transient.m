function [handles] = read_potentials_transient(user_input,handles)
timefilename = strcat(user_input.basename,user_input.file_potentials_transient,'.times');
fileID = fopen(timefilename,'r');
% number of electrodes
tline = fgetl(fileID);
electrodes = textscan(tline,'%s %d');
handles.no_electrodes = electrodes{2};
% number of processors
[~] = fgetl(fileID);
% read timestamps
tline = fgetl(fileID);
x = textscan(tline(6:end),'%f');
times = x{1};
fclose(fileID);
handles.times = times;


data3d_el_transient = NaN(handles.no_electrodes,handles.no_electrodes,length(times));

textprogressbar('Reading the transient potential data ');
for j = 1:length(times)
    time = times(j);
    filename = strcat(user_input.basename,user_input.file_potentials_transient,'_',num2str(time),'.data');
    data = importdata(filename);
    for i = 1:size(data.data,1)
        data3d_el_transient(data.data(i,1),data.data(i,2),j) = data.data(i,3);
    end
    textprogressbar(100*j/(length(times)+1));
end
textprogressbar(100);
fprintf('\n');
handles.data3d_el_transient = data3d_el_transient;
end

