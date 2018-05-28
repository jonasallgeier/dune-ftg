function [handles] = read_potentials_moments(user_input,handles)
    momentfilename = strcat(user_input.basename,user_input.file_potentials_moments,'.moments');
    fileID = fopen(momentfilename,'r');
    % number of electrodes
    tline = fgetl(fileID);
    electrodes = textscan(tline,'%s %d');
    no_electrodes = electrodes{2};
    if no_electrodes ~= handles.no_electrodes
        warning('Number of electrodes in transient and moment data  files is inconsistent.');
    end
    % number of processors
    [~] = fgetl(fileID);
    % read moments
    tline = fgetl(fileID);
    x = textscan(tline(8:end),'%d');
    moments = x{1};
    fclose(fileID);
    if length(moments) <= 1
        error('Implementation currently expects at least the 0th and 1st moment!');
    end
    handles.moments = moments;
    
    data3d_el_moments = NaN(handles.no_electrodes,handles.no_electrodes,length(moments));
    
    textprogressbar('Reading the potential moment data ');
    for j = 1:length(moments)
        moment = moments(j);
        filename = strcat(user_input.basename,user_input.file_potentials_moments,'_',num2str(moment),'.data');
        data = importdata(filename);
        for i = 1:size(data.data,1)
            data3d_el_moments(data.data(i,1),data.data(i,2),j) = data.data(i,3);
        end
        textprogressbar(100*j/(length(moments)+1));
    end
    textprogressbar(100);
    fprintf('\n');
    handles.data3d_el_moments = data3d_el_moments;
end

