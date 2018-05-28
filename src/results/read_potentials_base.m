function [handles] = read_potentials_base(user_input,handles)
    textprogressbar('Reading the base potential data ');
    basedata = NaN(handles.no_electrodes,handles.no_electrodes);
    filename = strcat(user_input.basename,user_input.file_potentials_base,'.data');
    data = importdata(filename);
    for i = 1:size(data.data,1)
        basedata(data.data(i,1),data.data(i,2)) = data.data(i,3);
        textprogressbar(100*i/(size(data.data,1)+1));
    end
    handles.basedata = basedata;
    textprogressbar(100);
    fprintf('\n');
end

