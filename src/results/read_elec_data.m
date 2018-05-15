clear variables
close all
home

times = importdata('electrode.times');
no_electrodes = times(1);
times = times(2:end);

data3d = NaN(no_electrodes,no_electrodes,length(times));

for i = 1:no_electrodes
    for j = 1:length(times)
        time = times(j);
        filename = strcat('electrode_',num2str(i),'_',num2str(time),'.data');
        data = importdata(filename);
        data = unique(data,'rows');
        data = sortrows(data);  % sort, so we get a vector of 1->no_electrodes
        data3d(i,data(:,1),j) = data(:,2);
    end
end

