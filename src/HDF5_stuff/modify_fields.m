clear variables
close all
home

% MIND THE LOGTRANSFORM! --> conductivities will be taken exp()

prename = 'ftg';
fieldname = 'conductivity';
vertical_profile = [10 -18.421; 18 -6.9077; 20 -13.816];

%fieldname = 'porosity';
%vertical_profile = [20 0.25];

%fieldname = 'storativity';
%vertical_profile = [20 0];


ini = ini2struct(strcat(prename,'.fieldList'));
extensions = ini.grid.extensions(2:end-1);
extensions = strsplit(extensions,' ');
z_ext = str2double(extensions{3});

cells = ini.grid.cells(2:end-1);
cells = strsplit(cells,' ');
z_cells = str2double(cells{3});

discretization = z_ext/z_cells;

z_vec = 0:discretization:z_ext;
property_vec = NaN(length(z_vec(2:end)),1);

layer_counter = 1;
for i = 1:z_cells
    lower = z_vec(i);
    upper = z_vec(i+1);
    property_vec(i) = vertical_profile(layer_counter,2);
    if upper >= vertical_profile(layer_counter,1)
        layer_counter = layer_counter+1;
    end
end
%[z_vec(2:end)',property_vec];

x_max = max(vertical_profile(:,2))+abs(0.1*max(vertical_profile(:,2)));
x_min = min(vertical_profile(:,2))-abs(0.1*min(vertical_profile(:,2)));
if (x_max == 0 )&& (x_min == 0)
    x_max = 1;
    x_min = -1;
end

subplot(1,2,1)
title('Input Layering')
hold on
for i = 1:size(vertical_profile,1)
    if i == 1
    plot([vertical_profile(i,2),vertical_profile(i,2)],[0,vertical_profile(i,1)],'k-','LineWidth',2.5);
    else
    plot([vertical_profile(i,2),vertical_profile(i,2)],[vertical_profile(i-1,1),vertical_profile(i,1)],'k-','LineWidth',2.5);
    end
end
for i = 1:length(z_vec)
    plot([x_min,x_max],[z_vec(i),z_vec(i)],'k-');
end
%if size(vertical_profile,1) > 1
xlim([x_min,x_max]);
%end

subplot(1,2,2)
title('Output Layering')
hold on
for i = 1:length(property_vec)
    plot([property_vec(i),property_vec(i)],[z_vec(i),z_vec(i+1)],'k-','LineWidth',2.5);
end
for i = 1:length(z_vec)
    plot([x_min,x_max],[z_vec(i),z_vec(i)],'k-');
end
%if size(vertical_profile,1) > 1
xlim([x_min,x_max]);
%end

filename = strcat(prename,'.',fieldname,'.stoch.h5');
% %file_id = H5F.open(filename);
% h5disp(filename);
data = h5read(filename,'/stochastic');
% 
% %order of data is x/y/z
if z_cells ~= size(data,3)
    error('Inconsistency in numbers of cells detected. Please check input files!');
end
for i = 1:size(data,3)
    data(:,:,i) = property_vec(i);
end
h5write(filename,'/stochastic',data);
 


fid = fopen(strcat(prename,'.',fieldname,'.trend'),'wt');
fprintf(fid, 'mean = 0');
fclose(fid);
