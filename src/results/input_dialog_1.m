function [Answer,Cancelled] = input_dialog_1()

Title = 'Select Data Types To Display';

%%%% SETTING DIALOG OPTIONS
%Options.WindowStyle = 'modal';
Options.Resize = 'on';
%Options.CancelButton = 'on';
%Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration

Prompt = {};
Formats = {};
DefAns = struct([]);

Prompt(1,:) = {['Please select the types of data to display!'],[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
Formats(1,1).span = [1 1]; % item is 1 field x 4 fields


Prompt(2,:) = {'Potentials (base)' 'doPotentialsBase',[]};
Formats(2,1).type = 'check';
DefAns(1).doPotentialsBase = true;


Prompt(3,:) = {' Potentials (transient)' 'doPotentialsTransient',[]};
Formats(3,1).type = 'check';
Formats(3,1).enable = 'off';
DefAns.doPotentialsTransient = true;

Prompt(4,:) = {'Concentration (transient)' 'doConcentrationTransient',[]};
Formats(4,1).type = 'check';
Formats(4,1).enable = 'off';
DefAns.doConcentrationTransient = false;

Prompt(5,:) = {'Concentration (moments)' 'doConcentrationMoments',[]};
Formats(5,1).type = 'check';
Formats(5,1).enable = 'off';
DefAns.doConcentrationMoments = false;

Prompt(6,:) = {'Potentials (moments)' 'doPotentialsMoments',[]};
Formats(6,1).type = 'check';
DefAns.doPotentialsMoments = false;

Formats(1).size = [100, 50];

%warning('The display of concentration data (transient and moments) is not implemented yet.');

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
end
