function [Answer, Cancelled] = input_dialog_2(Answer_in)

info = [Answer_in.doPotentialsBase;
        Answer_in.doPotentialsTransient;
        Answer_in.doConcentrationTransient;
        Answer_in.doConcentrationMoments;
        Answer_in.doPotentialsMoments];

Title = 'Define Input Files';

Prompt = {
    'Results data directory/basename','basename',''
    'Number of electrodes','n_el',''
    'Name of electrode coordinates file','file_coords',''
    'Number of combinations','n_comb',''
    'Name of electrode configuration file','file_combs',''
    'Name of potential results (base)','file_potentials_base',''
    'Name of potential results (transient)', 'file_potentials_transient',''
    'Name of concentration results (transient)', 'file_concentration_transient',''
    'Name of concentration results (moments)', 'file_concentration_moments',''
    'Name of potential results (moments)','file_potentials_moments',''
    };

DefAns.basename = './data/';
DefAns.n_el = '144';
DefAns.file_coords = 'electrode.configuration';
DefAns.n_comb = '4650';
DefAns.file_combs = '4_electrodes.txt';
DefAns.file_potentials_base = 'longtestrun_base';
DefAns.file_potentials_transient = 'longtestrun';
DefAns.file_concentration_transient = '_conc';
DefAns.file_concentration_moments = '_m_c';
DefAns.file_potentials_moments = 'moments_el';
Formats = struct([]);

fn = fieldnames(DefAns);
keep = [ones(5,1); info];
for i = 1:length(keep)
    if keep(i)==0
        DefAns = rmfield(DefAns,fn{i});
    end
end
Prompt(~logical(keep),:) = [];

Options.AlignControls = 'on';

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
end
