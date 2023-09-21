function [EMG_channel_White, EMG_channel_Blue] =determine_EMG(ALabels, AnalogSignals);

% get EMG data out of c3d file. 
% Order of Channels: first 8 channel EMG, second 16 channel EMG.
% Input: - Labels of AnalogData --> Output of readc3d-function
%        - AnalogSignals --> Output of readc3d-function
% Output: Matrix with all signals of white EMG and a matrix with signals of
% blue EMG
% 
% By Tuur van der Have, March 2022. 


% determine 8 channel EMG data.
    Awhite_1 = find(strcmp('Voltage.1',ALabels),1);
    Awhite_2 = find(strcmp('Voltage.2',ALabels),1);
    Awhite_3 = find(strcmp('Voltage.3',ALabels),1);
    Awhite_4 = find(strcmp('Voltage.4',ALabels),1);
    Awhite_5 = find(strcmp('Voltage.5',ALabels),1);
    Awhite_6 = find(strcmp('Voltage.6',ALabels),1);
    Awhite_7 = find(strcmp('Voltage.7',ALabels),1);
    Awhite_8 = find(strcmp('Voltage.8',ALabels),1);
    Awhite_9 = find(strcmp('Voltage.9',ALabels),1);
    Awhite_10 = find(strcmp('Voltage.10',ALabels),1);
    Awhite_11 = find(strcmp('Voltage.11',ALabels),1);
    Awhite_12 = find(strcmp('Voltage.12',ALabels),1);
    Awhite_13 = find(strcmp('Voltage.13',ALabels),1);
    Awhite_14 = find(strcmp('Voltage.14',ALabels),1);
    Awhite_15 = find(strcmp('Voltage.15',ALabels),1);
    Awhite_16 = find(strcmp('Voltage.16',ALabels),1);

for i = 1:16;
    eval([ 'C', '=', 'Awhite_', num2str(i), ';']);
    EMG_channel_White(:,i) = AnalogSignals(:,C(1,1));
end


% determine 16 channel EMG data.

    ABlue_1 = find(strcmp('Voltage.1',ALabels),1,'last');
    ABlue_2 = find(strcmp('Voltage.2',ALabels),1,'last');
    ABlue_3 = find(strcmp('Voltage.3',ALabels),1,'last');
    ABlue_4 = find(strcmp('Voltage.4',ALabels),1,'last');
    ABlue_5 = find(strcmp('Voltage.5',ALabels),1,'last');
    ABlue_6 = find(strcmp('Voltage.6',ALabels),1,'last');
    ABlue_7 = find(strcmp('Voltage.7',ALabels),1,'last');
    ABlue_8 = find(strcmp('Voltage.8',ALabels),1,'last');
    ABlue_9 = find(strcmp('Voltage.9',ALabels),1,'last');
    ABlue_10 = find(strcmp('Voltage.10',ALabels),1,'last');
    ABlue_11 = find(strcmp('Voltage.11',ALabels),1,'last');
    ABlue_12 = find(strcmp('Voltage.12',ALabels),1,'last');
    ABlue_13 = find(strcmp('Voltage.13',ALabels),1,'last');
    ABlue_14 = find(strcmp('Voltage.14',ALabels),1,'last');
    ABlue_15 = find(strcmp('Voltage.15',ALabels),1,'last');
    ABlue_16 = find(strcmp('Voltage.16',ALabels),1,'last');
    
 for j = 1:16;
    eval([ 'D', '=', 'ABlue_', num2str(j), ';']);
    EMG_channel_Blue(:,j) = AnalogSignals(:,D(1,1)) ;
 end

end
