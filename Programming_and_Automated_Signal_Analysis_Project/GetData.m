function [out] = GetData(age, subject, goodgaits, kinetics, kinematics)

cd('/Users/woutervancaekenberghe/Documents/KU Leuven/Master BReMS/Programming & automated signal analysis')

path = strcat('DataCodeAssignment2122/osimData_',age, ...
    '/',subject,'/RefWalk/');

out = struct;
side = {'IndLeft', 'IndRight'};

% Saving the right variable names

for k = 1:2
    out.(side{k}(4:end)).names = [kinetics kinematics];
end

% First the kinetics data
% Checking if the files exist for the subject

if and(~isempty(kinetics), exist(strcat(path,'GRF/Refwalk_GRF.mot'), 'file'))

    ind_names_GRF = zeros(1,length(kinetics));
    data = ReadMotFile(strcat(path,'GRF/Refwalk_GRF.mot'));
    for i = 1:length(kinetics)
        ind_names_GRF(i) = find(strcmp(data.names, kinetics{i}));
    end

% Distinguish the correct side 

    for j = 1:length(ind_names_GRF)
        
            if kinetics{j}(1) == '1'
                side = {'IndRight'};
            else
                side = {'IndLeft'};
            end
            
% Interpolating the data + finding the mean and standard deviation

        for k = 1:length(side)
            matrix = zeros(101, length(goodgaits.(side{k})));
            for i = 1:length(goodgaits.(side{k}))
                start = goodgaits.(side{k}){i}(1);
                stop = goodgaits.(side{k}){i}(2);
                interpolated_gait = interp1( ...
                    linspace(start, (stop-1), (stop-start+1)).', ...
                    data.data(start:stop, ind_names_GRF(j)), ...
                    linspace(start, (stop-1), 101).');
                matrix(:, i) = interpolated_gait;
                out.(side{k}(4:end)).data{i,j} = interpolated_gait;
            end
            m = mean(matrix,2);
            sd = std(matrix, 0, 2);
        
            out.(side{k}(4:end)).mean{1,j} = m;
            out.(side{k}(4:end)).stdev{1,j} = sd;

        end
    end

end

% Same for the kinematics data

if and(~isempty(kinematics), exist(strcat(path,'KS/KS_Refwalk.mot'), 'file'))

    ind_names_KS = zeros(1,length(kinematics));
    data = ReadMotFile(strcat(path,'KS/KS_Refwalk.mot'));
    for i = 1:length(kinematics)
        ind_names_KS(i) = find(strcmp(data.names, kinematics{i}));
    end


    for j = 1:length(ind_names_KS)

            if kinematics{j}(end) == 'r'
                side = {'IndRight'};
            elseif kinematics{j}(end) == 'l'
                side = {'IndLeft'};
            else
                side = {'IndLeft', 'IndRight'};
            end

        for k = 1:length(side)

            matrix = zeros(101, length(goodgaits.(side{k})));
            for i = 1:length(goodgaits.(side{k}))
                start = floor(goodgaits.(side{k}){i}(1)/5);
                stop = floor(goodgaits.(side{k}){i}(2)/5);
                interpolated_gait = interp1( ...
                    linspace(start, (stop-1), (stop-start+1)).', ...
                    data.data(start:stop, ind_names_KS(j)), ...
                    linspace(start, (stop-1), 101).');
                matrix(:, i) = interpolated_gait;
                out.(side{k}(4:end)).data{i,j+length(kinetics)} = interpolated_gait;
            end
            m = mean(matrix,2);
            sd = std(matrix, 0, 2);
        
            out.(side{k}(4:end)).mean{1,j+length(kinetics)} = m;
            out.(side{k}(4:end)).stdev{1,j+length(kinetics)} = sd;
        end
    end

end

       
