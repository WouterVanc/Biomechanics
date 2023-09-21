% Compare loadvars for condyles 

% Compare loading variables OA - CO

clear all;
close all;
clc

set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultAxesFontName', 'Arial')

load Final_loadvars.mat

loadvars_perpeak = fieldnames(loading_results.oa.p1);
loadvars_imp = fieldnames(loading_results.oa.imp);
loadvars_stance = fieldnames(loading_results.oa.stance);

c = 0;
for label = [5:8,14,17]
    c=c+1;
    loadvars_perpeak_lat{c,1} = loadvars_perpeak{label};
end

cc = 0;
for label2 = [9:12,15,18]
    cc = cc+1;
    loadvars_perpeak_med{cc,1} = loadvars_perpeak{label2};
end

ccc = 0;
for label3 = [5:8]
    ccc = ccc+1;
    loadvars_imp_lat{ccc,1} = loadvars_imp{label3};
end

cccc = 0;
for label4 = [9:12]
    cccc = cccc+1;
    loadvars_imp_med{cccc,1} = loadvars_imp{label4};
end

loadvars_stance_med = loadvars_stance{3};
loadvars_stance_lat = loadvars_stance{2};

% Results directories
path_out = fullfile(cd,'Results_load_condyles_OA');
path_out_p1 = fullfile(path_out, 'p1');
path_out_p2 = fullfile(path_out, 'p2');
path_out_imp = fullfile(path_out, 'imp');
path_out_stance = fullfile(path_out, 'stance');

if ~exist("path_out", 'dir')
    mkdir(path_out)
end

if ~exist("path_out_p1", 'dir')
    mkdir(path_out_p1)
end

if ~exist("path_out_p2", 'dir')
    mkdir(path_out_p2)
end

if ~exist("path_out_imp", 'dir')
    mkdir(path_out_imp)
end

if ~exist("path_out_stance", 'dir')
    mkdir(path_out_stance)
end

columnames = {'loadingvar', 'Medial mean', 'Medial std', 'Lateral mean', 'Lateral std', 'pval'};
rownames = [loadvars_perpeak; loadvars_perpeak; loadvars_imp; loadvars_stance];

c=0;
for t = [1:4,13,16,19:22,31,34,37:40,49]
    c=c+1;
    rownames2{c,1} = rownames{t};
end

for q = 1:length(loadvars_perpeak_med)
    if q>4
        factor = 1000000;
        suffix = '[MPa]';
    else
        factor = 1;
        suffix = '[N]';
    end
    %%%%% Peak 1 %%%%%
    % mean and std
    mean_lat = mean(loading_results.oa.p1.(loadvars_perpeak_lat{q}));
    std_lat = std(loading_results.oa.p1.(loadvars_perpeak_lat{q}));
    mean_med = mean(loading_results.oa.p1.(loadvars_perpeak_med{q}));
    std_med = std(loading_results.oa.p1.(loadvars_perpeak_med{q}));

    pval_matrix_p1(q,3) = mean_lat/factor;
    pval_matrix_p1(q,4) = std_lat/factor;
    pval_matrix_p1(q,1) = mean_med/factor;
    pval_matrix_p1(q,2) = std_med/factor;
    
    [p,h] = signrank(loading_results.oa.p1.(loadvars_perpeak_med{q}), loading_results.oa.p1.(loadvars_perpeak_lat{q}));
    pval_matrix_p1(q,5) = p;

    % Boxplot figure
    y(:,1) = loading_results.oa.p1.(loadvars_perpeak_med{q})/factor;
    y_temp = loading_results.oa.p1.(loadvars_perpeak_lat{q})/factor;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'Medial', 'Lateral'})
    ylabel([loadvars_perpeak_lat{q}(1:end-4) ' ' suffix], Interpreter="none")
    title('Peak 1')
    dim = [.73 .01 .3 .3];
    if p<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p,4)) symbol];
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['P1 ' (loadvars_perpeak_lat{q}(1:end-4)) '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

    %%%%% Peak 2 %%%%%
    % mean and std
    mean_lat = mean(loading_results.oa.p2.(loadvars_perpeak_lat{q}));
    std_lat = std(loading_results.oa.p2.(loadvars_perpeak_lat{q}));
    mean_med = mean(loading_results.oa.p2.(loadvars_perpeak_med{q}));
    std_med = std(loading_results.oa.p2.(loadvars_perpeak_med{q}));

    pval_matrix_p2(q,1) = mean_med/factor;
    pval_matrix_p2(q,2) = std_med/factor;
    pval_matrix_p2(q,3) = mean_lat/factor;
    pval_matrix_p2(q,4) = std_lat/factor;

    [p2,h2] = signrank(loading_results.oa.p2.(loadvars_perpeak_med{q}), loading_results.oa.p2.(loadvars_perpeak_lat{q}));
    pval_matrix_p2(q,5) = p2;

    % Boxplot figure
    y(:,1) = loading_results.oa.p2.(loadvars_perpeak_med{q})/factor;
    y_temp = loading_results.oa.p2.(loadvars_perpeak_lat{q})/factor;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'Medial', 'Lateral'})
    ylabel([loadvars_perpeak_lat{q}(1:end-4) ' ' suffix], Interpreter="none")
    title('Peak 2')
    dim = [.73 .01 .3 .3];
    if p2<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p2,4)) symbol];
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['P1 ' (loadvars_perpeak_lat{q}(1:end-4)) '.png']);
    saveas(gcf, fullfile(path_out_p2,plotname))

end

for q = 1:length(loadvars_imp_lat)
    suffix = '[N*s]';

    mean_lat = mean(loading_results.oa.imp.(loadvars_imp_lat{q}));
    std_lat = std(loading_results.oa.imp.(loadvars_imp_lat{q}));
    mean_med = mean(loading_results.oa.imp.(loadvars_imp_med{q}));
    std_med = std(loading_results.oa.imp.(loadvars_imp_med{q}));

    pval_matrix_imp(q,1) = mean_med;
    pval_matrix_imp(q,2) = std_med;
    pval_matrix_imp(q,3) = mean_lat;
    pval_matrix_imp(q,4) = std_lat;

    [p3,h] = ranksum(loading_results.oa.imp.(loadvars_imp_lat{q}),loading_results.oa.imp.(loadvars_imp_med{q}));

    pval_matrix_imp(q,5) = p3;
                                                                                    
    % Boxplot figure
    y(:,1) = loading_results.oa.imp.(loadvars_imp_med{q});
    y_temp = loading_results.oa.imp.(loadvars_imp_lat{q});
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'Medial', 'Lateral'})
    ylabel([ 'Impulse ' loadvars_imp_lat{q}(1:end-4) ' ' suffix], Interpreter="none")
    title('Stance phase')
    dim = [.73 .01 .3 .3];
    if p3<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p3,4)) symbol];
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');


    plotname = (['IMP ' (loadvars_imp_lat{q}(1:end-4)) '.png']);
    saveas(gcf, fullfile(path_out_imp,plotname))
end

for q = 1
    suffix = '[MPa]';
    factor = 1000000;
    mean_oa = mean(loading_results.oa.stance.(loadvars_stance_lat));
    std_oa = std(loading_results.oa.stance.(loadvars_stance_lat));
    mean_co = mean(loading_results.oa.stance.(loadvars_stance_med));
    std_co = std(loading_results.oa.stance.(loadvars_stance_med));
    
    pval_matrix_stance(q,1) = mean_oa/factor;
    pval_matrix_stance(q,2) = std_oa/factor;
    pval_matrix_stance(q,3) = mean_co/factor;
    pval_matrix_stance(q,4) = std_co/factor;
    
    [p4,h] = ranksum(loading_results.oa.stance.(loadvars_stance_lat),loading_results.co.stance.(loadvars_stance_med));
    pval_matrix_stance(q,5) = p4;

    % Boxplot figure
    y(:,1) = loading_results.oa.stance.(loadvars_stance_med)/1000000;
    y_temp = loading_results.oa.stance.(loadvars_stance_lat)/1000000;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'Medial', 'Lateral'})
    ylabel([ 'Stance Phase Pressure ' suffix], Interpreter="none")
    title('Stance phase')
    dim = [.73 .01 .3 .3];
    if p4<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(round(p4,4)) symbol];
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['Stance ' (loadvars_stance_lat(1:end-4)) '.png']);
    saveas(gcf, fullfile(path_out_stance,plotname))
end
    
% Combine data
pvalmatrix_temp = [pval_matrix_p1; pval_matrix_p2; pval_matrix_imp; pval_matrix_stance];
pvalmatrix_total = [rownames2 num2cell(pvalmatrix_temp)];

pvalmatrix_total_ns = pvalmatrix_total;

for col = 6
    for row = 1:17
        if ~iscellstr(pvalmatrix_total_ns(row,col)) && cell2mat(pvalmatrix_total_ns(row,col)) > 0.05
            pvalmatrix_total_ns(row,col) = {'n.s'};
        end
    end
end

% Convert to table 
sTable_pval = cell2table(pvalmatrix_total, 'VariableNames',columnames);
sTable_pval_ns = cell2table(pvalmatrix_total_ns, 'VariableNames',columnames);

% Write to Excel file
filename = 'Loadvars_condyles_OA.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'OA') 
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'OA_vs_Control_ns')

close all;






