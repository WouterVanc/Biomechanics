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

% Results directories
path_out = fullfile(cd,'Results_load_Comp');
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

% Table info
rownames = [loadvars_perpeak; loadvars_perpeak; loadvars_imp; loadvars_stance];
columnames = {'Loadingvars', 'OA mean','OA std', 'Control mean', 'Control std', 'Pval'};

for q = 1:length(loadvars_perpeak)
    if q>12
        factor = 1000000;
    else
        factor = 1;
    end

    if q<=12
        suffix = '[N]';
    else 
        suffix = '[MPa]';
    end 
    % Peak 1
    mean_oa = mean(loading_results.oa.p1.(loadvars_perpeak{q}));
    std_oa = std(loading_results.oa.p1.(loadvars_perpeak{q}));
    mean_co = mean(loading_results.co.p1.(loadvars_perpeak{q}));
    std_co = std(loading_results.co.p1.(loadvars_perpeak{q}));

    pval_matrix_p1(q,1) = mean_oa/factor;
    pval_matrix_p1(q,2) = std_oa/factor;
    pval_matrix_p1(q,3) = mean_co/factor;
    pval_matrix_p1(q,4) = std_co/factor;

    [p,h] = ranksum(loading_results.oa.p1.(loadvars_perpeak{q}),loading_results.co.p1.(loadvars_perpeak{q}));

    pval_matrix_p1(q,5) = p;

    % Boxplot figure
    y(:,1) = loading_results.oa.p1.(loadvars_perpeak{q});
    y_temp = loading_results.co.p1.(loadvars_perpeak{q});
    y_temp(6:10,1) = NaN;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'OA', 'Control'})
    ylabel([loadvars_perpeak{q} ' ' suffix], Interpreter="none")
    title('Peak 1')
    dim = [.73 .01 .3 .3];
    if p<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(p) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['P1 ' (loadvars_perpeak{q}) '.png']);
    saveas(gcf, fullfile(path_out_p1,plotname))

    close gcf
       
    % Peak 2
    mean_oa2 = mean(loading_results.oa.p2.(loadvars_perpeak{q}));
    std_oa2 = std(loading_results.oa.p2.(loadvars_perpeak{q}));
    mean_co2 = mean(loading_results.co.p2.(loadvars_perpeak{q}));
    std_co2 = std(loading_results.co.p2.(loadvars_perpeak{q}));

    pval_matrix_p2(q,1) = mean_oa2/factor;
    pval_matrix_p2(q,2) = std_oa2/factor;
    pval_matrix_p2(q,3) = mean_co2/factor;
    pval_matrix_p2(q,4) = std_co2/factor;

    [p2,h2] = ranksum(loading_results.oa.p2.(loadvars_perpeak{q}),loading_results.co.p2.(loadvars_perpeak{q}));

    pval_matrix_p2(q,5) = p2;

    % Boxplot figure
    y(:,1) = loading_results.oa.p2.(loadvars_perpeak{q});
    y_temp = loading_results.co.p2.(loadvars_perpeak{q});
    y_temp(6:10,1) = NaN;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'OA', 'Control'})
    ylabel([loadvars_perpeak{q} ' ' suffix], Interpreter="none")
    title('Peak 2')
    dim = [.73 .01 .3 .3];
    if p2<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(p2) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['P2 ' (loadvars_perpeak{q}) '.png']);
    saveas(gcf, fullfile(path_out_p2,plotname))

    close gcf
end

for q = 1:length(loadvars_imp)
    if q>12
        factor = 1000000;
    else
        factor = 1;
    end
    mean_oa = mean(loading_results.oa.imp.(loadvars_imp{q}));
    std_oa = std(loading_results.oa.imp.(loadvars_imp{q}));
    mean_co = mean(loading_results.co.imp.(loadvars_imp{q}));
    std_co = std(loading_results.co.imp.(loadvars_imp{q}));

    pval_matrix_imp(q,1) = mean_oa/factor;
    pval_matrix_imp(q,2) = std_oa/factor;
    pval_matrix_imp(q,3) = mean_co/factor;
    pval_matrix_imp(q,4) = std_co/factor;

    [p3,h] = ranksum(loading_results.oa.imp.(loadvars_imp{q}),loading_results.co.imp.(loadvars_imp{q}));

    pval_matrix_imp(q,5) = p3;

    % Boxplot figure
    y(:,1) = loading_results.oa.imp.(loadvars_imp{q});
    y_temp = loading_results.co.imp.(loadvars_imp{q});
    y_temp(6:10,1) = NaN;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'OA', 'Control'})
    ylabel(['Impulse ' loadvars_imp{q} ' ' suffix], Interpreter="none")
    title('Stance phase')
    dim = [.73 .01 .3 .3];
    if p3<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(p3) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['IMP ' (loadvars_imp{q}) '.png']);
    saveas(gcf, fullfile(path_out_imp,plotname))

    close gcf
end

for q = 1:length(loadvars_stance)
    factor = 1000000;
    mean_oa = mean(loading_results.oa.stance.(loadvars_stance{q}));
    std_oa = std(loading_results.oa.stance.(loadvars_stance{q}));
    mean_co = mean(loading_results.co.stance.(loadvars_stance{q}));
    std_co = std(loading_results.co.stance.(loadvars_stance{q}));
    
    pval_matrix_stance(q,1) = mean_oa/factor;
    pval_matrix_stance(q,2) = std_oa/factor;
    pval_matrix_stance(q,3) = mean_co/factor;
    pval_matrix_stance(q,4) = std_co/factor;
    
    [p4,h] = ranksum(loading_results.oa.stance.(loadvars_stance{q}),loading_results.co.stance.(loadvars_stance{q}));
    
    pval_matrix_stance(q,5) = p4;

    % Boxplot figure
    y(:,1) = loading_results.oa.stance.(loadvars_stance{q});
    y_temp = loading_results.co.stance.(loadvars_stance{q});
    y_temp(6:10,1) = NaN;
    y(:,2) = y_temp;
    
    f1 = figure;
    boxchart(y)
    xticklabels({'OA', 'Control'})
    ylabel([loadvars_stance{q} ' ' suffix], Interpreter="none")
    title('Stance phase')
    dim = [.73 .01 .3 .3];
    if p4<0.05
        symbol = '*';
    else
        symbol = '';
    end
    str = ['p = ' num2str(p4) symbol];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor','none');

    plotname = (['Stance ' (loadvars_stance{q}) '.png']);
    saveas(gcf, fullfile(path_out_stance,plotname))

    close gcf
end

% Combine data
pvalmatrix_temp = [pval_matrix_p1; pval_matrix_p2; pval_matrix_imp; pval_matrix_stance];
pvalmatrix_total = [rownames num2cell(pvalmatrix_temp)];

pvalmatrix_total_ns = pvalmatrix_total;

for col = 6
    for row = 1:length(pvalmatrix_total_ns)
        if ~iscellstr(pvalmatrix_total_ns(row,col)) && cell2mat(pvalmatrix_total_ns(row,col)) > 0.05
            pvalmatrix_total_ns(row,col) = {'n.s'};
        end
    end
end


% Convert to table 
sTable_pval = cell2table(pvalmatrix_total, 'VariableNames',columnames);
sTable_pval_ns = cell2table(pvalmatrix_total_ns, 'VariableNames',columnames);

% Write to Excel file
filename = 'CompLoadVars.xlsx';

writetable(sTable_pval,fullfile(path_out, filename),"Sheet",'OA_vs_Control')
writetable(sTable_pval_ns,fullfile(path_out, filename),"Sheet",'OA_vs_Control_ns')




