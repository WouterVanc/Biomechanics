function [] = plotData(age, subject, kinetics, kinematics, data, side)

dim = get(0,'Screensize');

% Define titels for the plots and files

str0 = 'figures/';
str1 = (age);
str2 = (subject);
str3 = (side);
figuretitel = (append(str0,str1,'_',str2,'_',str3));
figuretitelv2 = (append(str1,' ',str2,' ',str3));
Total = [kinetics kinematics];

% Creating a tiled layout

f = figure('Position', [1 1 1280 550]);
t1 = tiledlayout(2,1,'TileSpacing','Compact');
t2 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3.Layout.Tile = 2;
t1.Title.String = figuretitelv2;
t1.Title.Interpreter = 'none';
t1.Title.FontWeight = 'bold';
t1.Title.FontSize = 16;
t2.YLabel.String = 'Force [N]';
t3.YLabel.String = 'Joint Angle [°]';

% Iterating the kinetics on tile and kinematics on another tile

for i = 1:length(kinetics)

    x = linspace(0,100,101);
    m = data.(age).(subject).(side).mean{1,i};
    s = data.(age).(subject).(side).stdev{1,i};

    upper = (m+s).';
    lower = (m-s).';

    xx = [x fliplr(x)].';
    yy = [lower fliplr(upper)];

    nexttile(t2);
    hold on
    title(Total{i},'Interpreter','none')
    plot(x, m,'k',LineWidth=1)
    fill(xx, yy, 'k', FaceAlpha=0.2, EdgeAlpha=0.0)
    hold off
end

for i = length(kinetics)+1:length(Total)

    x = linspace(0,100,101);
    m = data.(age).(subject).(side).mean{1,i};
    s = data.(age).(subject).(side).stdev{1,i};

    upper = (m+s).';
    lower = (m-s).';

    xx = [x fliplr(x)].';
    yy = [lower fliplr(upper)];

    nexttile(t3);
    hold on
    title(Total{i},'Interpreter','none')
    plot(x, m,'k',LineWidth=1)
    fill(xx, yy, 'k', FaceAlpha=0.2, EdgeAlpha=0.0)
    hold off

end

savefig(f,figuretitel,"compact")

close(f)
            





















