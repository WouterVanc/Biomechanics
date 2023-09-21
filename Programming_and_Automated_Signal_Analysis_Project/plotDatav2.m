function [] = plotDatav2(age, kinetics, kinematics, CompareGroups, side)

dim = get(0,'Screensize');

% Define titels for the plots and files 

str0 = 'figures/';
str1 = (age);
str2 = 'MeanStdev';
str3 = (side);
figuretitel = (append(str0,str2,'_',str1,'_',str3));
figuretitelv2 = (append(str1,' ',str3));

% Creating a tiled layout

f = figure('Position', [1 1 1280 550]);
t1 = tiledlayout(2,1,'TileSpacing','Compact');
t2 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3.Layout.Tile = 2;
t1.Title.String = figuretitelv2;
t2.YLabel.String = 'Force[N]';
t3.YLabel.String = 'Joint Angle';
Total = [kinetics kinematics];

% Iterating the kinetics on tile and kinematics on another tile

for i = 1:length(kinetics)

    x = linspace(0,100,101);
    y = CompareGroups.(age).MeanStdev(:,i);
    
    nexttile(t2);
    hold on
    title(Total{i},'Interpreter','none')
    plot(x, y,'k',LineWidth=1)
    hold off
end

for i = length(kinetics)+1:length(Total)

    x = linspace(0,100,101);
    y = CompareGroups.(age).MeanStdev(:,i);

    nexttile(t3);
    hold on
    title(Total{i},'Interpreter','none')
    plot(x, y,'k',LineWidth=1)
    hold off
end

savefig(f,figuretitel,"compact")

close(f)

