%% 1d single wave animation
clc;
clear;

load('test_mat_file.mat');

cnt = evalin("base", "total");
MovieBuffer(cnt) = struct('cdata',[],'colormap',[]);
for i = 1: cnt
    x_name = "x_" + num2str(i);
    y_name = "y_" + num2str(i);
    plot(evalin("base", x_name), evalin("base", y_name));
    axis([-inf inf -1 2]);
    MovieBuffer(i) = getframe();
end
movie(MovieBuffer, 1, 30);

%% 1d wave animation comparable
clc;
clear;

load('test_mat_file.mat');
name1 = "fft";
name2 = "fd";

cnt = evalin("base", name1 + "_cnt");
MovieBuffer(cnt) = struct('cdata',[],'colormap',[]);
for i = 1: cnt
    x_name = name1 + "_x_" + num2str(i);
    y1_name = name1 + "_y_" + num2str(i);
    y2_name = name2 + "_y_" + num2str(i);
    plot(evalin("base", x_name), evalin("base", y1_name), "-r", ...
         evalin("base", x_name), evalin("base", y2_name), "-b");
    axis([-inf inf -1 1]);
    MovieBuffer(i) = getframe();
end
movie(MovieBuffer, 1, 30);

%% 2d surf
clc;
clear;

load('test_mat_file.mat');
name = "test";

cnt = evalin("base", "total");
MovieBuffer(cnt) = struct('cdata',[],'colormap',[]);
for i = 1: cnt
    x_name = name + "_x_" + num2str(i);
    y_name = name + "_y_" + num2str(i);
    z_name = name + "_z_" + num2str(i);
    X = evalin("base", x_name);
    Y = evalin("base", y_name);
    Z = evalin("base", z_name);
    figure(i);
    s = surf(X, Y, Z);
    colorbar;
    grid on;
    s.EdgeColor = 'none';
end


%% 2d single wave animation
clc;
clear;

load('test_mat_file.mat');
name = "fft";

cnt = evalin("base", "total");
MovieBuffer(cnt) = struct('cdata',[],'colormap',[]);
for i = 1: cnt
    x_name = name + "_x_" + num2str(i);
    y_name = name + "_y_" + num2str(i);
    z_name = name + "_z_" + num2str(i);
    X = evalin("base", x_name);
    Y = evalin("base", y_name);
    Z = evalin("base", z_name);
    s = surf(X, Y, Z);
    axis([-10 10 -10 10 0 1]);
    colorbar;
    s.EdgeColor = 'none';
    MovieBuffer(i) = getframe();
end
movie(MovieBuffer, 1, 10);

%% 1d plot

clc;
clear;

load('test_mat_file.mat');
name = "fft";

cnt = evalin("base", name + "_cnt");
for i = 1: cnt
    x_name = "x_" + num2str(i);
    y_name = "y_" + num2str(i);
    plot(evalin("base", x_name), evalin("base", y_name));
    hold on;
    axis([-inf inf -1 2]);
end

%% 3d scatter
clc;
clear;

load('test_mat_file.mat');
name = "fft";
cnt = evalin("base", name + "_cnt");
transparency = 0.05;
dotsize = 10;

MovieBuffer(cnt) = struct('cdata',[],'colormap',[]);
for i = 1: cnt
    points_name = name + "_3d_" + num2str(i);
    points = evalin("base", points_name);
    % figure(i);
    scatter3(points(1, :), points(2, :), points(3, :), dotsize, points(4, :) .* (-10), 'filled', 'MarkerFaceAlpha', transparency);
    axis([-10 10 -10 10 -10 10]);
    colormap(gca, "winter");
    MovieBuffer(i) = getframe();
end
movie(MovieBuffer, 1, 5);


%% scatter3 test

N = 25;
points = zeros(4, N * N * N);
color = zeros(1, N * N * N);
transparency = 0.2;
len = 0;

for i = 0: N - 1
    for j = 0: N - 1
        for k = 0: N - 1
            index = k * N * N + j * N + i;
            value = exp(-((i - N / 2) ^ 2 + (j - N / 2) ^ 2 + (k - N / 2) ^ 2) * 0.01);
            if value > 0.01
                len = len + 1;
                color(len) = -value;
                points(:, len) = [i; j; k; value];
            end
        end
    end
end

h = scatter3(points(1, :), points(2, :), points(3, :), 5, color, 'filled', 'MarkerFaceAlpha', transparency);
colormap(gca,"hot");

%% 


