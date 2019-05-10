% Examples for computing and using Identity Domains to infer personality
% traits from behavioral data
%
% Oren Forkosh, May 2018:  oren.forkosh@gmail.com
%
if verLessThan('matlab','9.2')
    warning('this code has been tested on Matlab 2017a')
end
runtime = tic;

%% load data
load('data_table.mat')

%% Set parameters
%  - number of IDs to compute
opt.nIDs    = 4; 
%  - behaviors normalization parameters
opt.behaviors = {'NormalizeBy', @(t) [ Auxiliary.getargout(3, @unique, t.GroupType), t.Day], 'IgnoreNan', true};
%  - choose subset of groups to be used for training (optional)
opt.trainset   = @(train) train.ConditionID <= 1 & train.Day <= 4;
%  - name of dim-reduction algorithm to use
opt.alg     = @IdentityDomains.GroupLDA;
%  - parameters for dim-reduction algorithm
opt.algaux  = {'IsGroup', false, 'Rotation', 'lda', 'Algorithm', @DimReduction.LDA, 'RelativeTo', [], 'AlgorithmParams', {'Gamma', 1}};

%% Compute and display Identity-Domains

[normalized, behaviors_list] = IdentityDomains.Normalize(data_table, opt.behaviors{:}); 
[w, ~, y] = opt.alg(normalized, opt.nIDs, behaviors_list, 'TrainMap', opt.trainset(normalized), opt.algaux{:});

factor = [1 -1 1 -1]; % choose the positive direction of the axis
w = bsxfun(@times, factor, w);
x = y * w;

figure(1); clf;
IdentityDomains.ShowBase(x, y, behaviors_list);

%% Show precomputed Identity-Domains in Personality space
[x, w, y, arc] = IdentityDomains.ComputeIDs(data_table);
X = Auxiliary.accumrows(data_table.MouseNumber, x, @mean, nan);

% - assign color to each point in pareto space according to distance from
%   archetype
arcmap = [245 165 3; 54 177 191; 242 56 90] / 255;
arcdist = exp(-pdist2(arc, X(:, [1 2])));
arcdist = bsxfun(@rdivide, arcdist, max(arcdist));
pointcolor = zeros(size(arcdist, 2), 3);
for c=1:size(arcdist, 2)
    pointcolor(c, :) = Colors.Mix(arcmap, 'weights', (arcdist(:, c) / sum(arcdist(:, c)))');
end

% - plot Identity Domains in personality space
figure(2); clf;
for mn=Auxiliary.torow(unique(data_table.MouseNumber))
    curr = x(data_table.MouseNumber == mn, :);
    p1 = curr(:, 1);
    p2 = curr(:, 2);
    currcolor = pointcolor(mn, :);
    plot(mean(p1), mean(p2), 'o', 'MarkerFaceColor', currcolor *.8, 'MArkerSize', 10, 'MarkerEdgeColor', currcolor);
    hold on;
    try
        k = convhull(double(p1), double(p2));
        Auxiliary.Polygon(p1(k), p2(k), currcolor, 'FaceAlpha', 1, 'EdgeColor',currcolor*.8);
    catch
    end
    plot(mean(p1), mean(p2), 'o', 'MarkerFaceColor', currcolor *.8, 'MArkerSize', 5, 'MarkerEdgeColor', currcolor);
end
ch = convhull(arc);
patch(arc(ch, 1), arc(ch, 2), 'w', 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .3)
for i=1:size(arc, 1)
    plot(arc(i, 1), arc(i, 2), 'o', 'MarkerFaceColor', arcmap(i, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
end
hold off

% - making things prettier...
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'XTick', [-2 2], 'YTick', [-2 2], 'XMinorTick', 'off')
range = [-4.8 4.8];
xlim(range)
ylim(range)
Auxiliary.Text('ID1', 'location', 'e', 'Color', Colors.PrettyBlue)
Auxiliary.Text('ID2', 'location', 'n', 'Color', Colors.PrettyGreen)
grid minor,axis square,box off

%%
fprintf('# execution took %f secs\n', toc(runtime))