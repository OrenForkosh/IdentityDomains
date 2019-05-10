% Toolbox for computing and using Identity Domains to infer personality
% traits from behavioral data
%
% Oren Forkosh, May 2018:  oren.forkosh@gmail.com
%
classdef IdentityDomains
    methods (Static = true)
        function [x, w, y, arc] = ComputeIDs(data_table, varargin)
            %%
            p = inputParser;
            addOptional(p, 'NormalizeBy', @(t) [ Auxiliary.getargout(3, @unique, t.GroupType), t.Day]); 
            addOptional(p, 'IgnoreNan', true); 
            addOptional(p, 'Features', 'all'); 
            addOptional(p, 'Arcs', []); 
            addOptional(p, 'W', []); 
            p.parse(varargin{:});
            opt = p.Results;
            %%
            props = Auxiliary.cpfield(opt, struct(),  {'IgnoreNan', 'NormalizeBy'});
            props = Auxiliary.struct2cell(props);
            %%
            load('ID');
            if ~isempty(opt.W)
                ID.w = opt.W;
            end
            t = IdentityDomains.Normalize(data_table, props{:});

            y = IdentityDomains.GetData(false, t, size(ID.w, 2), ID.props);
            x = y * ID.w; 
            w = ID.w; 
            %%
            arc = ID.Pareto.arc;
        end
        
        function [normalized, behaviors] = Normalize(data_table, varargin)
            % Normalize the behaviors.
            %   [normalized, behaviors] = Normalize(data_table) normalizes the
            %   brhaviors in 'data_table'. Returns the 'normalized' table 
            %   and a list of the behaviors.
            %
            %   [...] = Normalize(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
            %   specifies optional parameter name/value pairs to control 
            %   the computation and handling of special data types. (see code)
            %%
            p = inputParser;
            p.addOptional('IgnoreNan', false); % ignore NaN table entries
            p.addOptional('NormalizeBy', 'warp');
            p.addOptional('NormalizeMap', []);
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            behaviors = {'FractionOfTimeOutside', 'VisitsOutsideRate', 'ForagingCorrelation', 'ContactRate', 'ContactRateOutside', 'FractionOfTimeInContactOutside', 'MedianContactDuration', 'MeanContactDuration', 'DiffBetweenApproachesAndChases', 'FractionOfChasesPerContact', 'FractionOfEscapesPerContact', 'FractionOfFollowsPerContact', 'FractionOfBeingFollowedPerContact', 'FractionOfNAChasesPerContact', 'FractionOfNAEscapesPerContact', 'AggressiveChaseRateOutside', 'AggressiveEscapeRateOutside', 'FollowRateOutside', 'BeingFollowedRateOutside', 'NAChaseRateOutside', 'NAEscapeRateOutside', 'AggressiveChaseRate', 'AggressiveEscapeRate', 'FollowRate', 'BeingFollowedRate', 'NAChaseRate', 'NAEscapeRate', 'NumberOfApproachs', 'ApproachRateOutside', 'NumberOfApproachesPerMiceOut', 'FractionOfApproachesPerContact', 'ApproachRate', 'FractionOfBeingApproachedPerContact', 'BeingApproachedRateOutside', 'BeingApproachedRate', 'FractionOfApproachEscapeBehaviorPerAggression', 'Entropy', 'EntropyOutside', 'FractionOfTimeNearFoodOrWater', 'FoodOrWaterPerTimeOutside', 'FractionOfTimeInFeederOutside', 'FractionOfTimeInWaterOutside', 'ProximateVsDistantFood', 'ProximateVsDistantWater', 'FractionOfTimeAtHighPlace', 'HighPlacePerTimeOutside', 'FractionOfTimeInTheOpenOutside', 'FractionOfTimeInSmallNestOutside', 'FractionOfTimeOnRampOutside', 'FractionOfTimeInLabyrinthOutside', 'DistanceFromWallsInOpen', 'DistanceFromNest', 'FractionOfTimeAloneOutside', 'MedianSpeedOutside', 'MeanSpeedOutside', 'TangentialVelocity', 'AngularVelocity', 'DistanceOutside', 'GridEntropy6x6', 'GridMI6x6'};
                        
            %% normalize data table
            normalized = data_table;
            [~, propidx] = ismember(behaviors, normalized.Properties.VariableNames);
            proptable = data_table(:, propidx(propidx>0));
            orig = double(table2array(proptable));
            if isa(opt.NormalizeMap, 'function_handle')
                opt.NormalizeMap = opt.NormalizeMap(normalized);
            end
            if isempty(opt.NormalizeMap)
                opt.NormalizeMap = true(size(orig, 1), 1);
            end
            data = orig(opt.NormalizeMap, :);
            % normalize tables using warping
            if ischar(opt.NormalizeBy)
                switch opt.NormalizeBy
                    case {'warp', ''}
                        normdata = Auxiliary.nwarp(data);
                    case 'znorm'
                        normdata = Auxiliary.znorm(data);
                    case 'none'
                        normdata = data;
                    otherwise
                        normdata = zeros(sum(opt.NormalizeMap), size(proptable, 2));
                        [nrmlzvalue, ~, nrmlz] = unique(data_table.(opt.NormalizeBy)(opt.NormalizeMap));
                        for i=1:length(nrmlzvalue)
                            normdata(nrmlz == i, :) = Auxiliary.nwarp(data(nrmlz == i, :));
                        end
                end
            else
                normdata = zeros(sum(opt.NormalizeMap), size(proptable, 2));
                n = opt.NormalizeBy(data_table);
                [nrmlzvalue, ~, nrmlz] = unique(n(opt.NormalizeMap, :), 'rows');
                for i=1:length(nrmlzvalue)
                    normdata(nrmlz == i, :) = Auxiliary.nwarp(data(nrmlz == i, :));
                end
            end
            if ~isempty(opt.NormalizeMap)
                newdata = zeros(size(orig));
                newdata(opt.NormalizeMap, :) = normdata;
                for i=1:size(data, 2)
                    idx = knnsearch(normdata(:, i), orig(~opt.NormalizeMap, i));
                    newdata(~opt.NormalizeMap, i) = normdata(idx, i);
                end
                normdata = newdata;
            end
            
            if ~opt.IgnoreNan && any(isnan(Auxiliary.torow(table2array(proptable))))
                [i,~] = find(sum(isnan(table2array(proptable)), 1));
                e = proptable.Properties.RowNames(i);
                error('found NaNs in: %s', sprintf('\n\t%s', e{:}))
            end
            normalized(:, propidx(propidx>0)) = array2table(normdata);
        end
                
        function ShowBase(x, y, props)
            %% Show basis
            A = x\y;
            A = A';
            [props, order] = sort(props);
            A = A(order, :);
            if size(A, 1) == length(props) * 2
                subplot(1,2,1);
                Auxiliary.HintonPlot(A(index, :));
                a1 = gca;
                a1.Position(1) = 0.05;
                a1.Position(3) = 0.45;
                title('Individual');
                set(gca, 'XTick', 1:size(A, 2));
                set(gca, 'YTick', 1:length(props));
                subplot(1,2,2);
                Auxiliary.HintonPlot(A(length(props)+index, :));
                a2 = gca;
                a2.Position(1) = 0.55;
                a1.Position(3) = 0.45;
                set(gca, 'YTick', 1:length(props), 'YTickLabel', props(:), 'YAxisLocation', 'right');
                set(gca, 'XTick', 1:size(A, 2));
                title('Group');
            else
                Auxiliary.HintonPlot(A);
                style = {'FontName', 'Noto Sans|Open Sans|Helvetica|Arial'};
                set(gca, 'YTick', 1:length(props), 'YTickLabel', props, 'YAxisLocation', 'right', style{:});
                title({'\bfIdentity-Domains', '\fontsize{10}\rmfull list of behaviors'}, style{:})
            end
        end
        
        function [W, X, Y, e, Sw, Sb] = GroupLDA(fulltable, nf, prop, varargin)
            %%
            p = inputParser;
            p.addOptional('RelativeTo', 'all');
            p.addOptional('TrainGroups', {});
            p.addOptional('IsGroup', false);
            p.addOptional('Rotation', 'lda');
            p.addOptional('NormalizeToGroup', false);
            p.addOptional('GroupScale', 1);
            p.addOptional('Algorithm', @DimReduction.LDA);
            p.addOptional('AlgorithmParams', {});
            p.addOptional('TrainMap', []);
            p.addOptional('PCA', 0);
            %
            p.parse(varargin{:});
            opt = p.Results;
            %%
            alg = opt.Algorithm;
            if isequal(alg,  @DimReduction.LDA)
                aux = opt.AlgorithmParams;
            elseif isequal(alg,  @DimReduction.TraceLDA)
                aux = {'Rotation', opt.Rotation, opt.AlgorithmParams{:}};
            else
                aux = opt.AlgorithmParams;
            end
            
            %%
            [~, ~, labels] = unique(fulltable.MouseNumber);
            [~, ~, groups] = unique(fulltable.GroupNumber);
            [~, ~, epochs] = unique(fulltable.Day);
            [~, ~, grtype] = unique(fulltable.GroupType);
            %% extract individual data from table
            [propmap, propidx] = ismember(prop, fulltable.Properties.VariableNames);
            if length(propidx) < length(prop)
                error('could not find in table the following individual variables: %s', sprintf('''%s'' ', prop{propmap}));
            end
            
            if opt.PCA > 0
                orig = table2array(fulltable(:, propidx));
                m = Auxiliary.accumrows(labels, orig, @mean);
                Wpca = pca(m,  'NumComponents', opt.PCA);
                Yindep = orig * Wpca;
            else
                Yindep = table2array(fulltable(:, propidx));
            end
            %%
            if opt.NormalizeToGroup
                mg = Auxiliary.accumrows(groups, Yindep, @mean);
                mg = mg(groups, :);
                Yindep = Yindep - mg;
            end
            
            %% extract group data from table
            if islogical(opt.IsGroup) && opt.IsGroup
                Ygroupfull = Yindep;
                Ygroup = zeros(size(Yindep));
                for i=1:size(Ygroupfull, 1)
                    if opt.GroupScale == 1
                        Ygroup(i, :) = mean(Ygroupfull(labels ~= labels(i) & groups == groups(i) & epochs == epochs(i), :));
                    else
                        curr = Ygroupfull(labels ~= labels(i) & groups == groups(i) & epochs == epochs(i), :);
                        Ygroup(i, :) = mean(sign(curr) .* (abs(curr).^opt.GroupScale));
                    end
                end
                Y = [Yindep, Ygroup];
            elseif ~islogical(opt.IsGroup)
                Y = [Yindep, opt.IsGroup];
            else
                Y = Yindep;
            end
            if any(isnan(Y(:)))
                warning('NaN in data table');
            end
            
            %%
            if ~isempty(opt.TrainGroups)
                map = ismember(fulltable.GroupType, opt.TrainGroups);
            elseif ~isempty(opt.TrainMap)
                map = opt.TrainMap;
            else
                map = true(size(Y, 1), 1);
            end
            if ischar(opt.RelativeTo)
                switch lower(opt.RelativeTo)
                    case {'groupbyday', ''}
                        gb = Auxiliary.getargout(3, @unique, [groups(map, :), epochs(map, :)], 'rows');
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
                    case 'day'
                        gb = epochs(map, :);
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
                    case 'group'
                        gb = groups(map, :);
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
                    case 'experiment'
                        gb = grtype(map, :);
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
                    case 'all'
                        gb = [];
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, aux{:});
                    otherwise
                        g = fulltable.(opt.RelativeTo);
                        gb = g(map);
                        [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
                end
            elseif isa(opt.RelativeTo,'function_handle')
                g = opt.RelativeTo(fulltable(map, :));
                gb = g(map);
                [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
            else
                gb = opt.RelativeTo;
                [Wlda, e] = alg(Y(map, :), labels(map, :), nf, 'Groups', gb, aux{:});
            end
            if isa(Wlda, 'function_handle')
                X = Wlda(Y);
            else
                % sorting is only needed is using regularization (like 'gamma')
                x = Y * Wlda;
                if ~isequal(alg,  @DimReduction.GroupPCA)
                    [Sw, Sb] = DimReduction.Scatter(x(map, :), labels(map), gb);
                    [e, order] = sort(diag(Sw ./ Sb));
                    Wlda = Wlda(:, order);
                end
                %
                Wlda = bsxfun(@rdivide, Wlda, sqrt(sum(Wlda.^2)));
                X = Y * Wlda;
                %
            end
            
            if opt.PCA > 0
                if islogical(opt.IsGroup) && opt.IsGroup
                    W = [Wpca * Wlda(1:end/2, :); Wpca * Wlda(end/2+1:end, :)];
                elseif ~islogical(opt.IsGroup)
                    W = [Wpca * Wlda(1:end-size(opt.IsGroup, 2), :); Wpca * Wlda(end-size(opt.IsGroup, 2)+1:end, :)];
                else
                    W = Wpca * Wlda;
                end
            else
                W = Wlda;
            end
            
        end
        
        function [Y, nf, labels, groups] = GetData(isgroup, fulltable, nf, props)
            % Create matrix from behavioral table
            [~, ~, labels] = unique(fulltable.MouseNumber);
            [~, ~, groups] = unique(fulltable.GroupNumber);
            [~, ~, epochs] = unique(fulltable.Day);
            
            %% extract individual data from table
            [indepmap, indepidx] = ismember(props, fulltable.Properties.VariableNames);
            if length(indepidx) < length(props)
                error('could not find in table the following individual variables: %s', sprintf('''%s'' ', props{indepmap}));
            end
            Yindep = table2array(fulltable(:, indepidx));
            %% extract group data from table
            if islogical(isgroup) && isgroup
                [groupmap, groupidx] = ismember(props, fulltable.Properties.VariableNames);
                if length(groupidx) < length(props)
                    error('could not find in table the following group variables: %s', sprintf('''%s'' ', props{groupmap}));
                end
                Ygroupfull = table2array(fulltable(:, groupidx));
                Ygroup = zeros(size(Ygroupfull));
                for i=1:size(Ygroupfull, 1)
                    Ygroup(i, :) = mean(Ygroupfull(labels ~= labels(i) & groups == groups(i) & epochs == epochs(i), :));
                end
                %Ygroup = Ygroup - Yindep;
                Y = [Yindep, Ygroup];
            elseif ~islogical(isgroup) && ~isempty(isgroup)
                
                Y = [Yindep, isgroup];
            else
                Y = Yindep;
            end
            if any(isnan(Y(:)))
                warning('NaN in data table');
            end
        end
    end
end