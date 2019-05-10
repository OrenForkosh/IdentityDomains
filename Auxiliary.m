% A random collection of auxiliary functions
%
% Oren Forkosh, May 2018:  oren.forkosh@gmail.com
%
classdef Auxiliary
    methods(Static)
        function y = nwarp(x, varargin)
            % warp samples to match a normal distribution
            [rank, N] = Auxiliary.rank(x, varargin{:});
            vals = norminv(linspace(1/N, 1-1/N, N));
            y = vals(rank);
            y = reshape(y, size(x));
        end
        
        function y = znorm(x, varargin)
            cx =  bsxfun(@minus, x, nanmean(x, varargin{:}));
            y = bsxfun(@rdivide, cx, nanstd(x, 0, varargin{:}));
        end
        
        function [r, N] = rank(x, dim, nranks)
            % for each sample returns it's rank
            % for example, for x = [1 4 2 9 5 6], the function will return
            % rank(x) = [1 3 2 6 4 5]
            if nargin > 1
                x = permute(x, [dim, Auxiliary.exclude(ndims(x), dim)]);
            else
                dim = [];
                if size(x, 1) == 1
                    dim = 2;
                    x = permute(x, [dim, Auxiliary.exclude(ndims(x), dim)]);
                end
            end
            [sortx, o] = sort(x, 1);
            N = size(x, 1);
            %%
            seq = 1:N;
            r = zeros(size(x));
            idx = cell(1, ndims(r));
            [idx{:}] = ind2sub(size(r), 1:numel(r));
            idx{1} = o(:)';
            r(sub2ind(size(r), idx{:})) = repmat(seq, 1, numel(o)/N);
            %% map identical values to same rank
            seq = (1:size(x, 1))';
            for i=1:size(x, 2)
                c = cumsum(Auxiliary.padtop(diff(sortx(:, i)) ~= 0, 1, true));
                valid = Auxiliary.padtop(diff(sortx(:, i)) ~= 0, 1, true);
                s = round(Auxiliary.accumrows(cumsum(valid), seq, @mean));
                s = s(c);
                r(o(:, i), i) = s;
            end
            %%
            if ~isempty(dim)
                b = 1:dim-1;
                a = dim+1:ndims(x);
                r = permute(r, [b+1 1 a]);
            end
            %%
            if nargin >= 3
                r = floor((r - 1) / N * nranks) + 1;
            end
        end
        
        function j = exclude(n, i)
            % exclude(n, i) returns the sequence from 1 to n excluding i,
            % i.e. 1 ,.., i-1, i+1 ,..., n
            j = 1:n;
            if nargin > 1
                j(i) = [];
            end
        end
        
        function x = padtop(x, padsize, padval)
            if nargin < 2; padsize = 1; end
            if nargin < 3; padval = 0; end
            x = padarray(x, padsize, padval, 'pre');
        end
        
        function [idx, coord] = argmaxnd(x)
            [~, idx] = max(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), idx);
            coord = cell2mat(a);
        end
        
        function [idx, coord] = argminnd(x)
            [~, idx] = min(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), idx);
            coord = cell2mat(a);
        end
        
        function x = torow(x)
            x = x(:)';
        end
        
        function res = accumrows(subs, val, fun, fillval)
            % ACCUMROWS create matrix by accumulating rows of val specified by subs
            if nargin < 3
                fun = {};
            else
                fun = { fun };
            end
            if nargin < 4
                fillval = 0;
            end
            for i=1:size(val, 2)
                if i==1
                    a = accumarray(subs, val(:, i), [], fun{:}, fillval);
                    res = nan(length(a), size(val, 2));
                    res(:, 1) = a;
                else
                    res(:, i) = accumarray(subs, val(:, i), [], fun{:}, fillval);
                end
            end
        end
        
        function HintonPlot(m, varargin)
            % Generates a Hinton plot (for displaying matrix data)
            %%
            p = inputParser;
            addOptional(p, 'annotate', false, @islogical);
            addOptional(p, 'abs', true, @islogical);
            addOptional(p, 'minmax', false, @islogical);
            addOptional(p, 'mark', false(0), @islogical); % logical matrix showing which entries to highlight
            addOptional(p, 'colormap', [Colors.PrettyBlue; Colors.PrettyRed], @(x) isnumeric(x) && size(x, 2) == 3);
            addOptional(p, 'range', [], @(x) isnumeric(x) && length(x) == 2);
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            cla
            csz = size(opt.colormap, 1);
            if isempty(opt.range)
                opt.range = quantile(m(:), [.025, .975]);
                if opt.abs
                    opt.range = max(abs(opt.range)) * [-1 1];
                end
            end
            r = max(min(m / (opt.range(2)-opt.range(1)), 1), -1);
            cidx = floor((m - opt.range(1)) / (opt.range(2)-opt.range(1)) * (csz)) + 1;
            cidx = max(min(cidx, csz), 1);
            
            [mx, mxij] = Auxiliary.argmaxnd(m);
            [mn, mnij] = Auxiliary.argminnd(m);
            mx = m(mx);
            mn = m(mn);
            for i=1:size(m, 1)
                for j=1:size(m, 2)
                    wh = abs(r(i,j));
                    color = opt.colormap(cidx(i, j), :);
                    if isnan(m(i, j))
                        continue;
                    end
                    if ~isempty(opt.mark) && opt.mark(i, j)
                        Auxiliary.Rect(j-wh/2, i-wh/2, wh, wh, color, 'EdgeColor', 'k');
                    else
                        Auxiliary.Rect(j-wh/2, i-wh/2, wh, wh, color);
                    end
                    hold on;
                    if opt.annotate
                        c = Colors.RGB2HSV(color);
                        c(2) = 1 - c(2);
                        c(3) = 1 - c(3);
                        c = Colors.HSV2RGB(c);
                        text(j, i, sprintf('%3.2g', m(i,j)), 'horizontalalignment', 'center', 'Color', c);
                    end
                end
            end
            
            if opt.minmax
                for i=1:size(m, 1)
                    for j=1:size(m, 2)
                        if all([i,j] == mxij)
                            text(j, i, sprintf('%3.1g', mx), 'horizontalalignment', 'center');
                        elseif all([i,j] == mnij)
                            text(j, i, sprintf('%3.1g', mn), 'horizontalalignment', 'center');
                        end
                    end
                end
            end
            hold off;
            axis equal
            xlim([0, size(m,2)+1]);
            ylim([0, size(m,1)+1]);
            set(gca, 'YDir', 'reverse', 'YTick', 1:size(m,1), 'XTick', 1:size(m,2))
            Auxiliary.Prettify
            box on;
            grid on;
        end
        
        function Text(str, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('margin', 40);
            p.addOptional('location', 'nw', @(x) ismember(x, {'nw', 'sw', 'ne', 'se', 'n', 'w', 's', 'e'}));
            p.parse(varargin{:});
            opt = p.Results;
            %%
            nf = length(fieldnames(p.Unmatched));
            aux(2:2:nf*2) = struct2cell(p.Unmatched);
            aux(1:2:nf*2) = fieldnames(p.Unmatched);
            
            %%
            x = xlim;
            y = ylim;
            d = {'FontName', 'Calibri Light', 'FontSize', 12,  'color', [.3 .3 .3]};
            switch opt.location
                case 'n'
                    pos = [(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 's'
                    pos = [(x(2)-x(1))/opt.margin,  y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', d{:}, aux{:});
                case 'w'
                    pos = [x(1)+(x(2)-x(1))/opt.margin,  -2*(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', d{:}, aux{:});
                case 'e'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, -2*(y(2)-y(1))/opt.margin];
                    %pos = [x(2)-(x(2)-x(1))/opt.margin,  -(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', d{:}, aux{:});
                case 'nw'
                    pos = [x(1)+(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 'sw'
                    pos = [x(1)+(x(2)-x(1))/opt.margin, y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', d{:}, aux{:});
                case 'ne'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 'se'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', d{:}, aux{:});
            end
        end
        
        function Prettify
            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.005 .005] , ...
                'YMinorTick'  , 'off'      , ...
                'XColor'      , [.3 .3 .3], ...
                'YColor'      , [.3 .3 .3], ...
                'XMinorTick'  , 'on'      , ...
                'LineWidth'   , 1         );
            
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
            set(gcf, 'color', [1 1 1]);
        end
        
        function a = getargout(index, func, varargin)
            X = cell(1, index);
            if nargin > 2
                [X{:}] = func(varargin{:});
            else
                [X{:}] = func();
            end
            a = X{index};
        end
        
        function to = cpfield(from, to, names, ignore)
            if nargin == 2
                f = fieldnames(from);
                for i=1:length(f)
                    to = Auxiliary.cpfield(from, to, f{i});
                end
                return;
            end
            if nargin == 3
                ignore = false;
            end
            if ischar(names)
                if Auxiliary.isfield(from, names)
                    to = Auxiliary.setfield(to, names, Auxiliary.getfield(from, names));
                elseif ~ignore
                    throw(MException('Q:nonExistentField', 'Reference to non-existent field ''%s''', names));
                end
            else
                for i=1:length(names)
                    to = Auxiliary.cpfield(from, to, names{i}, ignore);
                end
            end
        end
        
        function val = getfield(var, name)
            fields = textscan(name, '%s', 'Delimiter', '.');
            val = getfield(var, fields{1}{:});
        end
        
        function var = setfield(var, varargin)
            for i=1:2:length(varargin)
                fields = textscan(varargin{i}, '%s', 'Delimiter', '.');
                var = setfield(var, fields{1}{:}, varargin{i+1});
            end
        end
        
        
        function tf = isfield(var, name)
            %% like matlab's isfield but also approves properties and
            %% can support several levels, like: Auxiliary.isfield(s, 'a.b.c.d');
            %%
            c = strsplit(name, '.');
            curr = var;
            tf = true;
            for i=1:length(c)
                if ~isfield(curr, c{i}) && ~isprop(curr, c{i})
                    tf = false;
                    break;
                end
                curr = curr.(c{i});
            end
        end
        
        function c = struct2cell(var)
            % convert a struct to a cell including fieldnames
            a = struct2cell(var);
            a(:, 2) = a(:, 1);
            a(:, 1) = fieldnames(var);
            a = a';
            c = a(:);
        end
        
        function s = IsTableVar(t, v)
            s = any(strcmp(t.Properties.VariableNames, v));
        end
        
        function t = InitColumnIfMissing(t, v, val)
            if ~Auxiliary.IsTableVar(t, v)
                if isnumeric(val)
                    t.(v) = repmat(val, size(t, 1), 1);
                else
                    t.(v) = repmat({val}, size(t, 1), 1);
                end
            end
        end
        
        function Rect(x,y,w,h,c,varargin)
            if ~ishold
                cla reset
            end
            patch([x x+w x+w x], [y y y+h y+h], c, 'EdgeColor', 'none', varargin{:});
        end
        
        function Polygon(x,y,c,varargin)
            if ~ishold
                cla reset
            end
            patch(x,y,c,'EdgeColor', 'none', varargin{:});
        end
        
    end
end
