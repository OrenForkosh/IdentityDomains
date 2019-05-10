% Collection of Dimensionality reduction methods
%
% Oren Forkosh, May 2018:  oren.forkosh@gmail.com
%
classdef DimReduction
    % DimReduction Collection of Dimensionality reduction methods
    methods (Static = true)
        function [W, e] = PCA(x, n)
            % PCA Principal component analysis
            %   [W, e] = PCA(x, n) computes the PCA for the data points in
            %   x. Returns the optimal projection W of
            %   dimension n, and the corresponding eigenvalues in e.
            x = bsxfun(@minus, x, mean(x));
            [W, e] = eig(cov(x));
            e(isnan(e)) = 0;
            [e, order] = sort(diag(e), 'descend');
            W = W(:, order(1:n));
        end

        function W = Whiten(x)
            % WHITEN whiten the data making the covariance matrix diagonal
            % aka ZCA or Mahalanobis whitening
            W = inv(sqrtm(cov(x)));
        end

        function W = WhitenPCACor(x)
            % see Kessy, Agnan, Alex Lewin, and Korbinian Strimmer. 
            % "Optimal Whitening and Decorrelation." arXiv:1512.00809 [stat], 
            % December 2, 2015. http://arxiv.org/abs/1512.00809.
            Cor = corr(x);
            Cov = cov(x);
            [G,S] = svd(Cor);
            Vinv = diag(1./sqrt(diag(Cov)));
            W = inv(sqrtm(S)) * G' * Vinv; %#ok<MINV>
        end
        
        function W = WhitenCATCAR(x)
            % see Kessy, Agnan, Alex Lewin, and Korbinian Strimmer. 
            % "Optimal Whitening and Decorrelation." arXiv:1512.00809 [stat], 
            % December 2, 2015. http://arxiv.org/abs/1512.00809.
            %%
            Cov = cov(x);
            Cor = corr(x);
            %V = diag(sqrt(diag(Cov)));
            Vinv = diag(1./sqrt(diag(Cov)));
            P = sqrtm(Cor);
            W = P \ Vinv;
        end
        
        function [d, D] = FisherRaoCriterion(x, labels)
            % measure variability between labels using Fisher-Rao's criterion
            D = DimReduction.ScatterWithin(x, labels)/DimReduction.ScatterBetween(x, 1);
            d = trace(D);
        end
        
        function [W, e, p] = TraceLDA(x, labels, n, varargin)
            % Ngo, Thanh T., Mohammed Bellalij, and Yousef Saad. "The trace
            % ratio optimization problem for dimensionality reduction." 
            % SIAM Journal on Matrix Analysis and Applications 31, no. 5 
            % (2010): 2950-2971.
            %%
            parser = inputParser;
            parser.addOptional('Rotation', '');
            parser.addOptional('Thresh', 1e-7, @isnumeric);
            parser.addOptional('Groups', []);
            parser.addOptional('MaxIters', 250);
            parser.addOptional('Gamma', 0);
            parser.parse(varargin{:});
            opt = parser.Results;
            %%
            if ~isempty(opt.Rotation) && strcmpi(opt.Rotation, 'lda')
                g = n + 1;
            else
                g = n;
            end
            [W, e] = DimReduction.LDA(x, labels, g);
            Sw = DimReduction.ScatterWithin(x, labels);
            if isempty(opt.Groups)
                Sb = DimReduction.ScatterBetween(x, Sw);
            else
                St = DimReduction.ScatterWithin(x, opt.Groups);
                Sb = St - Sw;
            end
            if opt.Gamma > 0
                Sw = Sw + opt.Gamma * eye(size(Sw));
            end
            %
            prevrho = -inf;
            for iter=1:opt.MaxIters
                %%
                rho = trace(W' * Sb * W) / trace(W' * Sw * W);
                if rho - prevrho < opt.Thresh
                    break;
                end
                prevrho = rho;
                G = Sb - rho * Sw;
                [W, e] = eig(G);
                [e, order] = sort(diag(e), 'descend');
                W = W(:, order(1:g));
                e = e(1:g);
            end
            %%
            if n > 1
                if ~isempty(opt.Rotation)
                    X = x * W;
                    fullX = X;
                    X = Auxiliary.accumrows(labels, X, @mean);
                    resort = true;
                else
                    resort = false;
                end

                switch lower(opt.Rotation)
                    case 'lda'
                        Wlda = DimReduction.LDA(fullX, labels, n, 'Groups', opt.Groups);
                        %Wlda = mapping.M;
                        W = W * Wlda;
                        %W = bsxfun(@rdivide, W, sqrt(sum(W.^2)));
                        W = bsxfun(@rdivide, W, sqrt(diag(W' * Sw * W)'));
                        resort = true;
                    case 'pca'
                        Wpca = pca(X);
                        W = W * Wpca;
                        resort = false;
                    case 'ica'
                        epsilon = 0.0001;
                        Wica = [];
                        while size(Wica, 1) < n
                            [~, ~, Wica] = fastica(X', 'verbose', 'off', 'epsilon', epsilon);
                            epsilon = epsilon * 2;
                        end
                        W = W * Wica';
                        
                    case 'catcar'
                        Wcatcar = DimReduction.WhitenCATCAR(X);
                        W = W * Wcatcar;

                    case 'zca'
                        Wzca = DimReduction.Whiten(X);
                        W = W * Wzca;

                    case 'pcacor'
                        Wpcacor = DimReduction.WhitenPCACor(X);
                        W = W * Wpcacor;
                    case {'varimax', 'parsimax', 'promax'}
                        W = rotatefactors(W, 'Method', opt.Rotation);
                        % Put rotated loadings in standard form, again, and make T match
                        [~, ord] = sort(sum(W.^2));
                        W = W(:,fliplr(ord));
                        signs = repmat(sign(sum(W)),size(W, 1),1);
                        W = W .* signs;
                        resort = true;
                    case ''
                    otherwise
                end
                if resort
                    y = x * W;
                    [e, order] = sort(diag(DimReduction.ScatterBetween(y, 1)) ./ diag(DimReduction.ScatterWithin(y, labels)), 'descend');
                    W = W(:, order);
                end
            end
            %%
            if g ~= n
                W = W(:,1:n);
                e = e(1:n);
            end
            
            %% p-value
            p = nan;
        end

        function [W, e] = TraceLDAReg(x, labels, n, varargin)
            % Ngo, Thanh T., Mohammed Bellalij, and Yousef Saad. "The trace
            % ratio optimization problem for dimensionality reduction." 
            % SIAM Journal on Matrix Analysis and Applications 31, no. 5 
            % (2010): 2950-2971.
            %%
            p = inputParser;
            p.addOptional('Thresh', 1e-7, @isnumeric);
            p.addOptional('Groups', []);
            p.addOptional('MaxIters', 250);
            p.addOptional('Gamma', 0);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            g = n; %size(x, 2)-1;
            [W, e] = DimReduction.LDA(x, labels, g);
            Sw = DimReduction.ScatterWithin(x, labels);
            if isempty(opt.Groups)
                Sb = DimReduction.ScatterBetween(x, Sw);
            else
                St = DimReduction.ScatterWithin(x, opt.Groups);
                Sb = St - Sw;
            end
            if opt.Gamma > 0
                Sw = Sw + opt.Gamma * eye(size(Sw));
            end
            %
            prevrho = -inf;
            for iter=1:opt.MaxIters
                %%
                rho = trace(W' * Sb * W) / trace(W' * Sw * W);
                if rho - prevrho < opt.Thresh
                    break;
                end
                prevrho = rho;
                G = Sb - rho * Sw + 10 * Sw;
                [W, e] = eig(G);
                [e, order] = sort(diag(e), 'descend');
                W = W(:, order(1:g));
                e = e(1:g);
            end
            %%
            if g ~= n
                W = W(:,1:n);
                e = e(1:n);
            end
        end
        
        
        function [W, e] = GroupTraceLDA(x, labels, n, varargin)
            % Ngo, Thanh T., Mohammed Bellalij, and Yousef Saad. "The trace
            % ratio optimization problem for dimensionality reduction." 
            % SIAM Journal on Matrix Analysis and Applications 31, no. 5 
            % (2010): 2950-2971.
            %%
            p = inputParser;
            p.addOptional('Thresh', 1e-7, @isnumeric);
            p.addOptional('Groups', []);
            p.addOptional('MaxIters', 250);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            g = size(x, 2)-1;
            [W, e] = DimReduction.LDA(x, labels, g);
            [~, ~, Sw] = DimReduction.ScatterWithin(x, labels);
            if isempty(opt.Groups)
                error('list of groups is empty');
            else
                [~, ~, St] = DimReduction.ScatterWithin(x, opt.Groups);
                Sb = St(:, :, Auxiliary.accumrows(labels, opt.Groups, @mode)) - Sw;
            end
            %
            pG = -inf;
            for iter=1:opt.MaxIters
                %%
                G = 0;
                for i=1:size(Sw, 3)
                    rho = trace(W' * Sb(:, :, i) * W) / trace(W' * Sw(:, :, i) * W);
                    G = G + Sb(:, :, i) - rho * Sw(:, :, i);
                end
                if norm(G-pG,'fro') < opt.Thresh
                    break;
                end
                pG = G;
                [W, e] = eig(G);
                [e, order] = sort(diag(e), 'descend');
                W = W(:, order(1:g));
                e = e(1:g);
            end
            [e, order] = sort(diag(e), 'descend');
            W = W(:, order(1:n));
            e = e(1:n);
        end
        
        function [W, r] = CanonCorr(x, y, n, thresh)
            % CanonCorr Canonial correlation analysis
            %   [W, r] = CanonCorr(x, y, n) computes the CCA for the data points in x
            %   and data points in y. Returns the optimal projection W of
            %   dimension n, and the corresponding correlations in 4.
            %   [W, r] = CanonCorr(x, y, n, thresh) stops the correlation
            %   computation if the correlation is lower then the threshold.
            %%
            if nargin < 4
                thresh = -1;
            end
            if nargin < 3 || isempty(n)
                n = size(y, 2);
            end
            B = eye(size(y, 2));
            W2 = zeros(size(y, 2), n);
            W1 = zeros(1, n);
            r = zeros(1, n);
            Y = y;
            for i=1:n
                %%
                [w1, w2, r(i)] = canoncorr(x, Y); 
                if r(i) < thresh
                    r = r(1:i-1);
                    W2 = W2(:, 1:i-1);
                    W1 = W1(1:i-1);
                    break;
                end
                % w = w / norm(w);
                W2(:, i) = B * w2;
                W1(i) = w1;
                B = null(W2');
                Y = y * B;
            end
            %%
            W = [W1; W2];
        end
        
        function K = Kernel(u,v,varargin)
            % computes the Kernel of two vectors
            p = inputParser;
            p.addOptional('type', 'mlp'); % 'rbf', 'poly', 'quadratic', 'mlp'
            p.addOptional('sigma', 10); % for rbf - gaussian std
            p.addOptional('order', 2); % for poly - polynomial order
            p.addOptional('P1',   1); % for mlp
            p.addOptional('P2',   0); % for mlp
            p.parse(varargin{:});
            opt = p.Results;
            %%
            switch opt.type
                case 'linear'
                    K = (u*v');
                case 'quadratic'
                    dotproduct = (u*v');
                    K = dotproduct.*(1 + dotproduct);
                case 'rbf'
                    K = exp(-(1/(2*opt.sigma^2))*(repmat(sqrt(sum(u.^2,2).^2),1,size(v,1))...
                        -2*(u*v')+repmat(sqrt(sum(v.^2,2)'.^2),size(u,1),1)));
                case 'poly'
                    dotproduct = (u*v');
                    K = dotproduct;
                    for i = 2:opt.order
                        K = K.*(1 + dotproduct);
                    end
                case 'mlp'
                    K = tanh(opt.P1*(u*v')+opt.P2);
            end
        end
        
        function [Kernel, e, W] = KernelLDA(X, labels, ndim, varargin)
            % Computes Kernel LDA
            %%
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('Groups', []);  
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            if isempty(opt.Groups)
                opt.Groups = ones(size(labels));
                ng = 1;
            else
                [~, ~, opt.Groups] = unique(opt.Groups);
                ng = max(opt.Groups);
            end
            
            %%
            % ensures zero mean
            X = bsxfun(@minus, X, mean(X));
            
            % Make sure labels are nice
            [~, ~, labels] = unique(labels, 'rows');
            
            %%
            M = 0;
            N = 0;
            nlabels = 0;
            for g=1:ng
                map = opt.Groups == g;
                [~, ~, local] = unique(labels(map));
                nlocal = max(local);
                %%
                % Compute kernel matrix
                p = Auxiliary.struct2cell(p.Unmatched);
                K = DimReduction.Kernel(X(map, :), X, p{:});
                
                %%
                Mj = Auxiliary.accumrows(local, K, @mean);
                lj = histc(local, 1:nlocal);
                Ms = mean(K);
                m = bsxfun(@minus, Mj, Ms);
                m = bsxfun(@times, bsxfun(@minus, Mj, Ms), lj)' * m;
                n = (K - Mj(local, :))' * K;
                
                %%
                M = M + m * nlocal;
                N = N + n * nlocal;
                %%
                nlabels = nlabels + nlocal;
            end
            M = M / nlabels;
            N = N / nlabels;
            
            %%
            % Compute the eigenvalues and eigenvectors of (N^-1)M
            % We would like to do something like:
            %    [v,ed] = eig(B,W);
            % but in order to insure the condition v'*N*v=I we need:
            epsilon = 0;
            pp = 1;
            while pp > 0
                [R,pp] = chol(N + eye(size(N, 1)) * epsilon);
                epsilon = max(epsilon * 10, 1e-8);
            end
            S = R' \ M / R;
            S = (S+S')/2;     % remove asymmetry caused by roundoff
            [vv,ed] = eig(S);
            W = R\vv;
            
            [e,ei] = sort(diag(ed), 'descend');            % put in descending order
            W = W(:, ei(1:ndim));
            e = e(1:ndim);
            
            Kernel = @(q) DimReduction.Kernel(X, q, p{:})' * W;
        end
        
        function [W, e, s] = LDA(x, labels, n, varargin)
            % LDA Linear discriminant analysis
            %   [W,e] = lda(x, labels, n) computes the LDA for the data points in x
            %   using the class labels. Returns the optimal projection W of
            %   dimension n, and the corresponding eigenvalues in e.
            %%
            p = inputParser;
            p.addOptional('Groups', []);  
            p.addOptional('Whiten', false);
            p.addOptional('Fix', false);
            p.addOptional('Gamma', []);
            p.parse(varargin{:});
            opt = p.Results;
            
            %% whitening
            if opt.Whiten
                P = sqrtm(cov(x));
                x = x / P;
            end
            %%
            
            x = bsxfun(@minus, x, mean(x)); % ensures zero mean
            [Sw, Sb] = DimReduction.Scatter(x, labels, opt.Groups, opt.Gamma);
            
            s.Sb = Sb;
            s.Sw = Sw;
            
            [R,p] = chol(Sw);
            epsilon = 1e-8;
            while (p > 0)
                [R,p] = chol(Sw + eye(size(Sw, 1)) * epsilon);
                epsilon = epsilon * 10;
                warning('matrix singular; adding normalization');
            end
            S = R' \ Sb / R;
            S = (S+S')/2;     % remove asymmetry caused by roundoff
            [vv,e] = eig(S);
            W = R\vv;
            
            if opt.Fix
                y = x * W;
                [e, order] = sort(diag(DimReduction.ScatterBetween(y, 1)) ./ diag(DimReduction.ScatterWithin(y, labels)), 'descend');
                W = W(:, order(1:n));
                W = bsxfun(@rdivide, W, sqrt(sum(W.^2)));
            else
                % sort by eigenvalue
                e(isnan(e)) = 0;
                [e, order] = sort(diag(e), 'descend');
                W = W(:, order(1:n));
                e = e(1:n);
                
            end
            if opt.Whiten
                W = P * W;
            end
        end
        
        function [W, e] = GroupPCA(x, labels, n, varargin)
            % GroupPCA Principal component analysis
            %%
            [~, ~, ic] = unique(labels);
            x = Auxiliary.accumrows(ic, x, @mean);
            [W, e] = DimReduction.PCA(x, n);
        end
        
        function Sw = ScatterWithin(x, labels)            
            % ScatterWithin Computes the scatter within classes for a
            % multiclass LDA model
            %   [Sw, Mw] = ScatterBetweeen(x, labels) computes the scatter of the data
            %   x using the class labels. Mw is the matrix of class
            %   centroids.
            %%
            % fix labels
            [classes, ~, labels] = unique(labels);
            nc = length(classes);
            
            % Now compute the Within sum of squares matrix
            Sw = zeros(size(x,2));
            for j=1:nc
                r = find(labels == j);
                nr = length(r);
                if (nr > 1)
                    z = x(r,:);
                    xm = nanmean(z);
                    z = z - xm(ones(nr,1),:);
                    
                    Sw = Sw + z' * z; %nancov(z, 1) * size(z,1);
                end
            end
        end

        function [Sb, St] = ScatterBetween(x, varargin)
            % ScatterBetweeen Computes the scatter between classes for a
            % multiclass LDA model
            %   ScatterBetweeen(x, Sw) computes the scatter of x using the
            %   scatter-withing matrix Sw
            %   ScatterBetweeen(x, labels) computes the scatter by first
            %   computing the scatter-within matrix using the lables (less
            %   efficient)
            
            if nargin < 1
                error('not enough input arguments');
            end
            if isvector(varargin{1}) && size(x, 1) == length(varargin{1})
                labels = varargin{1};
                Sw = DimReduction.ScatterWithin(x, labels);
            else
                Sw = varargin{1};
            end
            % Start by computing the Total sum of squares matrix
            xm = nanmean(x);
            n = size(x,1);
            x = x - xm(ones(n,1),:);       % done with the original x
            St = nancov(x, 1) * size(x,1);
            
            Sb = St - Sw;
        end
        
        function [Sw, Sb] = Scatter(x,labels,groups,gamma)
            % Computed the scatter-within and -between clusters. Used for
            % the computer of the Fisher-Rao criterion (for LDA)
            if nargin < 4
                gamma = [];
            end
            if nargin < 3
                groups = [];
            end
            if ~isempty(groups)
                [~, ~, groups] = unique(groups);
                ng = max(groups);
            else
                ng = 1;
                groups = ones(size(x,1), 1);
            end
            
            Sb = 0;
            Sw = 0;
            nlabels = 0;
            for g=1:ng
                map = groups == g;
                [~, ~, local] = unique(labels(map));
                sw = DimReduction.ScatterWithin(x(map, :), local);
                if ~isempty(gamma)
                    sw = sw + eye(size(sw,1)) * gamma;
                end
                sb = DimReduction.ScatterBetween(x(map, :), sw);
                nlocal = max(local);
                Sw = Sw + sw * nlocal;
                Sb = Sb + sb * nlocal;
                nlabels = nlabels + nlocal;
            end
            Sw = Sw / nlabels;
            Sb = Sb / nlabels;
        end
                
    end
end
