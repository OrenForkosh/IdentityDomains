% Various colors and tools to handle them
%
% Oren Forkosh, May 2018:  oren.forkosh@gmail.com
%
classdef Colors
    methods (Static)
        function c = Black;       c = Colors.FromHex('000000'); end
        function c = White;       c = Colors.FromHex('FFFFFF'); end
        function c = Red;         c = Colors.FromHex('FF0000'); end
        function c = Yellow;      c = Colors.FromHex('FFFF00'); end
        function c = Blue;        c = Colors.FromHex('0000FF'); end
        function c = PrettyRed;   c = [220 073 089] / 255; end
        function c = PrettyGreen; c = [053 178 087] / 255; end
        function c = PrettyBlue;  c = [000 177 229] / 255; end
        function c = DarkGray;    c = [076 076 076] / 255; end
        function c = LightGray;   c = [178 178 178] / 255; end
        
        function c = DutchTeal;   c = Colors.FromHex('1693A5'); end
        function c = HotPink;     c = Colors.FromHex('FF0066'); end
        function c = Serenity;    c = Colors.FromHex('ADD8C7'); end
        function c = Gold;        c = Colors.FromHex('FBB829'); end
        function c = HauntedMilk; c = Colors.FromHex('CDD7B6'); end
        function c = Slate;       c = Colors.FromHex('556270'); end
        function c = Frogs;       c = Colors.FromHex('C3FF68'); end
        function c = Vanilla;     c = Colors.FromHex('FCFBE3'); end
        function c = Bloons;      c = Colors.FromHex('D31996'); end
        function c = VitaminC;    c = Colors.FromHex('FF9900'); end
        
        function c = RetroBlue;   c = Colors.FromHex('80A2CA'); end
        function c = RetroGreen;  c = Colors.FromHex('C3D254'); end
        function c = RetroRed;    c = Colors.FromHex('E45E57'); end
        function c = RetroOrange; c = Colors.FromHex('EF8444'); end
        function c = RetroPurple; c = Colors.FromHex('A07EBB'); end
        function c = RetroYellow; c = Colors.FromHex('E4D859'); end
        
        function rgb = RGB2HSV(c)
            if size(c, 3) == 3
                rgb = rgb2hsv(c);
            else
                rgb = reshape(rgb2hsv(reshape(c, [size(c,1) 1 3])), [size(c,1) 3]);
            end
        end
        
        function hsv = HSV2RGB(c)
            if size(c, 3) == 3
                hsv = hsv2rgb(c);
            else
                hsv = reshape(hsv2rgb(reshape(c, [size(c,1) 1 3])), [size(c,1) 3]);
            end
        end
        
        function c = Brighten(c, beta)
            if nargin < 2
                beta = .5;
            end
            c = 1 - (beta * (1 - c));
        end

        function c = Mix(varargin)
            if any(strcmpi(varargin, 'weights'))
                stridx = find(strcmpi(varargin, 'weights'));
                weights = varargin{stridx + 1};
                varargin([stridx stridx+1]) = [];
                m = cat(1, varargin{:});
            else
                m = cat(1, varargin{:});
                weights = ones(size(m, 1), 1);
            end
             %cmyk = Colors.RGB2CMYK(m);
             %c = Colors.CMYK2RGB(min(weights(:)' * cmyk , 1));
            M = weights(:)' * Colors.RGB2HSV(m);
            C = Colors.RGB2HSV(sqrt(weights(:)' * m.^2));
            C(3) = max(M(:, 3));
            c = Colors.HSV2RGB(C);
            c = max(min(c, 1), 0);
        end        
    end
end