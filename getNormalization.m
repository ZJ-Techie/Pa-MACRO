function Y = getNormalization(X, type)

if nargin < 2
    type = 'std';
end

[~,p] = size(X);
Y = X;

if strcmpi(type, 'std')
%     disp('Normlaztion method: standard');
    for i = 1 : p
        Xv = X(:,i);
        Xvn = (Xv-mean(Xv))/std(Xv);
        Y(:,i) = Xvn;
        %     Y(:,j) = Xvn/sqrt(sum(Xvn.^2));
    end
elseif strcmpi(type, 'maxmin')
%     disp('Normlaztion method: maxmin');
    for i = 1 : p
        Xv = X(:,i);
        Xvn = (Xv-min(Xv))/(max(Xv)-min(Xv));
        Y(:,i) = Xvn;
    end
elseif strcmpi(type, 'dmax')
    for i = 1 : p
        Xv = X(:,i);
        Xvn = Xv/max(Xv);
        Y(:,i) = Xvn;
    end
elseif strcmpi(type, 'unit')
    for i = 1 : p
        Xv = X(:,i);
        Xvn = Xv/norm(Xv);
        Y(:,i) = Xvn;
    end
elseif strcmpi(type, 'centered')
%     disp('Normlaztion method: centered');
    for i = 1 : p
        Xv = X(:,i);
        Xvn = Xv - mean(Xv);
        Y(:,i) = Xvn;
    end
elseif strcmpi(type, 'normlize')
%     disp('Normlaztion method: normlized');
    for i = 1:p
        Xv = X(:, i);
        if all(Xv == 0)
            Y(:, i) = Xv;
        else
            Xv = Xv - mean(Xv);
            Xvn = Xv / (norm(Xv)+eps);
            Y(:, i) = Xvn;
        end
    end

%     for i = 1 : p
%         Xv = X(:,i);
%         Xv = Xv - mean(Xv);
%         Xvn = Xv/norm(Xv);
%         Y(:,i) = Xvn;
%     end
end