function [Uall, Vall, Wall, Q] = PaMACRO(data, opts)   
    %------------------------------------------
    X = data.X;
    C = data.C;
    unX = data.unX;
    unC = data.unC;
    E = data.E;
    Y = data.Y;
    Z1 = data.Z1;
    Z2 = data.Z2;

    a = size(X, 1);
    b = size(X, 2);
    c = size(C, 2);
    d = size(E, 2);
    I = size(Y, 2);

    X = getNormalization(X, 'normlize');
    C = getNormalization(C, 'normlize');
    Y = getNormalization(Y, 'normlize');
    unX = getNormalization(unX, 'normlize');
    unC = getNormalization(unC, 'normlize');
    E = getNormalization(E, 'normlize');
    Z1 = getNormalization(Z1, 'normlize');

    U = ones(b, 1); 
    U1 = U / norm(U, 2);
    U2 = U / norm(U, 2);
    U3 = U / norm(U, 2);
    
    V = ones(c, 1);
    V1 = V / norm(V, 2);
    V2 = V / norm(V, 2);
    V3 = V / norm(V, 2);
    
    W = ones(I, 1); 
    W1 = W / norm(W, 2);
    W2 = W / norm(W, 2);
    W3 = W / norm(W, 2);  

    max_iter = 30;
    q = 0;
    
    % set parameters
    lamdaU1L1 =  opts.lambda.u1;
    lamdaU2L1 =  opts.lambda.u1;
    lamdaU3L1 =  opts.lambda.u1;

    lamdaUL21 =  opts.lambda.u2;
    lamdaf21 =   opts.lambda.u3;

    lamdaV1L1 =  opts.lambda.v1;
    lamdaV2L1 =  opts.lambda.v1;
    lamdaV3L11 = opts.lambda.v11;
    lamdaVL21 =  opts.lambda.v2;

    lamdaW1L1 =  opts.lambda.w1;
    lamdaW2L1 =  opts.lambda.w1;
    lamdaW3L1 =  opts.lambda.w1;
    lamdaWL21 =  opts.lambda.w2;

    CC = C' * C;
    UCC = unC' * unC;
    
    YY = Y' * Y;
    
    p = b;
    block = 250;
    Xall = X;
    unXall = unX;
    nblock = ceil(p / block);
   
    Q_ini = rand(b, d);
    zz = vectorization(Q_ini);
    zz = zz./norm(zz);
    Q_ini = de_vectorization(zz, b, d);
    Q = Q_ini;
    step = 0.001;

  while (q < max_iter) % default 100 times of iteration

      q = q + 1;
      U1t = [];
      U2t = [];
      U3t = [];
      uq = [];

      sU1 = 0;
      sU2 = 0;
      sU3 = 0;

      U1_old = U1;
      U2_old = U2;
      U3_old = U3;

    for iu = 1 : nblock
        if iu * block <= p
            X = Xall(:, 1 + (iu - 1) * block : iu * block);
            unX = unXall(:, 1 + (iu-1) * block : iu * block);
            sub_U1 = U1(1 + (iu-1) * block : iu * block);
            sub_U2 = U2(1 + (iu-1) * block : iu * block);
            sub_U3 = U3(1 + (iu-1) * block : iu * block);
            Qb = Q(1 + (iu - 1) * block : iu * block, :);
            XX = X' * X;
            UXX = unX' * unX;   
        else
            X = Xall(: , 1 + (iu - 1) * block : end);
            unX = unXall( :, 1 + (iu - 1) * block : end);
            sub_U1 = U1(1 + (iu - 1) * block : end);
            sub_U2 = U2(1 + (iu - 1) * block : end);
            sub_U3 = U3(1 + (iu - 1) * block : end);
            Qb = Q(1 + (iu - 1) * block : end, :);
            XX = X' * X;
            UXX = unX' * unX;   
        end

        % solve u1 u2 u3
        % updata Du
        DU21 = updateDs(sub_U1, sub_U2, sub_U3);
        DU21 = diag(DU21);

        DUFGL21 = updateD_FGL21(sub_U1, sub_U2, sub_U3);
        DUFGL21 = diag(DUFGL21);

        DU1 = updateD(sub_U1);
        DU1 = diag(DU1);
        
        DU2 = updateD(sub_U2);
        DU2 = diag(DU2);
        
        DU3 = updateD(sub_U3);
        DU3 = diag(DU3);

        bsmall = size(X , 2);
        
        % Update U
        sub_U1 = (XX + lamdaU1L1 * DU1 + lamdaf21 * DUFGL21 + lamdaUL21 * DU21 ) \ (X' * (Y * W1 - C * V1 - diag(Xall * Q * E')));
        sub_U2 = (XX + lamdaU2L1 * DU2 + lamdaf21 * DUFGL21 + lamdaUL21 * DU21 ) \ (X' * C * V2);
        sub_U3 = (UXX + lamdaU3L1 * DU3 + lamdaf21 * DUFGL21 + lamdaUL21 * DU21 ) \ (unX' * unC * V3);
      
        U1t = [U1t; sub_U1];
        U2t = [U2t; sub_U2];
        U3t = [U3t; sub_U3];
        
        sU1 = sU1 + sub_U1' * XX * sub_U1;
        sU2 = sU2 + sub_U2' * XX * sub_U2;
        sU3 = sU3 + sub_U3' * UXX * sub_U3;

       % update Q
       m = size(X, 1);
       pb = size(X, 2);
       grad_Q = zeros(pb, d);
        for i = 1 : a
                xi = X(i, :);
                ci = C(i, :);
                ei = E(i, :);
                yWi = Y(i, :) * W1;
                Xtemp = (xi * sub_U1 +  xi * Qb * ei' + ci * V1 - yWi);
                grad_Q = grad_Q + (xi' * Xtemp) * ei;
        end
        
        Qb = Qb - step * grad_Q;
        Qvec = vectorization(Qb);
        Qvec = Qvec./norm(Qvec);
        lambdaQ = 0.01;
        Qs = soft_Threshold (Qvec, lambdaQ);
        Qs = Qs ./norm(Qs);
        Qb = de_vectorization(Qs, pb, d);
        uq = [uq; Qb];
    end

%     % scale U1 U2 U3
%     U1 = U1t./norm(U1,2);
%     U2 = U2t./norm(U2,2);
%     U3 = U3t./norm(U3,2);

    Q = uq;
    XU1 = Xall * U1;
    U1 = U1t./norm(XU1, 2);
    XU2 = Xall * U2;
    U2 = U2t./norm(XU2, 2);
    XU3 = unXall * U3;
    U3 = U3t./norm(XU3, 2);

    % update Dv
    d1 = updateD(V1);
    DV1 = diag(d1);

    d2 = updateD(V2);
    DV2 = diag(d2);

    d3 = updateD(V3);
    DV3 = diag(d3);

    DV = updateDs(V1, V2, V3);
    DsV = diag(DV);  

    V1 = (CC + lamdaV1L1 * DV1 + lamdaVL21 * DsV ) \ (C' * (Y * W1 - Xall * U1 - diag(Xall * Q * E')));
    V2 = (CC + lamdaV2L1 * DV2 + lamdaVL21 * DsV ) \ (C' * Xall * U2);
    V3 = (UCC + lamdaV3L11 * DV3 + lamdaVL21 * DsV) \  (unC' * unXall * U3);

    CV1 = C * V1;
    CV2 = C * V2;
    CV3 = C * V3;

    V1 = V1./norm(CV1, 2);
    V2 = V2./norm(CV2, 2);
    V3 = V3./norm(CV3, 2);

    % update DW
    d1 = updateD(W1);
    DW1 = diag(d1);

    d2 = updateD(W2);
    DW2 = diag(d2);

    d3 = updateD(W3);
    DW3 = diag(d3);

    DW = updateDs(W1, W2, W3);
    DsW = diag(DW);

    lr = 0.01;
    max_iter = 50;

    % Update W Imaging
    W1 = (YY + lamdaW1L1 * DW1 + lamdaWL21 * DsW  ) \ (Y' * (Xall * U1 + diag( Xall * Q * E' ) + C * V1));
    W2 = (YY + lamdaW2L1 * DW2 + lamdaWL21 * DsW  ) \ (Y' * Z1);
    W3 = LR(Y, Z2, lamdaW3L1, lamdaWL21, DW3, DsW, lr, max_iter);

    YW1 = Y * W1;
    YW2 = Y * W2;
    YW3 = Y * W3;

    W1 = W1./norm(W1, 2);
    W2 = W2./norm(W2, 2);
    W3 = W3./norm(W3, 2);

    Uall = [U1, U2, U3]; % snp
    Vall = [V1, V2, V3]; % protein
    Wall = [W1, W2, W3]; % imaging

    end
end

function D = updateD(W, group)

    [n_features, n_tasks] = size(W);
        for i = 1 : n_features
            d(i) = sqrt(sum(W(i, :) .^ 2) + eps);
        end
     D = 0.5 ./ d;

end


function [D] = updateD_FGL21(u1, u2, u3)

    ulen = length(u1);
    for i = 1 : ulen
            if i == 1
                d(i) = sqrt(u1(i).^2+u2(i).^2+u3(i).^2+u1(i+1).^2+u2(i+1).^2+u3(i+1).^2+eps);
                d(i) = 0.5 ./ d(i);
    
            elseif i==ulen
                d(i) = sqrt(u1(i-1).^2+u2(i-1).^2+u3(i-1).^2+u1(i).^2+u2(i).^2+u3(i).^2+eps);
                d(i) = 0.5 ./ d(i);
    
            else
                d(i) = 0.5./(sqrt(u1(i-1).^2+u2(i-1).^2+u3(i-1).^2+u1(i).^2+u2(i).^2+u3(i).^2+eps))+0.5./(sqrt(u1(i).^2+u2(i).^2+u3(i).^2+u1(i+1).^2+u2(i+1).^2+u3(i+1).^2+eps));
            end
        D = d;
    end

end



function W3 = LR(Y, Z2, lambdaW3L1, lambdaWL21, DW3, DsW, lr, max_iter)

    [n_samples, n_features] = size(Y);
    
    n_classes = max(Z2);
    
    Z2_onehot = full(sparse(1 : n_samples, Z2, 1, n_samples, n_classes));
    
    W3 = zeros(n_features, n_classes);
    
    for iter = 1 : max_iter

        logits = Y * W3;
        prob = softmax(logits')';
        grad = (Y') * (prob - Z2_onehot) / n_samples;
        grad = grad + lambdaW3L1 * (DW3 * W3);
        grad = grad + lambdaWL21 * (DsW * W3);
        W3 = W3 - lr * grad;

    end
end

% ----------- Softmax 函数 -----------
function s = softmax(x)

    x = x - max(x, [], 1);  % 稳定性处理
    ex = exp(x);
    s = ex ./ sum(ex, 1);

end


function z = vectorization(Q)
 
    p = size(Q, 1);
    d = size(Q, 2);
    z = zeros(p * d, 1);
    z(1:end) = reshape(Q, p * d, 1);

end



function [Q] = de_vectorization(z, p, d)

    Q = reshape(z(1 : end, 1), p, d);

end


function x = soft_Threshold(x, lambda)

    p = length(x);
    for i = 1 : p
        if x(i) > lambda
        x(i) = x(i) - lambda;
        elseif x(i) < -lambda
        x(i) = x(i) + lambda;
        else
        x(i) = 0;
        end
    end

end



function [D, pen_Value] = updateDs(v1, v2, v3)
vlen = length(v1);
    for i = 1 : vlen
        d(i) = sqrt(v1(i).^2+v2(i).^2+v3(i).^2+eps);
    end
D = 0.5 ./ d;
pen_Value = sum(d);
end