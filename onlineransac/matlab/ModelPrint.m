function ModelPrint(m)
% Print in console the model

    fprintf('Model[%s] ',m.type);
    if ~m.defined
        fprintf('UNDEF');
    else
        fprintf('[');
        if strcmp(m.type,'LL3')
            fprintf('a=%.9f, b=%.9f, c=%.9f',m.coeffs.a,m.coeffs.b,m.coeffs.c);
        elseif strcmp(m.type,'LN3')
            fprintf('g=%.9f, m=%.9f, s=%.9f',m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);
        elseif strcmp(m.type,'EXP2')
            fprintf('a=%.9f, b=%.9f',m.coeffs.alpha,m.coeffs.beta);
        else
            error('Unknown model type');
        end
        fprintf(']');
    end
    fprintf('\n');

end