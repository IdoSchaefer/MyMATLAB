clf
syms F V r real;
funchoice = input('Choose: a potential/force function(V/F):', 's');
if funchoice == 'V'
        getV = input('Enter a potential function in the terms of r: ', 's');
        V = eval(getV);
        F = -diff(V);
    else
        getF = input('Enter a force function in the terms of r: ', 's');
        F = eval(getF);
end
N = input('Enter the number of the bodies: ');
mdqsn = input('Do you want to insert mass defects?(y, n)?', 's');
if mdqsn ~= 'y' && mdqsn ~= 'Y'
    m = input('Enter the mass:');
else
    sm = input('Enter the standard mass: ');
    m = ones(1, N)*sm;
    md = input('Enter the first kind of mass defect:');
    nmd = input('Enter the number of the body with the mass defect (0 to continue): ');
    while md>0 && nmd<=N && nmd>=1            
        while nmd<=N && nmd>=1
            m(nmd) = md;
            nmd = input('Enter the number of the body with the mass defect (0 to continue): ');
        end
        md = input('Enter another kind of mass defect (0 to continue): '); 
    end
end
tfinal = input('Enter the time interval:');
boundary = input('Choose boundary conditions: open, fixed, or cyclic (o, f, c):', 's');
r0qsn = input('Choose: set the standard distance as the equilibrium distance of one spring, \n or as your choice(e, c)?', 's');
if r0qsn == 'c'
    r0 = input('Enter the standard distance:');
else
    solveF = subs(solve(F));
    r0 = min(solveF);
end
initx = 0:r0:(N-1)*r0;
initv = zeros(1, N);
initchoice = input('Choose: initialize values as: your choice, a sinus, a gaussian, a rectangular pulse, \n another r function or the displacement q function(c, s, g, p, r, q)?', 's');
if nnz(initchoice == ['s', 'g', 'p', 'r', 'q'])==0
        num = input('Enter the number of the body:');
        while num>=1 && num<=N
            initx(num) = input('Enter the initial x:');
            initv(num) = input('Enter the initial v:');
            num = input('Enter the number of the body:');
        end
    elseif initchoice == 'g' || initchoice == 'p' || initchoice == 'r'
        if initchoice == 'g'            
            A = input('Enter the gaussian''s height: ');
            sigma = input('Enter the standard deviation, with the spring number as the variable:');
            n0 = input('Enter the center of the gaussian, with the spring number as the variable:');
            ri = 1:(N-1);
            rvalue = r0 - A*exp(-(ri - n0).^2/(2*sigma^2));
        elseif initchoice == 'p'
            A = input('Enter the rectangular''s height: ');
            B = input('Enter the rectangular''s width, with the spring number as the variable: ');
            n0 = input('Enter the center of the rectangular, with the spring number as the variable:');
            rvalue = ones(1, N-1)*r0;
            rrec = r0 - A;
            rvalue(ceil(n0 - (B-1)/2):floor(n0 + (B-1)/2)) = rrec;
        else
            f = input('Enter a function of ri (the index of the spring): r = ', 's');
            ri = 1:(N-1);
            rvalue = r0 + eval(f);
        end
        rfit = ((N-1)*r0 - sum(rvalue))/(N-1);
        rvalue = rvalue + rfit;
        for xi = 2:N
            initx(xi) = initx(xi-1) + rvalue(xi - 1);
        end
    elseif initchoice == 'q'
        f = input('Enter a function of n (the index of the body): q = ', 's');
        n = 1:N;
        q = eval(f);
        initx = initx + q;
    elseif initchoice == 's'
        A = input('Enter the amplitude:');
        if boundary == 'c'
            Ncycles = input('How many cycles?');
            q = A*sin(2*pi*Ncycles/(N*r0)*(0:r0:(N-1)*r0));
        else
            mode = input('Which mode?');
            q = A*sin(pi*mode/((N-1)*r0)*(0:r0:(N-1)*r0));
        end
        initx = initx + q;
end