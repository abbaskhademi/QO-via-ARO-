n = 99;     % size n
dens = 0.5; % density
dvert = 2;
myseed = 1;

for igen = 1:5
    r = (4 * myseed + 1) / 16384 / 16384;
    myseed = myseed + 1;

    Fc = zeros(n+1, n+1);
    Fl = zeros(n+1, n+1);
    F = zeros(n+1, n+1);

    for i = 1:n
        for j = (i+1):(n+1)
            [r, num] = myrandom(r, 0, 1);
            if num < dens
                [r, num] = myrandom(r, 0, 10);
                Fc(i, j) = num;
                Fc(j, i) = num;
            else
                [r, num] = myrandom(r, -10, 0);
                Fc(i, j) = num;
                Fc(j, i) = num;
            end
        end
    end

    for i = 1:(n+1)
        [r, num] = myrandom(r, 0, dvert);
        Fl(i, i) = num;
    end

    for i = 1:n
        for j = (i+1):(n+1)
            Fl(i, j) = 0.5 * (Fl(i, i) + Fl(j, j));
            Fl(j, i) = Fl(i, j);
        end
    end

    for i = 1:(n+1)
        for j = i:(n+1)
            F(i, j) = Fl(i, j) - Fc(i, j);
            F(j, i) = F(i, j);
        end
    end
    Q = F;
    filename = sprintf('z7Problem_%dx%d(%.1f)_%d.mat', n+1, n+1, dens, igen);
    save(filename, 'Q');
end

function [r, num] = myrandom(r, a, b)
    r = mod(r * 41475557, 1);
    num = a + r * (b - a);
end
