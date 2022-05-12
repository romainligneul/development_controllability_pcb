function [shuffled_x] = shuffles(x)
%SHUFFLES
    shuffled_x=x(randperm(length(x)));
end

