function [stars] = getStars(p)
if p < 0.01
    stars = '***';
elseif p < 0.05
    stars = '**';
elseif p < 0.10
    stars = '*';
else
    stars = '';
end