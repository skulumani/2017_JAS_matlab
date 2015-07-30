% 17 september 2014 - using for crtbp problem as well

function [value,isterminal,direction]=events_xcross(t,x)

if abs(t) > 1.e-1, % wait a short time before checking to avoid initial crossing
      isterminal = 1;
    else
      isterminal = 0;
end
    
value=x(2);  % stop at x axis cross (y is 0 or x(2) = 0)                       % define events function

direction=-1;
end
