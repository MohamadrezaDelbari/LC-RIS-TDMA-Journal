function omega = fminfinder(x0, x1, lambda, phi,tmax,xmax,xmin)
    %%%% Simple version of function used in conference version
    %f = @(x, x0, x1, lambda, phi) -(sign(x-x0)-2).*(x - x0).^2-(sign(x1-x)-2).*(x - x1).^2 - lambda .* cos(x - phi);
    %%%% Function without any change
    f = @(x, x0, x1, lambda, phi) 9*(sign(x-x0)+1)/2*log((xmax-x0)/(xmax-x))-29*(sign(x-x0)-1)/2*log((xmin-x0)/(xmin-x)) - lambda .* cos(x - phi)+9*(sign(x1-x)+1)/2*log((xmax-x)/(xmax-x1))-29*(sign(x1-x)-1)/2*log((xmin-x)/(xmin-x1));




    % Initialize omega
    %omega = phi;

    % Loop over each element using arrayfun (to handle vector inputs)
    omega = arrayfun(@(x0, x1, l, p ,tmax, xmin , xmax) findMin(f, x0, x1, l, p ,tmax, xmin , xmax), x0,x1, lambda*ones(length(x0),1), phi,tmax*ones(length(x0),1), xmin*ones(length(x0),1), xmax*ones(length(x0),1));
end

function om = findMin(f, x0, x1, lambda, phi,tmax,xmin,xmax)
om=phi;
    if tmax<f(phi,x0,x1,0,phi)
    a = fminbnd(@(x) f(x, x0,x1, lambda, phi), xmin,xmax);
        if f(a, x0,x1, lambda, phi) < f(om, x0,x1, lambda, phi)
            om=a;
        end
    end
end












% Previous function. It worked.

% function om = findMin(f, delta, x0, lambda, phi)
%     if phi - delta < 0
%         % If condition to adjust the interval wrapping around 2*pi
%         a1 = fminbnd(@(x) f(x, x0, lambda, phi), phi - delta + 2*pi,2*pi);
%         a2 = fminbnd(@(x) f(x, x0, lambda, phi), 0,phi+delta);
%         % Select minimum value
%         if f(a1, x0, lambda, phi) < f(a2, x0, lambda, phi)
%             om = a1;
%         else
%             om = a2;
%         end
%     elseif phi+delta>2*pi
%         % If condition to adjust the interval wrapping around 2*pi
%         a1 = fminbnd(@(x) f(x, x0, lambda, phi), phi - delta ,2*pi);
%         a2 = fminbnd(@(x) f(x, x0, lambda, phi), 0,phi+delta-2*pi);
%         % Select minimum value
%         if f(a1, x0, lambda, phi) < f(a2, x0, lambda, phi)
%             om = a1;
%         else
%             om = a2;
%         end
%     else
%         % Normal case
%         om = fminbnd(@(x) f(x, x0, lambda, phi), phi - delta, phi + delta);
%     end
% end
