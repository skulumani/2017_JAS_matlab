function statedot=pcrtbp_ode(t,state,mu)

%CRTBP ODE Function for Circular Restricted Three Body Problem
%
% foo(x,y,z) example usage
%
% Multi-line paragraphs of descriptive text go here. It's fine for them to
% span lines. It's treated as preformatted text; help() and doc() will not
% re-wrap lines. In the editor, you can highlight paragraphs, right-click,
% and choose "Wrap selected comments" to re-flow the text.
%
%   Purpose: 
%       - Describe purpose of code
%
%   Inputs: 
%       - List all inputs into function
%
%   Outputs: 
%       - List/describe outputs of function
%
%   Dependencies: 
%       - list dependent m files required for function
%
%   Author: 
%       - Shankar Kulumani 31 August 2014
%           - list revisions
%
%   References
%       - list reference materials  
%
% More detailed help is in the <a href="matlab: help foo>extended_help">extended help</a>.
%
% Examples:
% foo(1,2,3)



% seperate out the states

% 
% x = state(1);
% y = state(2);
% % z = state(3);
% 
% xd = state(3);
% yd = state(4);
% % zd = state(6);
% 
% % constants/parameters
% 
% d1 = sqrt((x+mu)^2 + y^2 );
% d2 = sqrt( (x+mu-1)^2+ y^2 );
% 
% ux = x - ((1-mu)*(x+mu))/(d1^3) - mu*(x+mu-1)/(d2^3);
% uy = y - (1-mu)/d1^3 * y - mu *y/d2^3;
% % uz = -(1-mu)/d1^3 *z - mu*z/d2^3;
% 
% statedot = [xd;yd;...
%     2*yd + ux;...
%     -2*xd + uy];

statedot = pcrtbp_ode_update(state,mu,zeros(2,1));

function extended_help
%EXTENDED_HELP Some additional technical details and examples
%
% Here is where you would put additional examples, technical discussions,
% documentation on obscure features and options, and so on.

error('This is a placeholder function just for helptext');