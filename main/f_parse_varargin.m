function parsed_inputs = f_parse_varargin( inputs, varargin )
% Parses varargin along with a list of inputs, and
% fills in a parsed_inputs structure.
%
% authors: bohan
%
% "inputs" is a cell array of name-value pairs that holds the
% names and default values of inputs. Inputs with a default
% value of 'none' will be required to have values specified by
% the user, meaning that those inputs must be specified by
% varargin. Varargin will of course override default values.
%
% Inputs:
%   inputs
%       Type: cell array
%       Desc: Cell array of name-value pairs of required
%             inputs.
%             The array looks like { 'input1_name',
%             <default_val1>, 'input2_name', <default_val2>, ...
%             'inputn_name', <default_valn> }
%   varargin
%       Name-value pairs of inputs
% 
% Outputs:
%   parsed_inputs
%       Type: struct
%       Desc: Struct with fields named after inputs, with
%             values set by varargin or by the defaults
%
% Example:
%       inputs          = { 'notes', 'none' };
%       varargin        = { 'notes', 'hello' };
%       parsed_inputs   = f_parse_varargin( inputs, varargin )
%       parsed_inputs = 
%           struct with fields:
%             notes: 'hello'

% parse the inputs
p = struct();
for ii = 1:2:length(varargin)-1
    p.(varargin{ii}) = varargin{ii+1}; % = setfield( p, varargin{ii}, varargin{ii+1} );
end

% check existence of required inputs
for ii = 1:2:( length(inputs)-1 )

    if ~isfield( p, inputs{ii} )
        % this input was not found

        if ischar( inputs{ii+1} ) && strcmp( inputs{ii+1}, 'none' )
            % this input has no default, throw error
            error( 'Input ''%s'' was not set and requires a user-defined value. Try again.', inputs{ii} );
        else
            % this input has a default, set the default
            fprintf( 'Input ''%s'' was not set, setting to default value ''%s''\n', inputs{ii}, num2str(inputs{ii+1}) );
            p.(inputs{ii}) = inputs{ii+1};
        end

    end     % end if ~isfield

end     % end for loop

parsed_inputs = p;

end     % end parse varrgin

