function EXPR = str2sym(S)
%STR2SYM Evaluate string representing symbolic expression.
%   str2sym(S) evaluates S as a symbolic expression, where S can be
%   a string, character vector, or cell array of character vectors.
%
%   str2sym(S) behavior:
%     Only executes functions on path when converting the string.
%     Does not consider variables and values from workspace.
%     Does not add variables in the input to the workspace.
%     Does not support programming statements like IF.
%     Treats assignments ('=') as equations ('==').
%
%   For numeric strings, use SYM for exact symbolic numbers and
%   VPA for floating-point numbers.
%
%   See also SUBS, SYM, SYMFUN, SYMS, SYMVAR, VPA.
%
%   Examples:
%
%      >> str2sym("b^2 + pi/2 + sin(x) = x^3 - 2")
%      ans =
%      b^2 + pi/2 + sin(x) == x^3 - 2
%
%      >> str2sym('int(x^3,x)')
%      ans =
%      x^4/4
%
%      >> str2sym('y(0) = 42')
%      ans =
%      y(0) == 42
%
%      >> str2sym({'42', 'a', 'y(0)'})
%      ans =
%      [ 42, a, y(0)]
%
%      >> str2sym('[42, a, y(0)]')
%      ans =
%      [ 42, a, y(0)]
%
%   Copyright 2017 The MathWorks, Inc.

narginchk(1,1);
if isa(S,'char') || iscellstr(S)
    S = string(S);
elseif ~isstring(S)
    error(message("symbolic:str2sym:InvalidFirstArgument"));
end

% Convert newlines and tabs to spaces (newline are not supported by evalin)
if any(contains(S(:),newline))
    warning(message("symbolic:str2sym:NewlinesIgnored"));
end
S = regexprep(S, "\s", " ");

% Treat assigmnments as equations by replacing "=" by "==".
% --- But do not replace "=" in "<=", ">=", "~=", and "==".
S = regexprep(S, "(?<![<>~=])=(?![<>~=])", "=="); 

% Keep size of input for output
EXPR = cell(size(S));

% Iterate linearly over string array and convert to SYM.
for idx = 1:numel(S)
    % Behave like sym function on missing.
    if ismissing(S(idx))
        EXPR{idx} = sym(string(missing));
        continue;
    end

    % Empty string results in empty sym.
    if replace(S(idx)," ","") == ""
        EXPR{idx} = sym([]);
        continue;
    end

    if ~isempty(regexp(S(idx), "^\s*%{\s*$", "once")) || ~isempty(regexp(S(idx), "^\s*%}\s*$", "once"))
        % Note that newlines were replaced by blanks before. The 
        % only valid occurance of "%{" and "%}" to start and end 
        % a block comment is in a new line with nothing else but
        % spaces. Otherwise it is handled as single-line comment.
        error(message("symbolic:str2sym:CommentsNotSupported"));
    end
    
    % Parse string and identify objects that must be treated as SYMs.
    mt = mtree(S(idx), "-com");

    % Handle syntax errors in input
    if ~isnull(mtfind(mt, "Kind", "ERR"))
        error(message("symbolic:str2sym:UnableToConvert", string(mt)));
    end

    % Handle left curly bracket node in input
    if ~isnull(mtfind(mt, "Kind", "LC")) 
        error(message("symbolic:str2sym:CellArrayNotSupported"));
    end

    % Collect objects that need special treatment
    objects = [];

    % Comments
    comments = mtfind(mt, "Kind", "COMMENT");
    if ~isnull(comments)
        comments = struct( ...
            'type',    "comment", ...
            'name',    "", ...
            'start',   num2cell(position(comments))', ...
            'end',     num2cell(endposition(comments))', ...
            'fnoargs', false, ...
            'dcall',   false, ...
            'isi',     false ...
        );
        objects = [objects comments]; %#ok<AGROW>
    end

    % Numbers
    numbers = mtfind(mt, "Kind", {"INT", "DOUBLE"});
    if ~isnull(numbers)
        numbers = struct( ...
            'type',    "number", ...
            'name',    "", ...
            'start',   num2cell(position(numbers))', ...
            'end',     num2cell(endposition(numbers))', ...
            'fnoargs', false, ...
            'dcall',   false, ...
            'isi',     false ...
        );
        objects = [objects numbers]; %#ok<AGROW>
    end

    % Functions calls, variables and special values
    calls = mtfind(mt, "Kind", {"CALL","DCALL"});
    if ~isnull(calls)
        funvars = struct( ...
            'type',    "function", ...
            'name',    strings(calls.Left), ...
            'start',   num2cell(position(calls.Left))', ...
            'end',     num2cell(endposition(calls.Left))', ...
            'fnoargs', false, ...
            'dcall',   false, ...
            'isi',     false ...
        );

        % Depending on type of call, name might need special treatment.
        wrapName = true(1,count(calls));
        cidx = indices(calls);
        for i = 1:count(calls)
            % Note that DCALL nodes always have arguments.
            % The arguments for a DCALL call are always STRINGs. 
            funvars(i).dcall = iskind(select(calls,cidx(i)), "DCALL"); 
            if ~funvars(i).dcall
                funvars(i).fnoargs = ~isempty(regexp(S(idx),"^.{" + funvars(i).end + "}\s*\(\s*\)", "once"));
                if isnull(Right(select(calls,cidx(i)))) && ~funvars(i).fnoargs
                    % Call has no arguments (= right branch of tree)  and
                    % found no name() so it's a variable or special value
                    name = funvars(i).name;
                    if name == "i" || name == "j"
                        funvars(i).type = "special";
                        funvars(i).isi = true;
                    elseif name == "pi" || name == "inf" || name == "Inf" || name == "nan" || name == "NaN" || name == "missing"
                        funvars(i).type = "special";
                        funvars(i).isi = false;
                    else
                        funvars(i).type = "variable";
                    end
                    continue;
                end
            end
            % We know: It's a function call with at least one argument.
            % Respect functions on search path: do not add sym wrapper.
            wh = string(evalin("caller", "which('" + funvars(i).name + "')"));
            if ~isempty(wh) && ~strcmp(wh, "variable")
                % Note: 'which' is *not* case-sensitive so we need to double-check!
                % e.g.: S:\Bsymbolic\matlab\toolbox\symbolic\symbolic\@sym\sym.m    % sym method
                wh = regexprep(wh, "\s*%.*$", ""); % for robustness
                % e.g.: built-in (S:\Bsymbolic\matlab\toolbox\matlab\datatypes\double)
                wh = regexprep(wh, "\)$", "");
                % e.g.: S:\Bsymbolic\matlab\toolbox\symbolic\symbolic\str2sym.m
                wh = regexprep(wh, "\.m$", "");
                % e.g.: plus is a built-in method           % string method
                if startsWith(wh, funvars(i).name) || endsWith(wh, funvars(i).name)
                    wrapName(i) = false;
                end
            end
        end
        objects = [objects funvars(wrapName)]; %#ok<AGROW>
    end

    % Manipulate given string: take action on each identified object.
    T = S(idx);
    for i = 1:numel(objects)
        switch objects(i).type
            case "comment"
                % Remove comment
                [T, objects] = removeStringBetween(T, objects(i).start, objects(i).end, objects);
            case {"number", "variable"}
                % Rewrite <o> to sym("<o>")
                [T, objects] = insertStringBefore(T, objects(i).start, "sym(""", objects);
                [T, objects] = insertStringAfter (T, objects(i).end,   """)",    objects);
            case "special"
                % Rewrite <o> to sym(<o>)
                if objects(i).isi
                    [T, objects] = insertStringBefore(T, objects(i).start, "1", objects);
                    objects(i).start = objects(i).start - 1;
                end
                [T, objects] = insertStringBefore(T, objects(i).start, "sym(", objects);
                [T, objects] = insertStringAfter (T, objects(i).end,   ")",    objects);
            case "function"
                % Rewrite 'f(<any>)' to 'feval(symengine, sym("f"), <any>)'
                if objects(i).dcall
                    continue;
                end
                pos = objects(i).end;
                while(extractBetween(T,pos,pos) ~= "(")
                    pos = pos + 1;
                end
                T = replaceBetween(T, pos, pos, ",");
                if objects(i).fnoargs
                    % Rewrite 'f()' to 'feval(symengine, sym("f"), "null()")'
                    [T, objects] = insertStringAfter(T, pos,   """null()""",                     objects);
                end
                [T, objects] = insertStringBefore(T, objects(i).start, "feval(symengine,sym(""", objects);
                [T, objects] = insertStringAfter (T, objects(i).end,   """)",                    objects);
        end
    end

    if replace(T," ","") == ""
        EXPR{idx} = sym([]);
        continue;
    end

    % Evaluate T in 'caller' workspace to create symbolic object.
    try
        EXPR{idx} = evalin("caller", T);
    catch Exception
        if any(contains(["MATLAB:TooManyOutputs", "MATLAB:maxlhs"], Exception.identifier))
            % str2sym should always return a sym value. 
            % If there is no return value, e.g. for str2sym('disp(42)'),
            % then return an empty sym.
            EXPR{idx} = sym([]);
        else
            error(message("symbolic:str2sym:UnableToConvert", Exception.message));
        end
    end
end

try
    EXPR = cell2sym(EXPR);
catch Exception
    error(message("symbolic:str2sym:UnableToConvert", Exception.message));
end
end

% Helper: removeStringBetween
function [text, objs] = removeStringBetween(text, left, right, objs)
    text = replaceBetween(text, left, right, "");
    objs = adaptStringPositions(objs, left, left-right-1);
end

% Helper: insertStringBefore
function [text, objs] = insertStringBefore(text, pos, textBefore, objs)
    text = insertBefore(text, pos, textBefore);
    objs = adaptStringPositions(objs, pos, strlength(textBefore));
end

% Helper: insertStringAfter
function [text, objs] = insertStringAfter(text, pos, textAfter, objs)
    text = insertAfter(text, pos, textAfter);
    objs = adaptStringPositions(objs, pos+1, strlength(textAfter));
end

% Helper: adaptStringPositions
function objs = adaptStringPositions(objs, pos, len)
for i=1:numel(objs)
    if pos <= objs(i).start
        objs(i).start = objs(i).start + len;
    end
    if pos <= objs(i).end
        objs(i).end = objs(i).end + len;
    end
end
end
