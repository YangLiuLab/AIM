%% compile_mex.m - Compile all AIM MEX files
%
% macOS: Requires Homebrew GCC: brew install gcc

sources = {'IntersectionMax_core.cpp', 'IntersectionMaxZ_core.cpp', ...
           'PointIntersect2D.cpp', 'PointIntersect1D.cpp'};

%% Check source files
missing = {};
for k = 1:length(sources)
    if ~isfile(sources{k}), missing{end+1} = sources{k}; end %#ok<SAGROW>
end
if ~isempty(missing)
    fprintf('Missing:\n');
    for k = 1:length(missing), fprintf('  %s\n', missing{k}); end
    error('Place all .cpp files in: %s', pwd);
end

%% Back up old MEX
old = dir(['*.' mexext]);
for k = 1:length(old)
    name = old(k).name;
    if ~endsWith(name, '.bak')
        fprintf('Backing up: %s\n', name);
        movefile(name, [name '.bak']);
    end
end

%% Detect interleaved complex (R2018a+)
v = ver('MATLAB');
year = str2double(v.Version(1:find(v.Version=='.',1)-1));
if year >= 9  % R2018a = 9.4
    complex_flag = '-DMX_HAS_INTERLEAVED_COMPLEX';
else
    complex_flag = '';
end

%% Compile
ext = mexext;

if ispc
    fprintf('Platform: Windows (MSVC + OpenMP)\n\n');
    for k = 1:length(sources)
        fprintf('Compiling %s ...', sources{k});
        mex(sources{k}, 'COMPFLAGS="$COMPFLAGS /openmp"');
        fprintf(' OK\n');
    end

elseif ismac
    gcc_found = '';
    for v_gcc = 15:-1:11
        for prefix = {'/opt/homebrew/bin', '/usr/local/bin'}
            p = sprintf('%s/g++-%d', prefix{1}, v_gcc);
            if isfile(p), gcc_found = p; break; end
        end
        if ~isempty(gcc_found), break; end
    end

    mex_inc = [matlabroot '/extern/include'];
    mex_lib = [matlabroot '/bin/' computer('arch')];

    if contains(computer('arch'), 'maca64')
        arch_flag = '-arch arm64';
    else
        arch_flag = '-arch x86_64';
    end

    if ~isempty(gcc_found)
        fprintf('Platform: macOS with %s + OpenMP\n\n', gcc_found);
        for k = 1:length(sources)
            src = sources{k};
            obj = strrep(src, '.cpp', '.o');
            out = strrep(src, '.cpp', ['.' ext]);
            fprintf('Compiling %s ...', src);

            cmd = sprintf('%s -c %s -o %s -I"%s" -fPIC -fopenmp -O3 -DMATLAB_MEX_FILE %s %s 2>&1', ...
                gcc_found, src, obj, mex_inc, complex_flag, arch_flag);
            [s, r] = system(cmd);
            if s ~= 0, error('\nCompile failed:\n%s\n%s', cmd, r); end

            cmd = sprintf('%s -shared %s -o %s -fopenmp -L"%s" -lmx -lmex -lmat %s 2>&1', ...
                gcc_found, obj, out, mex_lib, arch_flag);
            [s, r] = system(cmd);
            if s ~= 0, error('\nLink failed:\n%s\n%s', cmd, r); end

            delete(obj);
            fprintf(' OK\n');
        end
    else
        fprintf('No Homebrew GCC. Compiling WITHOUT OpenMP.\n');
        for k = 1:length(sources)
            fprintf('Compiling %s ...', sources{k});
            mex(sources{k});
            fprintf(' OK\n');
        end
    end

else
    fprintf('Platform: Linux (GCC + OpenMP)\n\n');
    for k = 1:length(sources)
        fprintf('Compiling %s ...', sources{k});
        mex(sources{k}, 'CXXFLAGS=$CXXFLAGS -fopenmp', 'LDFLAGS=$LDFLAGS -fopenmp');
        fprintf(' OK\n');
    end
end

%% Verify
fprintf('\n=== Results ===\n');
compiled = dir(['*.' ext]);
count = 0;
for k = 1:length(compiled)
    if ~endsWith(compiled(k).name, '.bak')
        fprintf('  OK: %s\n', compiled(k).name);
        count = count + 1;
    end
end
if count == length(sources)
    fprintf('\nAll %d MEX files compiled. Run "clear mex" then your script.\n', count);
else
    fprintf('\nWARNING: Expected %d, got %d.\n', length(sources), count);
end
